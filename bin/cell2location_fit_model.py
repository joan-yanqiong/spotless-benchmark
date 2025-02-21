#!/usr/bin/env python3
# isort: off
import logging
import argparse
import os
import time
import warnings
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import scanpy as sc
import numpy as np

import cell2location
import scvi

from matplotlib import rcParams

rcParams["pdf.fonttype"] = 42
logger = logging.getLogger()

warnings.filterwarnings("ignore")


def get_args():
    # Script description
    description = """Build model"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("sp_data_path", type=str, help="path of spatial data")

    parser.add_argument("model_path", type=str, help="path to regression model")

    parser.add_argument("cuda_device", type=str, help="index of cuda device ID or cpu")

    parser.add_argument(
        "-o",
        "--out_dir",
        default=os.getcwd(),
        type=str,
        help="model and proportions output directory",
    )

    parser.add_argument(
        "-e",
        "--epochs",
        default=30000,
        type=int,
        help="number of epochs to fit the model",
    )

    parser.add_argument(
        "-p",
        "--posterior_sampling",
        default=1000,
        type=int,
        help="number of samples to take from the posterior distribution",
    )

    parser.add_argument(
        "-n",
        "--n_cells_per_location",
        default=8,
        type=int,
        help="estimated number of cells per spot",
    )

    parser.add_argument(
        "-d",
        "--detection_alpha",
        default=200,
        type=int,
        help="within-experiment variation in RNA detection sensitivity",
    )

    arg = parser.parse_args()
    arg.out_dir = abspath(arg.out_dir)

    if arg.out_dir is None:
        arg.out_dir = Path(os.path.dirname(arg.sc_data_path), "cell2location_results")
    else:
        arg.out_dir = arg.out_dir

    if not os.path.exists(arg.out_dir):
        os.makedirs(arg.out_dir)

    return arg


def fit_model(
    cuda_device,
    sp_data_path,
    n_cells_per_location=8,
    detection_alpha=200,
    posterior_sampling=1000,
    epochs=30000,
):
    logging.info("Reading in spatial data from " + sp_data_path + "...")
    adata = sc.read_h5ad(sp_data_path)
    adata.var["SYMBOL"] = adata.var_names

    # mitochondria-encoded (MT) genes should be removed for spatial mapping
    adata.var["mt"] = [gene.startswith("mt-") for gene in adata.var["SYMBOL"]]
    adata = adata[:, ~adata.var["mt"].values]

    adata_vis = adata.copy()
    adata_vis.raw = adata_vis

    logging.info("Reading in the model...")
    adata_scrna_raw = sc.read(args.model_path)

    # Export estimated expression in each cluster
    if "means_per_cluster_mu_fg" in adata_scrna_raw.varm.keys():
        inf_aver = adata_scrna_raw.varm["means_per_cluster_mu_fg"][
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in adata_scrna_raw.uns["mod"]["factor_names"]
            ]
        ].copy()
    else:
        inf_aver = adata_scrna_raw.var[
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in adata_scrna_raw.uns["mod"]["factor_names"]
            ]
        ].copy()
    inf_aver.columns = adata_scrna_raw.uns["mod"]["factor_names"]

    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    scvi.data.setup_anndata(adata=adata_vis)

    # Create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=n_cells_per_location,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection (using default here):
        detection_alpha=detection_alpha,
    )

    mod.train(
        max_epochs=epochs,
        # train using full data (batch_size=None)
        batch_size=None,
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=1,
        use_gpu=cuda_device.isdigit(),
    )

    # Export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={
            "num_samples": posterior_sampling,
            "batch_size": mod.adata.n_obs,
            "use_gpu": cuda_device.isdigit(),
        },
    )

    # Export proportion file, but first rename columns and divide by rowSums
    props = adata_vis.obsm["q05_cell_abundance_w_sf"]
    props = props.rename(
        columns={x: x.replace("q05cell_abundance_w_sf_", "") for x in props.columns}
    )
    props = props.div(props.sum(axis=1), axis="index")

    return mod, adata_vis, props


def main(args):
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")

    cuda_device = args.cuda_device

    assert cuda_device.isdigit() or cuda_device == "cpu", "invalid device input"

    logging.info("Parameters\n==========")
    logging.info(
        "Detection alpha: {}\nCells per location: {}".format(
            args.detection_alpha, args.n_cells_per_location
        )
    )
    logging.info(
        "Training epochs: {}\nPosterior sampling: {}".format(
            args.epochs, args.posterior_sampling
        )
    )
    logging.info("==========")

    ##### MAIN PART #####
    if cuda_device.isdigit():
        os.environ["CUDA_VISIBLE_DEVICES"] = cuda_device

    mod, adata_vis, props = fit_model(
        cuda_device=cuda_device,
        sp_data_path=args.sp_data_path,
        n_cells_per_location=args.n_cells_per_location,
        detection_alpha=args.detection_alpha,
        posterior_sampling=args.posterior_sampling,
        epochs=args.epochs,
    )

    # Save model and anndata object with results
    mod.save(args.out_dir, overwrite=True)
    adata_vis.write(Path(args.out_dir, "sp.h5ad"))
    props.to_csv(Path(args.out_dir, "proportions.tsv"), sep="\t")


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    logging.info(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
