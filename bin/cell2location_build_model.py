#!/usr/bin/env python3
# isort: off
import argparse
import os
import time
import logging
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

# silence scanpy that prints a lot of warnings
import warnings

import anndata
import scvi

# Run the model - use all data for training (validation not implemented yet, train_size=1)
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes
from matplotlib import rcParams


warnings.filterwarnings("ignore")
rcParams["pdf.fonttype"] = 42

logger = logging.getLogger()


def get_args():
    # Script description
    description = """Build model"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "sc_data_path", type=str, help="path to single cell h5ad count data"
    )

    parser.add_argument("cuda_device", type=str, help="index of cuda device ID or cpu")

    parser.add_argument(
        "-o", "--out_dir", default=None, type=str, help="directory for regression model"
    )

    parser.add_argument(
        "-a",
        "--annotation_column",
        default="celltype",
        type=str,
        help="column name for covariate",
    )

    parser.add_argument(
        "-s",
        "--sample_column",
        default=None,
        type=str,
        help="column containing sample id or batch information",
    )

    parser.add_argument(
        "-t",
        "--tech_column",
        default=None,
        nargs="+",
        type=str,
        help="multiplicative technical effects, such as platform effects",
    )

    parser.add_argument(
        "-e",
        "--epochs",
        default=250,
        type=int,
        help="number of epochs to train the model",
    )

    parser.add_argument(
        "-p",
        "--posterior_sampling",
        default=1000,
        type=int,
        help="number of samples to take from the posterior distribution",
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


def build_model(
    cuda_device,
    sc_data_path,
    sample_column,
    annotation_column,
    tech_column,
    epochs=250,
    posterior_sampling=1000,
):
    ##### MAIN CODE #####

    ## scRNA reference (raw counts)
    logging.info("Reading scRNA-seq data from " + sc_data_path + "...")
    adata_scrna_raw = anndata.read_h5ad(sc_data_path)
    adata_scrna_raw.var["SYMBOL"] = adata_scrna_raw.var_names

    # Filter genes
    logging.info(
        "Before filtering: {} cells and {} genes.".format(*adata_scrna_raw.shape)
    )

    selected = filter_genes(
        adata_scrna_raw,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    adata_scrna_raw = adata_scrna_raw[:, selected].copy()
    logging.info(
        "After filtering: {} cells and {} genes.".format(*adata_scrna_raw.shape)
    )

    logging.info("Preparing anndata for the regression model...")
    scvi.data.setup_anndata(
        adata=adata_scrna_raw,
        # 10X reaction / sample / batch
        batch_key=sample_column,
        # cell type, covariate used for constructing signatures
        labels_key=annotation_column,
        # multiplicative technical effects (platform, 3' vs 5', donor effect)
        categorical_covariate_keys=tech_column,
    )
    logging.info("Build model...")
    mod = RegressionModel(adata_scrna_raw)
    mod.train(
        max_epochs=epochs,
        batch_size=2500,
        train_size=1,
        lr=0.002,
        use_gpu=cuda_device.isdigit(),
    )

    logging.info(
        "Export the estimated cell abundance (summary of the posterior distribution)..."
    )
    adata_scrna_raw = mod.export_posterior(
        adata_scrna_raw,
        sample_kwargs={
            "num_samples": posterior_sampling,
            "batch_size": 2500,
            "use_gpu": cuda_device.isdigit(),
        },
    )
    return mod, adata_scrna_raw


def main(args):
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")

    cuda_device = args.cuda_device

    assert cuda_device.isdigit() or cuda_device == "cpu", "invalid device id"
    assert os.path.exists(args.sc_data_path), "sc file not found"

    if cuda_device.isdigit():
        os.environ["CUDA_VISIBLE_DEVICES"] = cuda_device

    logging.info("Parameters\n==========")
    logging.info(
        "Training epochs: {}\nPosterior sampling: {}".format(
            args.epochs, args.posterior_sampling
        )
    )
    logging.info("==========")

    mod, adata_scrna_raw = build_model(
        cuda_device,
        sc_data_path=args.sc_data_path,
        sample_column=args.sample_column,
        annotation_column=args.annotation_column,
        tech_column=args.tech_column,
        epochs=args.epochs,
        posterior_sampling=args.posterior_sampling,
    )

    # Save model and anndata object with results
    mod.save(args.out_dir, overwrite=True)

    try:
        adata_scrna_raw.write(Path(args.out_dir, "sc.h5ad"))
    except ValueError:
        logging.info(
            "There seems to be an issue with the conversion. Renaming columns..."
        )
        os.remove(Path(args.out_dir, "sc.h5ad"))
        adata_scrna_raw.__dict__["_raw"].__dict__["_var"] = (
            adata_scrna_raw.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
        adata_scrna_raw.write(Path(args.out_dir, "sc.h5ad"))


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    logging.info(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
