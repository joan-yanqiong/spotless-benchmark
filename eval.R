# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

require(GaitiLabUtils)
# Set working directory
set_wd()

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
# Comment/uncomment depending on whether you have an internal package based on the 'R' directory created with usethis::create_package()
# devtools::load_all("./", export_all = FALSE)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Check combining outputs",
    default_output = "output",
    default_log_file = NULL,
    default_log_dir = "logs"
)
params <- parser$parse_args()
if (interactive()) {
    # Provide arguments here for local runs
    params$output_dir <- file.path(getwd(), "output_testing")
}

# Set up logging
logr <- init_logging(log_level = params$log_level, log_file = NULL)
obj_logger <- init_obj_logging(log_file = NULL)
log_info(ifelse(
    interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(params$output_dir)

# Load additional libraries

# sample_id <- "test_spatial_sample1"

# c2l_deconv <- as.matrix(
#     fread(file.path(
#         params$output_dir,
#         "cell2location",
#         paste0("cell2location__", sample_id, ".txt")
#     )),
#     rownames = 1
# )

# spatialdwls_deconv <- as.matrix(
#     fread(file.path(
#         params$output_dir,
#         "spatialdwls",
#         paste0("spatialdwls__", sample_id, ".txt")
#     )),
#     rownames = 1
# )

# rctd_deconv <- as.matrix(
#     fread(file.path(
#         params$output_dir,
#         "rctd",
#         paste0("rctd__", sample_id, ".txt")
#     )),
#     rownames = 1
# )

mat <- as.matrix(
    fread("output_testing/spatialdwls/spatialdwls__test_spatial_sample1.txt"),
    rownames = 1
)
