# global.R

# ---- Dependencies ----
# install.packages(c("shiny","ggplot2","Matrix","jsonlite","readr","jpeg","png","viridisLite"))
library(shiny)
library(ggplot2)
library(Matrix)
library(jsonlite)
library(readr)
library(jpeg)
library(png)
library(viridisLite)

# # ---- 路径配置 ----
# VISIUM_OUTS <- "../data/Multi-Section-data/ST-visium/2D/outs"

# counts_dir     <- file.path(VISIUM_OUTS, "filtered_feature_bc_matrix")
# positions_csv  <- file.path(VISIUM_OUTS, "spatial", "tissue_positions_list.csv")
# scalef_json    <- file.path(VISIUM_OUTS, "spatial", "scalefactors_json.json")
# img_path_png   <- file.path(VISIUM_OUTS, "spatial", "tissue_hires_image.png")
# img_path_jpg   <- file.path(VISIUM_OUTS, "spatial", "detected_tissue_image.jpg")
# img_path       <- if (file.exists(img_path_jpg)) img_path_jpg else img_path_png

# mtx_path       <- file.path(counts_dir, "matrix.mtx.gz")
# barcodes_path  <- file.path(counts_dir, "barcodes.tsv.gz")
# features_path  <- file.path(counts_dir, "features.tsv.gz")


# ---- 参数 ----
MAX_GENES <- 500
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b