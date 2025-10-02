# message("Loading counts ...")
# stopifnot(file.exists(mtx_path), file.exists(barcodes_path), file.exists(features_path))

# counts <- Matrix::readMM(mtx_path)
# barcodes <- readr::read_tsv(barcodes_path, col_names = FALSE, show_col_types = FALSE)[[1]]
# features <- readr::read_tsv(features_path, col_names = FALSE, show_col_types = FALSE)

# barcodes   <- as.character(barcodes)
# gene_names <- make.unique(as.character(features[[2]]))
# rownames(counts) <- gene_names
# colnames(counts) <- barcodes

# message("Loading positions ...")
# pos <- readr::read_csv(positions_csv, col_names = FALSE, show_col_types = FALSE)
# stopifnot(ncol(pos) >= 6)
# colnames(pos)[1:6] <- c("barcode","in_tissue","row","col","pxl_col_fullres","pxl_row_fullres")
# pos <- subset(pos, in_tissue == 1 & barcode %in% colnames(counts))

# message("Loading scalefactor ...")
# sf <- jsonlite::fromJSON(scalef_json)
# hires_scale <- as.numeric(sf$tissue_hires_scalef)

# message("Loading image ...")
# if (grepl("\\.png$", img_path, ignore.case = TRUE)) {
#   img <- png::readPNG(img_path)
# } else {
#   img <- jpeg::readJPEG(img_path)
# }

options(shiny.maxRequestSize = 1024*1024^2)  # 1 GB

# 读 counts + barcodes + features
load_counts <- function(mtx_path, barcodes_path, features_path) {
  stopifnot(file.exists(mtx_path), file.exists(barcodes_path), file.exists(features_path))

  counts   <- Matrix::readMM(mtx_path)
  barcodes <- readr::read_tsv(barcodes_path, col_names = FALSE, show_col_types = FALSE)[[1]]
  features <- readr::read_tsv(features_path, col_names = FALSE, show_col_types = FALSE)

  gene_names <- make.unique(as.character(features[[2]]))
  rownames(counts) <- gene_names
  colnames(counts) <- as.character(barcodes)

  list(counts = counts, gene_names = gene_names, barcodes = as.character(barcodes))
}

# 读 positions（tissue_positions_list.csv）
load_positions <- function(positions_csv) {
  stopifnot(file.exists(positions_csv))
  pos <- readr::read_csv(positions_csv, col_names = FALSE, show_col_types = FALSE)
  stopifnot(ncol(pos) >= 6)
  colnames(pos)[1:6] <- c("barcode","in_tissue","row","col","pxl_col_fullres","pxl_row_fullres")
  pos
}

# 读 scalefactors_json.json
load_scalef <- function(scalef_json) {
  stopifnot(file.exists(scalef_json))
  sf <- jsonlite::fromJSON(scalef_json)
  hires_scale <- as.numeric(sf$tissue_hires_scalef)
  stopifnot(is.finite(hires_scale), hires_scale > 0)
  hires_scale
}

# 读 HE 图像（png/jpg）
load_image <- function(path) {
  if (grepl("\\.png$", path, ignore.case = TRUE)) {
    img_arr <- png::readPNG(path)   # [h, w, c]
  } else if (grepl("\\.jpe?g$", path, ignore.case = TRUE)) {
    img_arr <- jpeg::readJPEG(path) # [h, w, c]
  } else {
    stop("Only PNG/JPG supported.")
  }
  img_h <- dim(img_arr)[1]
  img_w <- dim(img_arr)[2]
  img_rs <- as.raster(img_arr)      # <- 关键：转成 raster
  list(img = img_rs, img_w = img_w, img_h = img_h)
}


# fullres → hires 坐标映射
to_hires_coords <- function(pos, hires_scale) {
  pos$x_hires <- pos$pxl_col_fullres * hires_scale
  pos$y_hires <- pos$pxl_row_fullres * hires_scale
  pos
}



# img_w <- dim(img)[2]
# img_h <- dim(img)[1]

# # fullres → hires 坐标
# pos$x_hires <- pos$pxl_col_fullres * hires_scale
# pos$y_hires <- pos$pxl_row_fullres * hires_scale

# gene_choices <- head(gene_names, MAX_GENES)

# # Debug check
# cat(">>> types after loading:\n")
# print(str(gene_names))
# print(str(barcodes))
# print(str(rownames(counts)))
# print(str(colnames(counts)))

