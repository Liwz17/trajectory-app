# 取表达向量（来自 counts 或文件的 expr 列）
make_expr_from_counts <- function(gene, counts, barcodes, log1p = TRUE) {
  if (is.null(gene) || is.null(counts) || is.null(barcodes)) return(NULL)
  if (!(gene %in% rownames(counts))) return(NULL)
  v <- as.numeric(counts[gene, barcodes, drop = TRUE])
  if (isTRUE(log1p)) v <- log1p(v)
  v
}

# 合并：距离文件 df_dist + (expr_from_counts 或 df_dist$expr)
merge_dist_expr <- function(df_dist, expr_from_counts = NULL) {
  stopifnot("barcode" %in% names(df_dist))
  df <- df_dist
  if (!is.null(expr_from_counts)) {
    df$expr <- expr_from_counts
  } else {
    if (!("expr" %in% names(df_dist))) {
      stop("距离文件没有 expr 列，且当前未加载 counts 无法计算表达量。")
    }
  }
  # 清理非法
  bad <- !is.finite(df$expr)
  if (any(bad)) df$expr[bad] <- NA_real_
  df
}

# 按距离做分箱统计
bin_summary <- function(x, y, bins = 50) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(data.frame())
  x <- x[ok]; y <- y[ok]
  bx <- cut(x, breaks = bins, include.lowest = TRUE)
  agg <- aggregate(
    cbind(y, x) ~ bx,
    FUN = function(z) c(mean = mean(z), sd = sd(z), n = length(z))
  )
  # 展开多列矩阵
  out <- data.frame(
    bin     = agg$bx,
    x_mean  = agg$x[, "mean"],
    x_sd    = agg$x[, "sd"],
    n       = agg$y[, "n"],
    y_mean  = agg$y[, "mean"],
    y_sd    = agg$y[, "sd"]
  )
  out
}
