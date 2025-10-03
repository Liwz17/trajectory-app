# 依赖：RANN （或 FNN）
# install.packages("RANN")

# 输入：
# - pos: data.frame，必须含 x_hires, y_hires（像素坐标）和可选 barcode
# - selected_idx: 整数向量（1..nrow(pos)），表示被圈选的行号（tumor）
# 输出：
# - data.frame(pos 原列 + category1 + distance1)
# 依赖：RANN
# label_index_list: 命名 list，如 list(tumor = c(1,5,9), tls = c(2,3))
# 结果：对每个 label 生成一列 dist_<label>（该结构内为负，外部为正）
# 合并并去重：同一 label 下，同一 barcode 只保留一次
append_annotations <- function(current_df, barcodes, label) {
  stopifnot(is.character(barcodes), length(label) == 1)
  add <- data.frame(barcode = unique(barcodes), label = label, stringsAsFactors = FALSE)
  out <- rbind(current_df, add)
  out <- unique(out)   # 去重
  out[order(out$label, out$barcode), , drop = FALSE]
}

# 汇总统计
summarize_annotations <- function(df) {
  if (nrow(df) == 0) {
    return(data.frame(label = character(), n = integer()))
  }
  tb <- as.data.frame(table(df$label), stringsAsFactors = FALSE)
  colnames(tb) <- c("label", "n")
  tb[order(-tb$n), ]   # 按数量降序
}




compute_signed_distance_per_label <- function(pos, label_index_list) {
  stopifnot(is.data.frame(pos), all(c("x_hires","y_hires") %in% names(pos)))
  n <- nrow(pos)
  if (n == 0 || length(label_index_list) == 0) return(pos)

  coords <- cbind(pos$x_hires, pos$y_hires)
  out <- pos

  # 为每个 label 生成带符号距离
  for (lab in names(label_index_list)) {
    idx_lab <- sort(unique(as.integer(label_index_list[[lab]])))
    idx_lab <- idx_lab[idx_lab >= 1 & idx_lab <= n]
    if (length(idx_lab) == 0) {
      out[[paste0("dist_", lab)]] <- NA_real_
      next
    }
    mask_lab <- rep(FALSE, n); mask_lab[idx_lab] <- TRUE
    idx_non  <- which(!mask_lab)

    dist_vec <- rep(NA_real_, n)
    # 非该结构 → 最近该结构（正）
    if (length(idx_non) > 0) {
      d1 <- RANN::nn2(data  = coords[idx_lab, , drop = FALSE],
                      query = coords[idx_non, , drop = FALSE], k = 1)$nn.dists[, 1]
      dist_vec[idx_non] <- d1
    }
    # 该结构 → 最近非该结构（负）
    if (length(idx_lab) > 0 && length(idx_non) > 0) {
      d2 <- RANN::nn2(data  = coords[idx_non, , drop = FALSE],
                      query = coords[idx_lab, , drop = FALSE], k = 1)$nn.dists[, 1]
      dist_vec[idx_lab] <- -d2
    } else if (length(idx_non) == 0) {
      # 所有点都属于该结构：约定为 0（或 NA），这里给 0
      dist_vec[idx_lab] <- 0
    }

    # 列名尽量合法化
    colname <- paste0("dist_", make.names(lab))
    out[[colname]] <- dist_vec
  }

  out
}


setup_multi_distance <- function(output, rv, expr_vec, label_index_list) {
  multi_distance_df <- reactive({
    req(rv$pos)
    lst <- label_index_list()
    if (length(lst) == 0) return(rv$pos)
    compute_signed_distance_per_label(rv$pos, lst)
  })

  output$dl_multi_distance <- downloadHandler(
    filename = function() paste0("spots_with_signed_distance_multi_",
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      df <- multi_distance_df()
      df$expr <- expr_vec()
      dist_cols <- grep("^dist_", names(df), value = TRUE)
      base_cols <- intersect(c("barcode","x_hires","y_hires"), names(df))
      keep <- c(base_cols, dist_cols, "expr")
      utils::write.csv(df[, keep, drop = FALSE], file, row.names = FALSE)
    }
  )

  multi_distance_df
}
