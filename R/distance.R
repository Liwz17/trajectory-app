# 依赖：RANN （或 FNN）
# install.packages("RANN")

# 输入：
# - pos: data.frame，必须含 x_hires, y_hires（像素坐标）和可选 barcode
# - selected_idx: 整数向量（1..nrow(pos)），表示被圈选的行号（tumor）
# 输出：
# - data.frame(pos 原列 + category1 + distance1)
compute_signed_distance_to_selection <- function(pos, selected_idx) {
  n <- nrow(pos)
  if (n == 0) return(transform(pos, category1 = NA_character_, distance1 = NA_real_))

  tumor_idx     <- sort(unique(as.integer(selected_idx)))
  tumor_mask    <- rep(FALSE, n); tumor_mask[tumor_idx] <- TRUE
  non_tumor_idx <- which(!tumor_mask)

  out <- pos
  out$category1 <- ifelse(tumor_mask, "tumor", "non-tumor")
  out$distance1 <- NA_real_

  # 空集保护
  if (length(tumor_idx) == 0 && length(non_tumor_idx) == 0) return(out)
  coords <- cbind(out$x_hires, out$y_hires)

  # 非肿瘤 → 最近肿瘤（正号）
  if (length(non_tumor_idx) > 0 && length(tumor_idx) > 0) {
    nn_nt_to_t <- RANN::nn2(data = coords[tumor_idx, , drop = FALSE],
                            query = coords[non_tumor_idx, , drop = FALSE], k = 1)$nn.dists[, 1]
    out$distance1[non_tumor_idx] <- nn_nt_to_t  # 正数
  }

  # 肿瘤 → 最近非肿瘤（负号）
  if (length(tumor_idx) > 0 && length(non_tumor_idx) > 0) {
    nn_t_to_nt <- RANN::nn2(data = coords[non_tumor_idx, , drop = FALSE],
                            query = coords[tumor_idx, , drop = FALSE], k = 1)$nn.dists[, 1]
    out$distance1[tumor_idx] <- -nn_t_to_nt     # 负数
  }

  out
}
