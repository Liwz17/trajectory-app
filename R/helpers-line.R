make_line_ab <- function(Ax, Ay, Bx, By, label = NULL) {
  dx <- Bx - Ax; dy <- By - Ay
  L  <- sqrt(dx*dx + dy*dy)
  if (!is.finite(L) || L <= 0) return(NULL)
  list(Ax=Ax, Ay=Ay, Bx=Bx, By=By, dx=dx, dy=dy, L=L, label = label %||% "line")
}

# ===== 方向箭头（A -> B）=====
make_arrow_annot <- function(line,
                             arrowcolor = "red",
                             arrowsize  = 1.2,
                             arrowwidth = 2,
                             arrowhead  = 3) {
  list(
    x = line$Bx,  y = line$By,   # 箭头终点（B）
    ax = line$Ax, ay = line$Ay,  # 箭头起点（A）
    xref = "x", yref = "y", axref = "x", ayref = "y",
    showarrow = TRUE,
    arrowhead = arrowhead,
    arrowsize = arrowsize,
    arrowwidth = arrowwidth,
    arrowcolor = arrowcolor,
    opacity = 0.9
  )
}


# ===== 带宽走廊（±band 的半透明区域）=====
make_band_shape <- function(line, band,
                            fillcolor = "rgba(255,0,0,0.15)",
                            linecolor = "rgba(255,0,0,0.4)") {
  dx <- line$dx; dy <- line$dy; L <- line$L
  if (!is.finite(L) || L <= 0) return(NULL)

  # 法向单位向量（左法向，注意图像坐标方向如果有翻转，必要时改符号）
  nx <- -dy / L
  ny <-  dx / L

  # 四个顶点（沿法向 ±band 偏移）
  A_plus  <- c(line$Ax + band*nx, line$Ay + band*ny)
  A_minus <- c(line$Ax - band*nx, line$Ay - band*ny)
  B_plus  <- c(line$Bx + band*nx, line$By + band*ny)
  B_minus <- c(line$Bx - band*nx, line$By - band*ny)

  # Plotly path：M x y L x y ... Z
  path_str <- sprintf("M %f %f L %f %f L %f %f L %f %f Z",
                      A_plus[1],  A_plus[2],
                      B_plus[1],  B_plus[2],
                      B_minus[1], B_minus[2],
                      A_minus[1], A_minus[2])

  list(
    type = "path",
    path = path_str,
    fillcolor = fillcolor,
    line = list(color = linecolor, width = 1)
  )
}

project_points_to_line <- function(x, y, line) {
  Ax <- line$Ax; Ay <- line$Ay; dx <- line$dx; dy <- line$dy; L <- line$L
  t <- ((x - Ax) * dx + (y - Ay) * dy) / L
  d_perp <- ((x - Ax) * dy - (y - Ay) * dx) / L
  data.frame(t = t, d_perp = d_perp)
}
make_line_shape <- function(line) {
  list(type = "line", x0 = line$Ax, y0 = line$Ay, x1 = line$Bx, y1 = line$By,
        line = list(color = "red", width = 2))
}
# —— 沿线抽样：返回 t/d_perp/expr —— #
sample_along_line <- function(pos_df, expr_vec, line, band, end_tol = 0) {
  proj <- project_points_to_line(pos_df$x_hires, pos_df$y_hires, line)
  df <- cbind(pos_df[, c("barcode","x_hires","y_hires")], proj, stringsAsFactors = FALSE)

  L <- as.numeric(line$L)
  df <- df[
    is.finite(df$t) & is.finite(df$d_perp) &
    abs(df$d_perp) <= band &
    df$t >= (0 - end_tol) & df$t <= (L + end_tol),   # ✅ 关键：把 t 限制到 [0, L]
    , drop = FALSE
  ]

  df$expr <- as.numeric(expr_vec)[match(df$barcode, pos_df$barcode)]
  df <- df[order(df$t), , drop = FALSE]

  # 如果你想看相对位置（0=A, 1=B），顺带给出一个列（可选）
  df$t_norm <- df$t / L
  df
}
safe_gene_label <- function(x, alt = "expr") {
  if (is.null(x) || length(x) == 0) return(alt)
  x <- as.character(x)[1]
  if (!nzchar(x)) alt else x
}