# 需要的包（确保已安装：clusterProfiler、enrichplot、org.Hs.eg.db/org.Mm.eg.db、DOSE）
safely_load_orgdb <- function(org) {
  if (org == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("请先安装 org.Hs.eg.db")
    return("org.Hs.eg.db")
  } else if (org == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("请先安装 org.Mm.eg.db")
    return("org.Mm.eg.db")
  } else {
    stop("Unsupported organism: ", org)
  }
}

# 统一：把符号映射到 ENTREZID；返回 data.frame(symbol, ENTREZID)，去重
map_symbols_to_entrez <- function(symbols, organism = c("human","mouse")) {
  organism <- match.arg(organism)
  orgpkg <- safely_load_orgdb(organism)
  suppressMessages({
    df <- clusterProfiler::bitr(
      gsub("\\s+$","", symbols),
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = orgpkg
    )
  })
  unique(df[!is.na(df$ENTREZID), c("SYMBOL","ENTREZID")])
}

# 根据选择的 p 列取 Top N 基因；返回列表：top_df、universe_symbols
pick_top_genes_by_p <- function(res_df, p_col, top_n = 20) {
  stopifnot(p_col %in% names(res_df))
  x <- res_df
  # 只保留有 p 的行
  x <- x[is.finite(x[[p_col]]), , drop = FALSE]
  # 保护 0/负值
  x[[p_col]][x[[p_col]] <= 0] <- .Machine$double.xmin
  # FDR
  x$FDR <- p.adjust(x[[p_col]], method = "BH")
  x <- x[order(x[[p_col]], x$FDR), , drop = FALSE]
  # Top N（保留基因名和关键列）
  keep_cols <- intersect(c("gene","status","dev_expl", p_col, "FDR"), names(x))
  top_df <- head(x[, keep_cols, drop = FALSE], top_n)
  list(top_df = top_df, universe_symbols = x$gene)
}

# 运行 GO 富集（BP 为例）；返回 enrichResult
run_go_enrichment <- function(top_symbols, universe_symbols, organism = c("human","mouse")) {
  organism <- match.arg(organism)
  orgpkg <- safely_load_orgdb(organism)

  mapped_top <- map_symbols_to_entrez(top_symbols, organism)
  mapped_bg  <- map_symbols_to_entrez(universe_symbols, organism)

  if (nrow(mapped_top) == 0) stop("Top 基因没有成功映射到 ENTREZID。")
  if (nrow(mapped_bg)  == 0) stop("背景集合没有成功映射到 ENTREZID。")

  clusterProfiler::enrichGO(
    gene          = mapped_top$ENTREZID,
    universe      = unique(mapped_bg$ENTREZID),
    OrgDb         = orgpkg,
    keyType       = "ENTREZID",
    ont           = "BP",           # 也可切换 CC/MF
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,              # 显示所有，前端自己过滤
    qvalueCutoff  = 1,
    readable      = TRUE
  )
}