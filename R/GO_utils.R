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


# # ===== UI（最小化；如你已有 UI，可不使用这个） =====
# goUI <- function(id, title = "GO enrichment") {
#   ns <- NS(id)
#   tagList(
#     h4(title),
#     uiOutput(ns("dist_col_ui")),
#     numericInput(ns("topn_go"), "Top N genes", value = 200, min = 10, step = 10),
#     selectInput(ns("go_org"), "Organism", choices = c("org.Hs.eg.db","org.Mm.eg.db"), selected = "org.Hs.eg.db"),
#     actionButton(ns("run_go"), "Run GO"),
#     tags$hr(),
#     h5("Top genes"),
#     tableOutput(ns("tbl_top_genes")),
#     h5("Enriched GO"),
#     tableOutput(ns("tbl_go")),
#     plotOutput(ns("plot_go"), height = "420px"),
#     downloadButton(ns("dl_go_table"), "Download GO table")
#   )
# }

# ===== Server 模块 =====
goServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    gam_df_r <- reactive({
      req(input$gam_csv)
      df <- tryCatch(
        utils::read.csv(input$gam_csv$datapath, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
      )
      req(!is.null(df))
      nm <- trimws(names(df)); nm <- sub("^ï\\.+", "", nm); names(df) <- nm
      if (!"gene" %in% names(df)) stop("CSV missing 'gene' column")
      df
    })

    # 选择 p 值列（兼容你原来的候选规则）
    output$dist_col_ui <- renderUI({
      df <- gam_df_r(); req(df)
      nm <- names(df)
      candidates <- c("p_t", "p_int", grep("^p_dist_", nm, value = TRUE))
      pcols <- unique(candidates[candidates %in% nm])
      if (length(pcols) == 0) pcols <- grep("^p(_|$)", nm, value = TRUE)
      selectInput(session$ns("dist_col"), "Select p-value column",
                  choices = pcols, selected = pcols[1])
    })

    # 计算 Top 基因 + 背景集
    top_pick <- reactive({
      df <- gam_df_r(); req(df, input$dist_col, input$topn_go)
      x <- df[is.finite(suppressWarnings(as.numeric(df[[input$dist_col]]))), , drop = FALSE]
      # 数值化 + 最小正数
      pnum <- suppressWarnings(as.numeric(x[[input$dist_col]]))
      pnum[pnum <= 0 | !is.finite(pnum)] <- .Machine$double.xmin
      x[[input$dist_col]] <- pnum
      x$FDR <- p.adjust(pnum, method = "BH")
      x <- x[order(pnum, x$FDR), , drop = FALSE]
      keep_cols <- intersect(c("gene","status","dev_expl", input$dist_col, "FDR"), names(x))
      list(
        top_df = head(x[, keep_cols, drop = FALSE], input$topn_go),
        universe_symbols = x$gene
      )
    })

    output$tbl_top_genes <- renderTable({
      tp <- top_pick(); tp$top_df
    }, striped = TRUE, bordered = TRUE, hover = TRUE)

    go_res <- eventReactive(input$run_go, {
      tp <- top_pick(); req(input$go_org)
      run_go_enrichment(
        top_symbols      = tp$top_df$gene,
        universe_symbols = tp$universe_symbols,
        organism         = input$go_org
      )
    })

    output$tbl_go <- renderTable({
      eg <- go_res()
      df <- try(as.data.frame(eg), silent = TRUE)
      if (inherits(df, "try-error") || !nrow(df)) return(data.frame())
      keep <- intersect(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count"), names(df))
      df[, keep, drop = FALSE]
    }, striped = TRUE, bordered = TRUE, hover = TRUE)

    output$plot_go <- renderPlot({
      eg <- go_res()
      df <- try(as.data.frame(eg), silent = TRUE)
        if (!requireNamespace("enrichplot", quietly = TRUE)) {
        plot.new(); title("Package 'enrichplot' not installed"); return()
      }
      enrichplot::dotplot(eg, showCategory = min(15, nrow(df)))
    })

    output$dl_go_table <- downloadHandler(
      filename = function() paste0("go_enrich_", input$dist_col %||% "p", "_top", input$topn_go %||% "N", "_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        eg <- go_res()
        utils::write.csv(as.data.frame(eg), file, row.names = FALSE)
      }
    )
  })
}



goUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("GAM results → Top genes → GO enrichment"),
    fileInput(ns("gam_csv"), "Upload gam_results_*.csv (must contain 'gene' and 'p_dist_*' columns)", accept = ".csv"),
    fluidRow(
      column(4, uiOutput(ns("dist_col_ui"))),
      column(3, numericInput(ns("topn_go"), "Top N genes", 20, min = 5, step = 5)),
      column(3, selectInput(ns("go_org"), "Species", c(Human = "human", Mouse = "mouse"), "human")),
      column(2, actionButton(ns("run_go"), "Run GO enrichment", class = "btn-primary"))
    ),
    br(),
    h4("Top genes (ascending by selected p column)"),
    tableOutput(ns("tbl_top_genes")),
    br(),
    h4("GO enrichment results (BP)"),
    tableOutput(ns("tbl_go")),
    plotOutput(ns("plot_go"), height = "480px"),
    downloadButton(ns("dl_go_table"), "Download GO table")
  )
}
