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


# ===== 轨迹模块：读取 dist CSV → 合并表达 → 画 LOESS / 表格 / 下载 =====
# 依赖：rv$counts / rv$pos / input$gene / expr_vec()
# 可选外部函数：merge_dist_expr, bin_summary, make_expr_from_counts, safe_gene_label
# 若未提供上述函数，模块内有兜底实现（尽量不改变你现有行为）

trajServer <- function(id, rv, expr_vec) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a

    # ---- 兜底：列名清洗 ----
    .clean_names <- function(nm) {
      nm <- trimws(nm)
      sub("^ï\\.+", "", nm)
    }

    # ---- 兜底：表达从 counts 计算（单基因；log1p 在外层控制）----
    .make_expr_from_counts <- function(gene, counts, barcodes) {
      if (is.null(gene) || !nzchar(gene)) return(rep(NA_real_, length(barcodes)))
      rn <- rownames(counts)
      gn <- if (!is.null(rn)) rn else colnames(counts) # 以防你的 counts 是 genes x spots（默认如此）
      idx <- match(gene, gn)
      if (is.na(idx)) return(rep(NA_real_, length(barcodes)))
      ord <- match(barcodes, colnames(counts))
      if (anyNA(ord)) return(rep(NA_real_, length(barcodes)))
      as.numeric(counts[idx, ord, drop = TRUE])
    }

    # ---- 兜底：合并距离与表达 ----
    .merge_dist_expr <- function(df, expr_from_counts = NULL) {
      if (is.null(expr_from_counts)) {
        if (!"expr" %in% names(df)) df$expr <- NA_real_
        return(df)
      }
      df$expr <- expr_from_counts
      df
    }

    # ---- 兜底：分箱汇总 ----
    .bin_summary <- function(x, y, bins = 30L) {
      ok <- is.finite(x) & is.finite(y)
      x <- x[ok]; y <- y[ok]
      if (!length(x)) return(data.frame(bin=integer(), x_mid=numeric(), n=integer(), y_mean=numeric(), y_sd=numeric()))
      br <- pretty(range(x), n = bins)
      cutx <- cut(x, breaks = br, include.lowest = TRUE, labels = FALSE)
      agg <- function(z) c(mean = mean(z), sd = stats::sd(z), n = length(z))
      tmp <- do.call(rbind, tapply(y, cutx, agg))
      mids <- (br[-1] + br[-length(br)])/2
      data.frame(
        bin   = as.integer(names(tapply(y, cutx, length))),
        x_mid = mids[as.integer(names(tapply(y, cutx, length)))],
        n     = as.integer(tmp[, "n"]),
        y_mean= as.numeric(tmp[, "mean"]),
        y_sd  = as.numeric(tmp[, "sd"]),
        row.names = NULL
      )
    }

    # 若外部实现存在则优先用外部
    make_expr_from_counts <- get0("make_expr_from_counts", mode = "function") %||% .make_expr_from_counts
    merge_dist_expr       <- get0("merge_dist_expr",       mode = "function") %||% .merge_dist_expr
    bin_summary           <- get0("bin_summary",           mode = "function") %||% .bin_summary
    safe_gene_label       <- get0("safe_gene_label",       mode = "function") %||%
      function(g, fallback="expr") if (is.null(g) || !nzchar(g)) fallback else g

    # ========== 1) 读取距离 CSV ==========
    dist_df <- reactive({
      req(input$dist_csv)
      path <- input$dist_csv$datapath
      df <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
      names(df) <- .clean_names(names(df))
      df
    })

    # ========== 2) 选择距离列 ==========
    output$traj_dist_col_ui <- renderUI({
      df <- dist_df()
      dcols <- grep("^dist_", names(df), value = TRUE)
      selectInput(session$ns("traj_dist_col"), "Select distance column", choices = dcols, selected = dcols[1])
    })

    # ========== 3) 合并表达 ==========
    traj_merged <- reactive({
      df <- dist_df(); req(input$traj_dist_col)
      expr_cnt <- NULL
      if (!is.null(rv$counts) && !is.null(rv$pos) && isTruthy(input$gene)) {
        bc  <- df$barcode
        ord <- match(bc, colnames(rv$counts))
        if (!anyNA(ord)) {
          expr_cnt <- make_expr_from_counts(
            gene     = input$gene,
            counts   = rv$counts[, ord, drop = FALSE],
            barcodes = bc
          )
        }
      }
      merge_dist_expr(df, expr_from_counts = expr_cnt)
    })

    # ========== 4) 绘图 ==========
    output$traj_plot <- renderPlot({
      df <- traj_merged(); req(input$traj_dist_col)
      dvar <- input$traj_dist_col
      ok <- is.finite(df[[dvar]]) & is.finite(df$expr)

      # 分箱表缓存到 session（不污染 .GlobalEnv）
      tbl <- bin_summary(df[[dvar]][ok], df$expr[ok], bins = input$traj_bins %||% 30L)
      session$userData$traj_table_cache <- tbl

      gene_lab <- safe_gene_label(input$gene, "expr")
      ylab_txt <- if (isTRUE(input$log1p)) paste0("log1p(", gene_lab, ")") else gene_lab

      library(ggplot2)
      g <- ggplot(data.frame(x = df[[dvar]][ok], y = df$expr[ok]), aes(x, y))
      if (isTRUE(input$traj_show_points)) g <- g + geom_point(alpha = 0.25, size = 0.8)
      g <- g +
        geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = input$traj_span %||% 0.75) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        labs(x = dvar, y = ylab_txt, title = paste("Expression vs", dvar)) +
        theme_minimal(base_size = 12)
      print(g)
    })

    # ========== 5) 分箱表 ==========
    output$traj_table <- renderTable({
      tbl <- session$userData$traj_table_cache
      if (is.null(tbl)) data.frame() else tbl
    }, striped = TRUE, bordered = TRUE, hover = TRUE)

    # ========== 6) 下载：图 ==========
    output$dl_traj_plot <- downloadHandler(
      filename = function() {
        gene_lab <- safe_gene_label(input$gene, "expr")
        paste0("traj_", input$traj_dist_col %||% "dist", "_", gene_lab, "_",
               format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        df <- traj_merged(); req(input$traj_dist_col)
        dvar <- input$traj_dist_col
        ok <- is.finite(df[[dvar]]) & is.finite(df$expr)

        gene_lab <- safe_gene_label(input$gene, "expr")
        ylab_txt <- if (isTRUE(input$log1p)) paste0("log1p(", gene_lab, ")") else gene_lab

        library(ggplot2)
        png(file, width = 1600, height = 1000, res = 150)
        g <- ggplot(data.frame(x = df[[dvar]][ok], y = df$expr[ok]), aes(x, y))
        if (isTRUE(input$traj_show_points)) g <- g + geom_point(alpha = 0.25, size = 0.8)
        g <- g +
          geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = input$traj_span %||% 0.75) +
          geom_vline(xintercept = 0, linetype = "dashed") +
          labs(x = dvar, y = ylab_txt, title = paste("Expression vs", dvar)) +
          theme_minimal(base_size = 14)
        print(g)
        dev.off()
      }
    )

    # ========== 7) 下载：表 ==========
    output$dl_traj_table <- downloadHandler(
      filename = function() paste0("traj_bins_", input$traj_dist_col %||% "dist", "_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        tbl <- session$userData$traj_table_cache
        utils::write.csv(tbl %||% data.frame(), file, row.names = FALSE)
      }
    )
  })
}


trajUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Expression trajectory along distance (LOESS)"),
    fileInput(ns("dist_csv"), "Upload distance CSV (must have 'barcode' and 'dist_*' columns)", accept = ".csv"),
    fluidRow(
      column(4, uiOutput(ns("traj_dist_col_ui"))),
      column(3, numericInput(ns("traj_bins"), "Number of bins", value = 50, min = 10, step = 10)),
      column(3, sliderInput(ns("traj_span"), "LOESS span", min = 0.1, max = 1, value = 0.4, step = 0.05)),
      column(2, checkboxInput(ns("traj_show_points"), "Show scatter points", TRUE))
    ),
    plotOutput(ns("traj_plot"), height = "420px"),
    fluidRow(
      column(6, downloadButton(ns("dl_traj_plot"), "Download trajectory (PNG)")),
      column(6, downloadButton(ns("dl_traj_table"), "Download binning stats (CSV)"))
    ),
    tableOutput(ns("traj_table"))
  )
}
