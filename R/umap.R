umapUI <- function(id) {
  ns <- NS(id)
  tagList(
    hr(), h4("UMAP (by gene set)"),
    fluidRow(
      column(4,
        selectInput(ns("gene_source"), "Gene set source",
          choices = c("SPARK-X (q ≤ FDR)" = "sparkx",
                      "Manual select genes" = "manual",
                      "Top-N by total counts" = "heg",
                      "All nonzero genes" = "all"),
          selected = "sparkx"
        ),
        conditionalPanel(
          sprintf("input['%s'] == 'sparkx'", ns("gene_source")),
          numericInput(ns("fdr"), "FDR cutoff", value = 0.05, min = 0, step = 0.01)
        ),
        conditionalPanel(
          sprintf("input['%s'] == 'manual'", ns("gene_source")),
          selectizeInput(ns("genes_manual"), "Pick genes (multi)",
            choices = NULL, multiple = TRUE, options = list(placeholder = "Type to search...")
          )
        ),
        conditionalPanel(
          sprintf("input['%s'] == 'heg'", ns("gene_source")),
          numericInput(ns("topN"), "Top-N by total counts", value = 1000, min = 10, step = 50)
        ),
        checkboxInput(ns("log1p"), "log1p transform", value = TRUE),
        checkboxInput(ns("libnorm"), "Library-size normalize", value = TRUE),
        checkboxInput(ns("zscore"), "Z-score per gene", value = TRUE)
      ),
      column(4,
        numericInput(ns("n_neighbors"), "n_neighbors", value = 15, min = 2, step = 1),
        numericInput(ns("min_dist"), "min_dist", value = 0.1, min = 0, step = 0.05),
        numericInput(ns("n_components"), "n_components", value = 2, min = 2, step = 1),
        numericInput(ns("seed"), "random_state (seed)", value = 42, min = 0, step = 1),
        actionButton(ns("run_umap"), "Run UMAP", class = "btn-primary")
      ),
      column(4,
        downloadButton(ns("dl_umap_csv"), "Download UMAP CSV"),
        br(), br(),
        shiny::verbatimTextOutput(ns("umap_msg"))
      )
    ),
        # 在参数区加一个“上色基因”选择器
        selectizeInput(ns("color_gene"), "Color by gene (expression):",
        choices = NULL, multiple = FALSE,
        options = list(placeholder = "Pick a gene to color")
        ),

        # 把原来的 plotOutput 改成 plotlyOutput（高度随意）
        plotly::plotlyOutput(ns("umap_plot"), height = "480px"),

    DT::DTOutput(ns("umap_table"))
  )
}




umapServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 同步可选基因（手动/上色两处）
    observe({
      if (!is.null(rv$gene_names)) {
        updateSelectizeInput(session, "genes_manual", choices = rv$gene_names, server = TRUE)
        updateSelectizeInput(session, "color_gene",  choices = rv$gene_names, server = TRUE)
      }
    })

    .safe_msg <- function(txt) output$umap_msg <- renderText({ txt })

    # —— 选基因：从 SPARK-X / 手动 / HEG / ALL —— #
    pick_genes <- reactive({
      req(rv$counts, rv$gene_names)
      src <- input$gene_source

      if (src == "sparkx") {
        df  <- as.data.frame(rv$sparkx_res)
        nms <- names(df); ln <- tolower(nms)
        gene_col <- nms[match("gene", ln)]; q_col <- nms[match("qval", ln)]
        fdr <- as.numeric(input$fdr); if (!is.finite(fdr)) fdr <- 0.05
        g <- df[[gene_col]][df[[q_col]] <= fdr]
        g <- intersect(as.character(g), rv$gene_names)
        return(g)

      } else if (src == "manual") {
        g <- intersect(as.character(input$genes_manual), rv$gene_names)
        return(g)

      } else if (src == "heg") {
        topN  <- max(10L, as.integer(input$topN))
        score <- Matrix::rowSums(rv$counts)
        ord   <- order(score, decreasing = TRUE)
        g     <- rv$gene_names[ord[seq_len(min(topN, length(ord)))]]
        return(g)

      } else { # "all"
        score <- Matrix::rowSums(rv$counts)
        g     <- rv$gene_names[score > 0]
        return(g)
      }
    })

    # —— 构建矩阵（spots x genes），并做预处理 —— #
    build_matrix <- reactive({
      req(rv$counts, rv$barcodes)
      g   <- pick_genes()
      idx <- match(g, rv$gene_names)
      idx <- idx[!is.na(idx)]
      X   <- t(as.matrix(rv$counts[idx, , drop = FALSE]))  # spots x genes

      # 库大小归一化（按行/spot）
      if (isTRUE(input$libnorm)) {
        lib <- rowSums(X)
        sf  <- lib / stats::median(lib[lib > 0])
        sf[!is.finite(sf) | sf <= 0] <- 1
        X <- sweep(X, 1, sf, "/")
      }

      if (isTRUE(input$log1p)) {
        X <- log1p(X)
      }

      # 去掉方差为 0 的基因列
      v <- apply(X, 2, stats::var, na.rm = TRUE)
      keep <- which(is.finite(v) & v > 0)
      X <- X[, keep, drop = FALSE]

      # 每基因 z-score
      if (isTRUE(input$zscore)) {
        X <- scale(X)
        X[!is.finite(X)] <- 0
      }
      X
    })

    # —— 运行 UMAP —— #
    run_umap <- eventReactive(input$run_umap, {
      req(rv$barcodes)
      X <- build_matrix()

      set.seed(if (is.finite(as.numeric(input$seed))) as.numeric(input$seed) else 42)
      emb <- uwot::umap(
        X,
        n_neighbors  = as.integer(input$n_neighbors),
        min_dist     = as.numeric(input$min_dist),
        n_components = as.integer(input$n_components),
        verbose      = FALSE,
        ret_model    = FALSE
      )
      emb <- as.data.frame(emb)
      colnames(emb) <- paste0("UMAP", seq_len(ncol(emb)))
      emb$barcode <- rv$barcodes
      emb
    })

    # 写回 rv 以便其他模块复用
    observeEvent(run_umap(), {
      rv$umap_embed <- run_umap()
    })

    # —— 文本信息 —— #
    output$umap_msg <- renderText({
      g <- pick_genes()
      X <- build_matrix()
      sprintf("Selected genes: %d | Matrix dims (spots x genes) after filtering: %d x %d",
              length(g), nrow(X), ncol(X))
    })

    # —— 计算用于上色的表达（与 emb$barcode 对齐；仅做 libnorm/log1p，与 build_matrix 保持一致） —— #
    color_expr <- reactive({
      emb <- run_umap(); req(emb, rv$counts, rv$gene_names, rv$barcodes)
      g <- input$color_gene
      if (is.null(g) || !nzchar(g) || !(g %in% rv$gene_names)) {
        return(rep(NA_real_, nrow(emb)))
      }
      gi <- match(g, rv$gene_names)
      v  <- as.numeric(rv$counts[gi, , drop = TRUE])
      names(v) <- rv$barcodes
      v  <- v[emb$barcode]

      if (isTRUE(input$libnorm)) {
        lib <- Matrix::colSums(rv$counts)
        names(lib) <- rv$barcodes
        sf <- lib[emb$barcode] / stats::median(lib[lib > 0])
        sf[!is.finite(sf) | sf <= 0] <- 1
        v <- v / sf
      }
      if (isTRUE(input$log1p)) v <- log1p(v)
      v[!is.finite(v)] <- NA_real_
      v
    })

    # —— UMAP Plot（Plotly + Viridis 色系；其余配置最简） —— #
    output$umap_plot <- plotly::renderPlotly({
      emb  <- run_umap(); req(emb)
      expr <- color_expr()
      gene_label <- if (!is.null(input$color_gene) && nzchar(input$color_gene)) input$color_gene else "expression"

      df <- data.frame(
        x = emb$UMAP1, y = emb$UMAP2,
        barcode = emb$barcode,
        expr = expr
      )

      plotly::plot_ly(
        df, x = ~x, y = ~y, type = "scattergl", mode = "markers",
        marker = list(
          size = 4, opacity = 0.8,
          color = ~expr, colorscale = "Viridis",
          colorbar = list(title = gene_label)
        ),
        text = ~paste0("barcode: ", barcode,
                       ifelse(is.finite(expr),
                              paste0("<br>expr: ", signif(expr, 4)),
                              "<br>expr: NA")),
        hoverinfo = "text"
      )
    })

    # 表格预览
    output$umap_table <- DT::renderDT({
      emb <- run_umap(); req(emb)
      DT::datatable(emb, options = list(pageLength = 25), rownames = FALSE)
    })

    # 下载
    output$dl_umap_csv <- downloadHandler(
      filename = function() {
        src <- input$gene_source
        tag <- switch(src,
          sparkx = sprintf("sparkx_fdr%s", input$fdr),
          manual = "manual",
          heg    = sprintf("topN%d", as.integer(input$topN)),
          all    = "all",
          "genes"
        )
        paste0("umap_", tag, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        emb <- run_umap()
        utils::write.csv(emb, file, row.names = FALSE)
      }
    )
  })
}
