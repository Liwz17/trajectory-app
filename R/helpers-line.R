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


make_line_ab <- function(Ax, Ay, Bx, By, label=NULL) {
  dx <- Bx - Ax; dy <- By - Ay; L <- sqrt(dx*dx + dy*dy)
  if (!is.finite(L) || L <= 0) return(NULL)
  list(Ax=Ax, Ay=Ay, Bx=Bx, By=By, dx=dx, dy=dy, L=L, label=label %||% "line")
}

project_points_to_line <- function(x, y, line) {
  Ax <- line$Ax; Ay <- line$Ay; dx <- line$dx; dy <- line$dy; L <- line$L
  t <- ((x - Ax) * dx + (y - Ay) * dy) / L
  d_perp <- ((x - Ax) * dy - (y - Ay) * dx) / L
  data.frame(t=t, d_perp=d_perp)
}


sel_dragmode <- function(sel_mode) if (isTruthy(sel_mode) && sel_mode=="select") "select" else "lasso"

# 这里放 build_he_plotly_figure / make_band_shape / make_line_shape / make_arrow_annot


lineModuleServer <- function(id, rv, img_obj, expr_vec) {
  moduleServer(id, function(input, output, session) {
    line_current <- reactiveValues(A=NULL, B=NULL)
    lines <- reactiveVal(data.frame(id=character(), label=character(), Ax=numeric(), Ay=numeric(), Bx=numeric(), By=numeric()))

    observeEvent(plotly::event_data("plotly_click", source="main"), {
      req(isTRUE(input$line_draw))
      ed <- plotly::event_data("plotly_click", source="main")
      if (is.null(ed) || nrow(ed)==0) return()
      pt <- c(ed$x[1], ed$y[1])
      if (is.null(line_current$A)) line_current$A <- pt
      else if (is.null(line_current$B)) line_current$B <- pt
      else { line_current$A <- pt; line_current$B <- NULL }
    })
    observeEvent(input$line_commit, {
      req(line_current$A, line_current$B)
      id <- paste0("L", as.integer(Sys.time()))
      label <- if (nzchar(input$line_label_new)) input$line_label_new else sprintf("line_%s", format(Sys.time(), "%H%M%S"))
      df <- lines()
      df <- rbind(df, data.frame(id=id, label=label, Ax=line_current$A[1], Ay=line_current$A[2], Bx=line_current$B[1], By=line_current$B[2]))
      lines(df); line_current$A <- NULL; line_current$B <- NULL; updateTextInput(session, "line_label_new", value="")
      updateSelectInput(session, "line_active", choices=setNames(df$id, df$label), selected=id)
    })
    observeEvent(input$line_clear_all, {
      lines(data.frame(id=character(), label=character(), Ax=numeric(), Ay=numeric(), Bx=numeric(), By=numeric()))
      updateSelectInput(session, "line_active", choices=character(0), selected=NULL)
      line_current$A <- NULL; line_current$B <- NULL
    })

    active_line <- reactive({
      df <- lines(); req(nrow(df)>0, input$line_active)
      r <- df[df$id == input$line_active, , drop=FALSE]; req(nrow(r)==1)
      make_line_ab(r$Ax, r$Ay, r$Bx, r$By, label=r$label)
    })

    shapes <- reactive({
      sh <- list()
      if (!is.null(line_current$A) && !is.null(line_current$B)) {
        ln <- make_line_ab(line_current$A[1], line_current$A[2], line_current$B[1], line_current$B[2])
        if (!is.null(ln)) { sh <- c(sh, list(make_band_shape(ln, as.numeric(input$line_band))), list(make_line_shape(ln))) }
      }
      df <- lines()
      if (nrow(df)) {
        for (i in seq_len(nrow(df))) {
          ln <- make_line_ab(df$Ax[i], df$Ay[i], df$Bx[i], df$By[i], label=df$label[i])
          if (!is.null(ln)) { sh <- c(sh, list(make_band_shape(ln, as.numeric(input$line_band))), list(make_line_shape(ln))) }
        }
      }
      if (length(sh)) sh else NULL
    })

    annotations <- reactive({
      an <- list()
      if (!is.null(line_current$A) && !is.null(line_current$B)) {
        ln <- make_line_ab(line_current$A[1], line_current$A[2], line_current$B[1], line_current$B[2])
        if (!is.null(ln)) an <- c(an, list(make_arrow_annot(ln)))
      }
      df <- lines()
      if (nrow(df) && !is.null(input$line_active) && input$line_active %in% df$id) {
        r <- df[df$id == input$line_active, , drop=FALSE]
        ln <- make_line_ab(r$Ax, r$Ay, r$Bx, r$By, label=r$label)
        if (!is.null(ln)) an <- c(an, list(make_arrow_annot(ln)))
      }
      if (length(an)) an else NULL
    })

    return(list(
      active_line  = active_line,
      shapes       = shapes,
      annotations  = annotations
    ))
  })
}

  # 可选：沿线曲线子模块
  lineCurveServer <- function(id, rv, expr_vec, active_line) {
    moduleServer(id, function(input, output, session) {
      line_sample_df <- reactive({
        req(rv$pos, expr_vec(), input$line_band, active_line())
        sample_along_line(rv$pos, expr_vec(), active_line(), as.numeric(input$line_band), end_tol = 5)
      })
      output$line_expr_plot <- renderPlot({
        df <- line_sample_df(); ok <- is.finite(df$t) & is.finite(df$expr)
        ylab_txt <- if (isTRUE(input$log1p)) paste0("log1p(", safe_gene_label(input$gene, "expr"), ")") else safe_gene_label(input$gene, "expr")
        ggplot2::ggplot(df[ok,], ggplot2::aes(t, expr)) +
          ggplot2::geom_point(alpha=.25, size=.8) +
          ggplot2::geom_smooth(method="loess", formula = y ~ x, se=TRUE, span=input$line_span) +
          ggplot2::labs(x="Distance along line (pixels)", y=ylab_txt) +
          ggplot2::theme_minimal(base_size=12)
      })

      output$dl_line_expr_plot <- downloadHandler(
        filename = function() paste0("line_expr_", active_line()$label, "_", input$gene %||% "expr", "_",
                                    format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
        content = function(file) {
          df <- line_sample_df(); ok <- is.finite(df$t) & is.finite(df$expr)
          gene_label <- safe_gene_label(input$gene, "expr")
          ylab_txt   <- if (isTRUE(input$log1p)) paste0("log1p(", gene_label, ")") else gene_label
          png(file, width = 1600, height = 1000, res = 150)
          library(ggplot2)
          print(
            ggplot(df[ok,], aes(t, expr)) +
              geom_point(alpha = 0.25, size = 0.8) +
              geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = input$line_span) +
              labs(x = "Distance along line (pixels)", y = ylab_txt) +
              theme_minimal(base_size = 14)
          )
          dev.off()
        }
      )
      output$dl_line_expr_table <- downloadHandler(
        filename = function() paste0("line_sampling_", (active_line()$label %||% "line"), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
        content  = function(file) {
          df <- line_sample_df()
          utils::write.csv(df[, c("barcode","x_hires","y_hires","t","d_perp","expr")], file, row.names = FALSE)
        }
      )
    })
  }


# ---- 模块定义 ----
svgLineServer <- function(id, rv, active_line) {
  moduleServer(id, function(input, output, session) {
    
    analyze_genes_along_line <- function(counts, genes, covar_t, fam_nb, 
                                         k = 6L, top_n_expressed = 1000, min_nnz = 50) {
      distance_vars <- "t"
      analyze_spatial_genes(
        counts            = counts,
        genes             = genes,
        covar             = covar_t,
        distance_vars     = distance_vars,
        fam_nb            = fam_nb,
        use_interaction   = FALSE,
        k                 = k,
        min_nnz           = min_nnz,
        min_sum           = 0,
        min_nnz_gene      = 1L,
        top_n_expressed   = top_n_expressed
      )
    }
    
    observeEvent(input$run_svg, {
      req(rv$counts, rv$gene_names, rv$pos, active_line())
      
      withProgress(message = "Running SVG along line ...", value = 0, {
        ln  <- active_line()
        pr  <- project_points_to_line(rv$pos$x_hires, rv$pos$y_hires, ln)
        cov <- data.frame(barcode = rv$pos$barcode, t = pr$t, stringsAsFactors = FALSE)
        
        band <- as.numeric(input$line_band)
        keep <- abs(pr$d_perp) <= band & is.finite(pr$t)
        cov  <- cov[keep, , drop = FALSE]
        
        ord <- match(cov$barcode, colnames(rv$counts))
        cov  <- cov[!is.na(ord), , drop = FALSE]
        ord  <- ord[!is.na(ord)]
        counts_aligned <- rv$counts[, ord, drop = FALSE]
        
        # ---- library size normalization offset ----
        libsize <- Matrix::colSums(counts_aligned)
        sf <- as.numeric(libsize / stats::median(libsize[libsize > 0]))
        cov$log_off <- log(pmax(sf, 1e-8))
        
        fam_nb <- mgcv::nb(link="log")
        res <- analyze_genes_along_line(
          counts = counts_aligned,
          genes = rv$gene_names,
          covar_t = cov,
          fam_nb = fam_nb,
          k = as.integer(input$svg_k),
          top_n_expressed = as.integer(input$svg_topN),
          min_nnz = as.integer(input$svg_min_nnz)
        )
        
        rv$svg_line_result <- res
        showNotification(sprintf("SVG 完成：%d 行", nrow(res)), type = "message")
      })
    })
    
    output$dl_svg_results <- downloadHandler(
      filename = function() {
        # 1) 取 label（可能是 NULL、长度>1、或包含不适合文件名的字符）
        lab <- try(active_line()$label, silent = TRUE)

        # 2) 统一收敛成单个标量字符；为空则用 "line"
        if (inherits(lab, "try-error") || is.null(lab) || length(lab) == 0) {
          lab <- "line"
        } else {
          lab <- as.character(lab)[1]   # 只取第一个
        }

        # 3) 清洗为安全文件名（去掉/替换非法字符）
        lab <- gsub("[^A-Za-z0-9._-]+", "_", trimws(lab))

        paste0(
          "svg_line_", lab, "_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".csv"
        )
      },
      content = function(file) {
        # 如果结果还没算好，给出空表但不报错（也可改成 req(rv$svg_line_result) 强制有结果才下载）
        res <- rv$svg_line_result
        if (is.null(res)) {
          res <- data.frame()
        } else {
          # 确保是 data.frame（有些函数返回 tibble / data.table）
          res <- as.data.frame(res)
        }

        utils::write.csv(res, file, row.names = FALSE)
      }
    )

  })
}


svgLineUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("SVG (1D GAM) along line"),
    fluidRow(
      column(3, numericInput(ns("svg_topN"), "Top expressed genes (N)", value = 1000, min = 100, step = 100)),
      column(3, numericInput(ns("svg_min_nnz"), "Min nonzero per gene", value = 50, min = 1, step = 5)),
      column(3, numericInput(ns("svg_k"), "Spline k (df)", value = 6L, min = 3L, step = 1L)),
      column(3, checkboxInput(ns("svg_interact"), "Try interaction (ignored for 1D)", value = FALSE))
    ),
    fluidRow(
      column(6, actionButton(ns("run_svg"), "Run SVG (along line)", class = "btn-primary")),
      column(6, downloadButton(ns("dl_svg_results"), "Download SVG results (CSV)"))
    )
  )
}


# ===== A→B 画线 + 沿线曲线：UI =====
lineUI <- function(id) {
  ns <- NS(id)
  tagList(
    hr(), h3("User-drawn line (A→B) analysis"),
    fluidRow(
      column(3, checkboxInput(ns("line_draw"), "Enable draw-line mode", value = FALSE)),
      column(3, numericInput(ns("line_band"), "Band width (pixels)", value = 50, min = 1, step = 5)),
      column(3, numericInput(ns("line_span"), "LOESS span", value = 0.3, min = 0.05, max = 1, step = 0.05)),
      column(3, actionButton(ns("line_clear_all"), "Clear all lines"))
    ),
    fluidRow(
      column(6, textInput(ns("line_label_new"), "New line label (optional)", "")),
      column(3, actionButton(ns("line_commit"), "Commit current A→B")),
      column(3, actionButton(ns("line_reset_ab"), "Reset A/B"))
    ),
    verbatimTextOutput(ns("line_ab_info")),
    fluidRow(
      column(6, tableOutput(ns("line_list_tbl"))),
      column(6, selectInput(ns("line_active"), "Pick a line for analysis", choices = character(0)))
    ),

    # —— 沿线表达曲线 —— #
    hr(), h4("Expression curve along line"),
    plotOutput(ns("line_expr_plot"), height = "360px"),
    fluidRow(
      column(6, downloadButton(ns("dl_line_expr_plot"),  "Download curve (PNG)")),
      column(6, downloadButton(ns("dl_line_expr_table"), "Download sampled table (CSV)"))
    )
  )
}
