# ---
# Shiny module to run DESpace using in‑app annotations as spatial domains (clusters)
#
# Requirements (add to your app):
#   library(DESpace)
#   library(SpatialExperiment)
#   library(SummarizedExperiment)
#   library(dplyr)
#   library(DT)
#   library(ggplot2)
#
# Assumptions about your app state:
#   - `rv$spe` holds a SpatialExperiment (or can be coerced to one) with barcodes in colnames(rv$spe)
#   - `rv$pos` is a data.frame with a `barcode` column (already in your app)
#   - `annotations` is a reactiveVal data.frame with columns: barcode, label (from your existing init/setup)
#   - If you have multiple samples, ensure colData(rv$spe)$sample is present for replicate mode
#
# You can drop these functions into your server/ui files and call `despace_ui("despace1")` +
#   `despace_server("despace1", rv, annotations)` to wire it up.
# ---

#===========================#
# UI
#===========================#
despace_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabPanel(
      title = "DESpace",
      fluidRow(
        column(
          width = 3,
          h3("Domains (clusters) from annotations"),
          uiOutput(ns("anno_labels_ui")),
          tags$small("Tip: add labels via your existing annotation panel; we use them as domains."),
          tags$hr(),
          checkboxInput(ns("use_replicates"), "Replicates (multi-sample joint test)", value = FALSE),
          uiOutput(ns("sample_col_ui")),
          numericInput(ns("min_spots"), "Min spots per domain", value = 20, min = 1, step = 1),
          tags$hr(),
          actionButton(ns("run_svg"), "Run DESpace: DESpace_test", class = "btn-primary"),
          br(), br(),
          actionButton(ns("run_individual"), "Run DESpace: individual_test"),
          br(),
          uiOutput(ns("domain_pick_ui"))
        ),
        column(
          width = 9,
          h4("Overall (domain-aware) result — DESpace_test()"),
          DTOutput(ns("svg_table")),
          br(),
          downloadButton(ns("dl_svg"), "Download DESpace_test table"),
          tags$hr(),
          h4("Domain-specific ranking — individual_test()"),
          numericInput(ns("top_n"), "Top N per domain", value = 25, min = 1, step = 1),
          DTOutput(ns("indiv_table")),
          br(),
          downloadButton(ns("dl_indiv"), "Download per-domain table"),
          tags$hr(),
          h4("Inspect a gene across domains"),
          textInput(ns("gene_pick"), "Gene symbol / row name", value = ""),
          plotOutput(ns("gene_violin"), height = 300)
        )
      )
    )
  )
}

#===========================#
# Helpers
#===========================#
.ensure_spe <- function(obj) {
  # Coerce SummarizedExperiment to SpatialExperiment if needed; otherwise assume ready
  if (inherits(obj, "SpatialExperiment")) return(obj)
  if (inherits(obj, "SummarizedExperiment")) {
    return(SpatialExperiment::SpatialExperiment(
      assays = assays(obj),
      rowData = rowData(obj),
      colData = colData(obj)
    ))
  }
  stop("rv$spe must be a SpatialExperiment or SummarizedExperiment with counts.")
}

.annot_to_cluster <- function(spe, annotations_df, min_spots = 1, 
                              cluster_col = "domain", other_label = "other") {
  stopifnot(all(c("barcode", "label") %in% names(annotations_df)))
  bc <- colnames(spe)
  m <- match(bc, annotations_df$barcode)
  cluster <- annotations_df$label[m]
  # ✅ 把 NA（没 annotate 的 spot）改成 "other"
  cluster[is.na(cluster)] <- other_label
  spe2 <- spe  # ✅ 不 subset，全保留
  cluster2 <- droplevels(factor(cluster))
  # ✅ min_spots 仍然生效 —— 但可以选保留 other
  tbl <- table(cluster2)
  big_lvls <- names(tbl)[tbl >= min_spots]
  spe2 <- spe2[, cluster2 %in% big_lvls]
  cluster2 <- droplevels(factor(cluster2[cluster2 %in% big_lvls]))
  if (nlevels(cluster2) < 2) stop("Need at least 2 domains with >= min_spots.")
  colData(spe2)[[cluster_col]] <- cluster2
  list(spe = spe2, cluster_col = cluster_col, domains = levels(cluster2))
}


.make_gene_violin <- function(spe, cluster_col, gene) {
  if (!gene %in% rownames(spe)) return(ggplot() + ggtitle("Gene not found in rownames(spe)"))
  # use log1p(counts) to visualize
  cnt <- assays(spe)[["counts"]]
  if (is.null(cnt)) cnt <- assays(spe)[[1]]
  y <- as.numeric(cnt[gene, ])
  df <- data.frame(expr = log1p(y), cluster = colData(spe)[[cluster_col]])
  ggplot(df, aes(x = cluster, y = expr)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.size = 0.5) +
    labs(x = "Domain", y = "log1p(counts)") +
    theme_minimal(base_size = 12)
}

#===========================#
# Server
#===========================#
despace_server <- function(id, rv, annotations) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Store results
    rvs <- reactiveValues(svg_res = NULL, svg_tbl = NULL, indiv_tbl = NULL, spe2 = NULL, cluster_col = NULL, domains = NULL)

    # UI: show current labels/domains (from annotations)
    output$anno_labels_ui <- renderUI({
      ann <- annotations()
      if (nrow(ann) == 0) {
        tagList(tags$em("No annotations yet. Use your annotation panel to add labels."))
      } else {
        tbl <- ann %>% count(label, name = "spots") %>% arrange(desc(spots))
        HTML(paste0("<b>Annotated domains:</b><br>", paste(sprintf("%s (%d)", tbl$label, tbl$spots), collapse = "<br>")))
      }
    })

    # Replicates mode: sample column selector if present in colData(rv$spe)
    output$sample_col_ui <- renderUI({
      if (!isTRUE(input$use_replicates)) return(NULL)
      spe <- rv$spe
      if (is.null(spe)) return(tags$small("No spe in rv yet."))
      spe <- try(.ensure_spe(spe), silent = TRUE)
      if (inherits(spe, "try-error")) return(tags$small("spe not ready."))
      cdn <- colnames(colData(spe))
      if (is.null(cdn) || length(cdn) == 0) return(tags$small("colData has no columns."))
      selectInput(ns("sample_col"), "Sample column in colData(spe)", choices = cdn, selected = "sample")
    })

    # Build a domain-annotated SPE whenever we run
    .build_spe2 <- reactive({
      req(rv$spe)
      spe <- .ensure_spe(rv$spe)
      ann <- annotations()
      res <- .annot_to_cluster(spe, ann, min_spots = input$min_spots, cluster_col = "domain")
      rvs$spe2 <- res$spe
      rvs$cluster_col <- res$cluster_col
      rvs$domains <- res$domains
      res$spe
    })

    observeEvent(rvs$domains, {
      updateSelectInput(session, "domain_pick", choices = rvs$domains)
    })

    output$domain_pick_ui <- renderUI({
      if (is.null(rvs$domains)) return(NULL)
      selectInput(ns("domain_pick"), "Choose a domain for per-domain ranking", choices = rvs$domains)
    })

    #--------------------------
    # Run DESpace_test (overall, domain-aware)
    #--------------------------
    observeEvent(input$run_svg, {
      spe2 <- .build_spe2()
      withProgress(message = "Running DESpace::DESpace_test()", value = 0, {
        incProgress(0.1)
        cluster_col <- rvs$cluster_col
        # args for single vs replicate mode
        sample_col <- NULL
        replicates <- FALSE
        if (isTRUE(input$use_replicates) && !is.null(input$sample_col) && input$sample_col %in% colnames(colData(spe2))) {
          sample_col <- input$sample_col
          replicates <- TRUE
        }
        # Run
        svg_res <- try(DESpace::DESpace_test(
          spe = spe2,
          spatial_cluster = cluster_col,
          verbose = TRUE
        ), silent = TRUE)
        if (inherits(svg_res, "try-error")) {
          showNotification(paste("DESpace_test error:", as.character(svg_res)), type = "error")
          return(NULL)
        }
        rvs$svg_res <- svg_res
        # Top results helper if available, else rank by p_val
        top_tbl <- try(DESpace::top_results(svg_res), silent = TRUE)
        if (inherits(top_tbl, "try-error") || is.null(top_tbl)) {
          # Fall back: assume svg_res is a data.frame with p_val / padj columns
          df <- as.data.frame(svg_res)
          ord <- if ("padj" %in% names(df)) order(df$padj, df$p_val) else if ("p_val" %in% names(df)) order(df$p_val) else seq_len(nrow(df))
          top_tbl <- df[ord, , drop = FALSE]
        }
        rvs$svg_tbl <- as.data.frame(top_tbl)
        incProgress(0.9)
      })
    })

    output$svg_table <- renderDT({
      req(rvs$svg_tbl)
      datatable(rvs$svg_tbl, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    })

    output$dl_svg <- downloadHandler(
      filename = function() paste0("despace_DESpace_test_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv"),
      content = function(file) write.csv(rvs$svg_tbl %||% data.frame(), file, row.names = FALSE)
    )

    #--------------------------
    # Run individual_test (per-domain)
    #--------------------------
    observeEvent(input$run_individual, {
      req(input$domain_pick)
      spe2 <- .build_spe2()
      cluster_col <- rvs$cluster_col
      dom <- input$domain_pick
      withProgress(message = sprintf("Running DESpace::individual_test() for %s", dom), value = 0, {
        incProgress(0.1)
        sample_col <- NULL
        replicates <- FALSE
        if (isTRUE(input$use_replicates) && !is.null(input$sample_col) && input$sample_col %in% colnames(colData(spe2))) {
          sample_col <- input$sample_col
          replicates <- TRUE
        }
        # Subset to the chosen domain vs rest to avoid returning all domains at once
        # keep <- colData(spe2)[[cluster_col]] %in% c(dom, levels(colData(spe2)[[cluster_col]]))
        # spe_sub <- spe2[, keep]
        # Run. Some versions of DESpace::individual_test return a list; capture robustly.
        indiv_res <- try(DESpace::individual_test(
          spe = spe2,
          edgeR_y = rvs$svg_res$estimated_y,
          spatial_cluster = cluster_col,
        ), silent = TRUE)
        if (inherits(indiv_res, "try-error")) {
          showNotification(paste("individual_test error:", as.character(indiv_res)), type = "error")
          return(NULL)
        }
        # Extract table for the chosen domain
        tbl <- NULL
        if (is.list(indiv_res) && !is.data.frame(indiv_res)) {
          # try by name
          if (!is.null(names(indiv_res)) && dom %in% names(indiv_res)) tbl <- indiv_res[[dom]]
          # or first element
          if (is.null(tbl) && length(indiv_res) > 0) tbl <- indiv_res[[1]]
        } else if (is.data.frame(indiv_res)) {
          tbl <- indiv_res
        }
        if (is.null(tbl)) {
          showNotification("Could not parse individual_test result; showing raw structure in console.", type = "warning")
          print(indiv_res)
          return(NULL)
        }
        # Keep top N
        ord <- if ("padj" %in% names(tbl)) order(tbl$padj, tbl$p_val) else if ("p_val" %in% names(tbl)) order(tbl$p_val) else seq_len(nrow(tbl))
        tbl <- tbl[ord, , drop = FALSE]
        if (nrow(tbl) > input$top_n) tbl <- tbl[seq_len(input$top_n), , drop = FALSE]
        rvs$indiv_tbl <- as.data.frame(tbl)
        incProgress(0.9)
      })
    })

    output$indiv_table <- renderDT({
      req(rvs$indiv_tbl)
      datatable(rvs$indiv_tbl, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    })

    output$dl_indiv <- downloadHandler(
      filename = function() paste0("despace_individual_test_", input$domain_pick %||% "domain", "_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv"),
      content = function(file) write.csv(rvs$indiv_tbl %||% data.frame(), file, row.names = FALSE)
    )

    #--------------------------
    # Quick gene view across domains
    #--------------------------
    output$gene_violin <- renderPlot({
      req(rvs$spe2, rvs$cluster_col, input$gene_pick)
      .make_gene_violin(rvs$spe2, rvs$cluster_col, input$gene_pick)
    })
  })
}
