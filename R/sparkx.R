sparkxServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$run_sparkx, {
      req(rv$counts, rv$pos)
      counts <- isolate(rv$counts)
      coords <- as.matrix(isolate(rv$pos[, c("col","row")]))
      heg_n  <- 1000
      cores  <- max(1, as.integer(input$sparkx_cores))
      withProgress(message="Running SPARK-X...", value=0.1, {
        fut <- future({
          gene_score <- Matrix::rowSums(counts)
          keep_idx <- order(gene_score, decreasing=TRUE)[seq_len(min(heg_n, length(gene_score)))]
          counts_f <- as.matrix(counts[keep_idx, , drop=FALSE]); storage.mode(counts_f) <- "double"
          sparkx_formals <- names(formals(SPARK::sparkx))
          res <- if ("numCores" %in% sparkx_formals) SPARK::sparkx(counts_f, coords, numCores=cores, option="mixture")
                 else SPARK::sparkx(counts_f, coords, option="mixture")
          df <- tryCatch(as.data.frame(res$res_mtest), error=function(e) as.data.frame(res))
          nms <- tolower(names(df))
          if (!"gene" %in% nms) df$gene <- rownames(df)
          cand_p <- intersect(nms, c("combinedpval","pvalue","p_val","combined_pvalue","pval"))
          if (length(cand_p)) names(df)[match(cand_p[1], nms)] <- "pval"
          nms <- tolower(names(df))
          if (!"qval" %in% nms) {
            if ("adjustedpval" %in% nms) names(df)[match("adjustedpval", nms)] <- "qval"
            else if ("pval" %in% tolower(names(df))) df$qval <- p.adjust(df$pval, method="BH")
          }
          df$gene_total_counts <- gene_score[match(df$gene, names(gene_score))]
          df
        })
        fut %...>% (function(df){
          rv$sparkx_res <- df
          output$sparkx_table <- DT::renderDT({
            cut <- as.numeric(input$sparkx_fdr)
            df_show <- df
            if ("qval" %in% names(df_show)) {
              df_show <- df_show[order(df_show$qval, na.last=NA), ]
              if (!is.na(cut)) df_show <- subset(df_show, qval <= cut)
            }
            DT::datatable(head(df_show, 5000), options=list(pageLength=25), rownames=FALSE)
          })
          output$dl_sparkx <- downloadHandler(
            filename = function() sprintf("sparkx_results_HEG%d_fdr%s.csv", heg_n, input$sparkx_fdr),
            content  = function(file) write.csv(rv$sparkx_res, file, row.names=FALSE)
          )
          showNotification(sprintf("SPARK-X finished (HEG=%d).", heg_n))
        }) %...!% (function(e){
          showNotification(paste("SPARK-X failed:", e$message), type="error", duration=NULL)
        })
      })
    })
  })
}


# ---- sparkx UI 模块 ----
sparkxUI <- function(id) {
  ns <- NS(id)
  tagList(
    actionButton(ns("run_sparkx"), "Run SPARK-X", class = "btn-primary"),
    numericInput(ns("sparkx_cores"), "Cores",
                 value = max(1, parallel::detectCores()-1), min = 1),
    numericInput(ns("sparkx_fdr"), "FDR cutoff",
                 value = 0.05, min = 1e-6, step = 0.01),
    downloadButton(ns("dl_sparkx"), "Download SPARK-X results"),
    DTOutput(ns("sparkx_table"))
  )
}
