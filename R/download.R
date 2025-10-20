# R/05_downloads.R
setup_download_selected <- function(output, id, pos_r, expr_vec_r, selected_idx_r) {
  output[[id]] <- downloadHandler(
    filename = function() paste0("selected_spots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      pos <- pos_r(); idx <- selected_idx_r()
      if (length(idx)==0) {
        utils::write.csv(pos[0, c("barcode","x_hires","y_hires")], file, row.names = FALSE)
      } else {
        v <- expr_vec_r()
        out <- data.frame(
          barcode = pos$barcode[idx],
          x_hires = pos$x_hires[idx],
          y_hires = pos$y_hires[idx],
          expr    = v[idx],
          stringsAsFactors = FALSE
        )
        utils::write.csv(out, file, row.names = FALSE)
      }
    }
  )
}


# exportToolsUI_global <- function() {
#   tagList(
#     hr(), h4("Export & tools"),
#     fluidRow(
#       column(3, actionButton("clear_sel", "Clear selection")),
#       column(3, downloadButton("dl_sel", "Download selected CSV")),
#       column(3, verbatimTextOutput("sel_count")),
#       column(3, downloadButton("dl_multi_distance", "Download all structure distances")),
#       column(3, downloadButton(ns("dl_gam_results"), "Download GAM results"))
#     )
#   )
# }

exportToolsUI_global <- function() {
  tagList(
    hr(), h4("Export & tools"),
    fluidRow(
      column(3, actionButton("clear_sel", "Clear selection")),
      # column(3, downloadButton("dl_sel", "Download selected CSV")),
      column(3, verbatimTextOutput("sel_count")),
      column(3, downloadButton("dl_multi_distance", "Download all structure distances"))
    ),
    fluidRow(
      column(4, downloadButton(NS("gam", "dl_gam_results"), "Download GAM fitting results"))
    )
  )
}

