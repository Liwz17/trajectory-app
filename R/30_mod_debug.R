##############################
# mod_debug.R
##############################
# Lightweight debug panel (unchanged info, just moved here)

mod_debug_ui <- function(id) {
  ns <- NS(id)
  verbatimTextOutput(ns("dbg"))
}

mod_debug_server <- function(id, visium_state) {
  moduleServer(id, function(input, output, session) {
    output$dbg <- renderPrint({
      list(
        img_dim      = dim(visium_state$img_arr),
        hires_scale  = visium_state$hires_scale,
        range_x      = range(visium_state$pos$x_hires, na.rm = TRUE),
        range_y      = range(visium_state$pos$y_hires, na.rm = TRUE),
        head_coords  = head(visium_state$pos[, c("barcode","x_hires","y_hires")]),
        img_summary  = list(min_val = min(visium_state$img_arr, na.rm = TRUE),
                            max_val = max(visium_state$img_arr, na.rm = TRUE))
      )
    })
  })
}