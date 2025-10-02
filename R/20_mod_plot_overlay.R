##############################
# mod_plot_overlay.R
##############################
# HE background + points overlay; keeps your original plotting logic

mod_plot_overlay_ui <- function(id) {
  ns <- NS(id)
  plotOutput(ns("he_plot"), height = "750px", width = "100%")
}

mod_plot_overlay_server <- function(id, visium_state, ctrl) {
  moduleServer(id, function(input, output, session) {
    counts <- visium_state$counts
    pos    <- visium_state$pos
    img_df <- visium_state$img_df
    img_w  <- visium_state$img_dims$w
    img_h  <- visium_state$img_dims$h
    hires_scale <- visium_state$hires_scale

    expr_vec <- reactive({
      g <- as.character(ctrl$gene())
      v <- if (g %in% rownames(counts)) as.numeric(counts[g, pos$barcode]) else rep(0, nrow(pos))
      if (isTRUE(ctrl$log1p())) v <- log1p(v)
      v
    })

    output$he_plot <- renderPlot({
      v  <- expr_vec()
      df <- data.frame(x = pos$x_hires, y = pos$y_hires, expr = v)
      legend_name <- as.character(ctrl$gene() %||% "expression")

      ggplot(img_df, aes(x, y)) +
        geom_raster(aes(fill = col)) +
        scale_fill_identity() +
        geom_point(data = df, aes(x = x, y = y, color = expr), size = ctrl$ptsize(), alpha = ctrl$alpha()) +
        scale_color_viridis_c(name = legend_name) +
        coord_fixed(xlim = c(0, img_w), ylim = c(img_h, 0), expand = FALSE) +
        theme_void() +
        labs(title = paste0("HE overlay @ hires (scale=", round(hires_scale, 3), ")"))
    }, res = 120)
  })
}
