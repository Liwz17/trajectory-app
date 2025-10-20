# R/04_plot_he.R
build_he_plotly_figure <- function(pos, expr, img_obj, ptsize, alpha, gene_label) {
  stopifnot(!is.null(pos), !is.null(img_obj))
  W <- img_obj$img_w; H <- img_obj$img_h

  df <- data.frame(
    idx = seq_len(nrow(pos)),
    barcode = pos$barcode,
    x = pos$x_hires, y = pos$y_hires, expr = expr, stringsAsFactors = FALSE
  )
  bg_uri <- raster_to_data_uri(img_obj$img, W, H)

  plotly::plot_ly(
    df, x = ~x, y = ~y, type = "scattergl", mode = "markers",
    marker = list(size = ptsize, opacity = alpha, color = ~expr,
                  colorscale = "Viridis",
                  colorbar = list(title = gene_label %||% "expression")),
    text = ~paste0("barcode: ", barcode,
                   "<br>x: ", round(x,1),
                   "<br>y: ", round(y,1),
                   "<br>expr: ", signif(expr,4)),
    hoverinfo = "text",
    source = "main",
    customdata = ~idx
  ) %>%
    plotly::layout(
      xaxis = list(range = c(0, W), zeroline = FALSE, showgrid = FALSE),
      yaxis = list(range = c(H, 0), zeroline = FALSE, showgrid = FALSE,
                   scaleanchor = "x", scaleratio = 1),
      images = list(list(source = bg_uri, xref="x", yref="y",
                         x=0, y=0, sizex=W, sizey=H,
                         sizing="stretch", layer="below")),
      showlegend = FALSE
    )
}



# R/06_select_plotly.R
attach_plotly_selection_handlers <- function(output_id, session, selected_idx, sel_mode_input, clear_btn_input) {
  observeEvent(plotly::event_data("plotly_selected", source = "main"), {
    ed <- plotly::event_data("plotly_selected", source = "main")
    if (is.null(ed) || nrow(ed) == 0) selected_idx(integer(0)) else {
      idx <- ed$customdata
      if (is.null(idx)) idx <- ed$pointNumber + 1L
      selected_idx(sort(unique(as.integer(idx))))
    }
  })

  observeEvent(sel_mode_input(), {
    plotly::plotlyProxy(output_id, session) %>%
      plotly::plotlyProxyInvoke("relayout",
        list(dragmode = if (sel_mode_input()=="select") "select" else "lasso"))
  })

  observeEvent(clear_btn_input(), {
    selected_idx(integer(0))
    plotly::plotlyProxy(output_id, session) %>%
      plotly::plotlyProxyInvoke("restyle", list(selectedpoints = list(NULL)))
  })
}


get_current_image <- function(rot, rv) {
  ro <- tryCatch(rot$get(), error=function(e) NULL)
  if (!is.null(ro)) ro else if (!is.null(rv$img)) list(img=rv$img, img_w=rv$img_w, img_h=rv$img_h) else NULL
}
