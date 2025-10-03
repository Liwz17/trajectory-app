# R/03_module_select.R
# 依赖：plotly, base64enc, png
# 功能：
#  - 把 HE 背景作为底图（支持传入旋转后的图像和宽高）
#  - 在其上绘制 spot 点
#  - 用户用 lasso/box 选点，返回选中点的 data.frame（含 barcode/x_hires/y_hires/expr）
#  - 提供清空选择、下载 CSV

selectSpotsUI <- function(id, title = "在 HE 图上选择区域") {
  ns <- NS(id)
  tagList(
    tags$h4(title),
    radioButtons(ns("mode"), "选择工具", choices = c("套索 (lasso)" = "lasso", "矩形 (box)" = "select"), inline = TRUE),
    plotly::plotlyOutput(ns("p"), height = "750px"),
    fluidRow(
      column(4, actionButton(ns("clear"), "清空选择")),
      column(4, downloadButton(ns("dl"), "下载所选 CSV")),
      column(4, verbatimTextOutput(ns("nsel")))
    ),
    tags$small("提示：使用图上工具栏的套索/矩形工具进行选择；按住 Shift 可多次累积选择（plotly 行为）。")
  )
}

# 参数：
#  - pos_reactive(): reactive，返回 data.frame，至少包含：barcode, x_hires, y_hires
#  - expr_reactive(): reactive，返回 numeric 向量，对应 pos 中 rows 的表达值（已 log1p 与否由外部决定）
#  - img_reactive():  reactive，返回 list(img = raster, img_w = W, img_h = H)
# 返回：
#  list(
#    selected = reactive(data.frame(...)),  # 当前选择的 spot 表
#    plot_proxy = reactiveVal/NULL          # 如需对 plot 做后续联动可扩展
#  )
selectSpotsServer <- function(id, pos_reactive, expr_reactive, img_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    requireNamespace("plotly")
    requireNamespace("base64enc")
    requireNamespace("png")

    # 把 raster 背景转为 data URI，供 plotly layout(images=...) 使用
    raster_to_data_uri <- function(rs) {
      stopifnot(inherits(rs, "raster"))
      # raster -> array (H x W x 3/4), 再写入 PNG 原始字节流
      # as.raster -> matrix of "#RRGGBB" or with alpha; 我们直接画在 png 里
      m <- grDevices::col2rgb(rs, alpha = TRUE)
      # m shape: 4 x (H*W)，还原成 array(H,W,4)
      # 注意：as.raster 的列主序，需要我们知道 dim；无法直接从 raster取 dim，只能从外部 img_reactive 提供的 w/h
      # 因此这里不从 rs 还原，改用 magick 方案更稳；但为避免 magick 依赖，这里用一个简化路径：
      # 直接把 raster 绘制到临时 PNG
      tf <- tempfile(fileext = ".png")
      # 用 png() 画底图
      # 读取宽高由上层传入
      return(tf) # 我们换一个更稳妥方案：用 file 协议，不走 data URI
    }

    # 我们采用更稳的方式：把背景图写临时 PNG 文件，然后通过 file:// 路径给 plotly 用
    # 这里写一个 util：把 raster 写为临时 PNG，返回文件路径
    write_raster_to_png <- function(rs, width, height) {
      tf <- tempfile(fileext = ".png")
      # 用 png::writePNG 需要 0..1 array；raster->数值数组较繁琐
      # 更简单：用 magick 渲染 raster 再写文件（magick 在你的项目里已存在）
      im <- magick::image_read(rs)
      im <- magick::image_resize(im, paste0(width, "x", height, "!"))  # 强制到 width x height
      magick::image_write(im, path = tf, format = "png")
      tf
    }

    # 当前图数据：依赖上游 pos/expr/img（img 由旋转模块传入）
    plot_data <- reactive({
      pos <- pos_reactive(); req(pos)
      expr <- expr_reactive(); req(expr)
      stopifnot(nrow(pos) == length(expr))

      img_obj <- img_reactive(); req(img_obj$img, img_obj$img_w, img_obj$img_h)
      # 将 raster 写为临时 PNG，供 plotly 背景图使用
      png_path <- write_raster_to_png(img_obj$img, img_obj$img_w, img_obj$img_h)

      list(
        pos = pos,
        expr = expr,
        img_w = img_obj$img_w,
        img_h = img_obj$img_h,
        png_path = png_path
      )
    })

    # 画交互图
    output$p <- plotly::renderPlotly({
      d <- plot_data()
      pos <- d$pos
      expr <- d$expr
      W <- d$img_w; H <- d$img_h

      # plotly 的 y 轴默认向上增，为了和 ggplot 一致（图像坐标原点左上），我们反转 y 轴
      p <- plotly::plot_ly(
        type = "scattergl",
        mode = "markers",
        x = pos$x_hires,
        y = pos$y_hires,
        marker = list(size = 6, opacity = 0.9, color = expr, colorscale = "Viridis", colorbar = list(title = "expr")),
        hoverinfo = "text",
        text = paste0("barcode: ", pos$barcode,
                      "<br>x: ", round(pos$x_hires,1),
                      "<br>y: ", round(pos$y_hires,1),
                      "<br>expr: ", signif(expr, 4))
      ) %>%
        plotly::layout(
          xaxis = list(range = c(0, W), zeroline = FALSE, showgrid = FALSE),
          yaxis = list(range = c(H, 0), zeroline = FALSE, showgrid = FALSE, scaleanchor = "x", scaleratio = 1),
          images = list(
            list(
              source = d$png_path,  # plotly 支持本地文件路径（在 Shiny 下会被静态服用）
              xref = "x", yref = "y",
              x = 0, y = 0,
              sizex = W, sizey = H,
              sizing = "stretch",
              layer = "below"
            )
          ),
          dragmode = if (input$mode == "lasso") "lasso" else "select",
          margin = list(l = 0, r = 0, t = 0, b = 0),
          showlegend = FALSE
        )

      p
    })

    # 跟随切换选择模式（lasso/select）
    observeEvent(input$mode, {
      plotly::plotlyProxy(ns("p"), session) %>%
        plotly::plotlyProxyInvoke("relayout", list(dragmode = if (input$mode == "lasso") "lasso" else "select"))
    })

    # 取选中的点（event_data）
    # plotly_selected：完成选择时触发；返回被选点在 trace 中的 pointNumbers（从 0 开始）
    selected_idx <- reactiveVal(integer(0))

    observeEvent(plotly::event_data("plotly_selected", source = NULL), {
      ed <- plotly::event_data("plotly_selected")
      if (is.null(ed) || nrow(ed) == 0) {
        selected_idx(integer(0))
      } else {
        # 累积选择：合并已选与新选
        idx <- ed$pointNumber
        selected_idx(sort(unique(c(selected_idx(), idx))))
      }
    })

    # 清空
    observeEvent(input$clear, {
      selected_idx(integer(0))
      # 清除高亮（重绘或调用 restyle）
      plotly::plotlyProxy(ns("p"), session) %>%
        plotly::plotlyProxyInvoke("restyle", list(selectedpoints = list(NULL)))
    })

    # 导出所选
    output$dl <- downloadHandler(
      filename = function() {
        paste0("selected_spots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        sel <- sel_table()
        utils::write.csv(sel, file, row.names = FALSE)
      }
    )

    # 汇总表
    sel_table <- reactive({
      d <- plot_data()
      idx <- selected_idx()
      if (length(idx) == 0) {
        return(d$pos[0, c("barcode","x_hires","y_hires")])
      }
      pos <- d$pos
      expr <- d$expr
      # plotly pointNumber 从 0 开始，R 索引从 1 开始
      idx_r <- idx + 1L
      out <- data.frame(
        barcode = pos$barcode[idx_r],
        x_hires = pos$x_hires[idx_r],
        y_hires = pos$y_hires[idx_r],
        expr    = expr[idx_r],
        stringsAsFactors = FALSE
      )
      out
    })

    output$nsel <- renderText({
      paste0("已选 spots: ", nrow(sel_table()))
    })

    # 对外暴露：所选数据
    return(list(
      selected = sel_table
    ))
  })
}

setup_selection <- function(output, session, selected_idx, sel_mode_input, clear_btn_input) {
  attach_plotly_selection_handlers(
    output_id       = "he_plotly",
    session         = session,
    selected_idx    = selected_idx,
    sel_mode_input  = sel_mode_input,
    clear_btn_input = clear_btn_input
  )
  output$sel_count <- renderText(paste0("已选 spots: ", length(selected_idx())))
}




init_annotations <- function() {
  reactiveVal(data.frame(barcode = character(), label = character(), stringsAsFactors = FALSE))
}

setup_annotation_observers <- function(input, output, rv, selected_idx, annotations) {
  # 添加
  observeEvent(input$add_annotation, {
    req(rv$pos, length(selected_idx()) > 0, input$bio_label)
    bc <- rv$pos$barcode[selected_idx()]
    annotations( append_annotations(annotations(), bc, input$bio_label) )
  })

  # 清空
  observeEvent(input$clear_annotations, {
    annotations(data.frame(barcode = character(), label = character(), stringsAsFactors = FALSE))
  })

  # 展示
  output$anno_summary <- renderTable({
    summarize_annotations(annotations())
  }, striped = TRUE, bordered = TRUE, hover = TRUE)
}

# 转换成 label->索引 list
make_label_index_list <- function(rv, annotations) {
  req(rv$pos)
  ann <- annotations()
  if (nrow(ann) == 0) return(list())
  m <- match(ann$barcode, rv$pos$barcode)
  keep <- !is.na(m)
  if (!any(keep)) return(list())
  ann$idx <- m
  split(ann$idx, ann$label)
}
