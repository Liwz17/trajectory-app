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
