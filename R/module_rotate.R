# R/02_module_rotate.R
# 依赖：magick
# 安装：install.packages("magick")；若失败先装系统 ImageMagick

rotateImageUI <- function(id, title = "HE 图像旋转") {
  ns <- NS(id)
  tagList(
    tags$h4(title),
    sliderInput(ns("angle"), "旋转角度 (°)", min = -180, max = 180, value = 0, step = 1),
    fluidRow(
      column(4, actionButton(ns("ccw90"), "↺ 90°")),
      column(4, actionButton(ns("reset"), "重置")),
      column(4, actionButton(ns("cw90"),  "↻ 90°"))
    ),
    fluidRow(
      column(6, actionButton(ns("flip_h"), "水平翻转")),
      column(6, actionButton(ns("flip_v"), "垂直翻转"))
    ),
    helpText("只旋转或翻转背景 HE 图像，点坐标保持不变。")
  )
}


rotateImageServer <- function(id, image_path_reactive) {
  moduleServer(id, function(input, output, session) {
    requireNamespace("magick")

    # 内部维护角度和翻转状态
    angle <- reactiveVal(0)
    flip_h <- reactiveVal(FALSE)
    flip_v <- reactiveVal(FALSE)

    # 旋转角度控制
    observeEvent(input$angle,  ignoreInit = TRUE, { angle(input$angle) })
    observeEvent(input$ccw90,                     { angle((angle() - 90) %% 360) })
    observeEvent(input$cw90,                      { angle((angle() + 90) %% 360) })
    observeEvent(input$reset, {
      angle(0); flip_h(FALSE); flip_v(FALSE)
    })

    # 翻转控制
    observeEvent(input$flip_h, { flip_h(!flip_h()) })
    observeEvent(input$flip_v, { flip_v(!flip_v()) })

    # 产出：旋转+翻转后的 raster 与宽高
    rotated <- reactive({
      path <- image_path_reactive()
      req(path)

      im <- magick::image_read(path)
      # 旋转
      if (angle() != 0) {
        im <- magick::image_rotate(im, degrees = angle())
      }
      # 翻转
      if (isTRUE(flip_h())) im <- magick::image_flop(im)
      if (isTRUE(flip_v())) im <- magick::image_flip(im)

      rs <- as.raster(im)
      inf <- magick::image_info(im)
      list(img = rs, img_w = inf$width[1], img_h = inf$height[1],
           angle = angle(), flip_h = flip_h(), flip_v = flip_v())
    })

    # 对外暴露 getter
    list(get = rotated,
         angle = reactive(angle()),
         flip_h = reactive(flip_h()),
         flip_v = reactive(flip_v()))
  })
}
