# R/02_module_rotate.R
# Dependency: magick
# Install: install.packages("magick"); if it fails, install system ImageMagick first

rotateImageUI <- function(id, title = "HE Image Rotation") {
  ns <- NS(id)
  tagList(
    tags$h4(title),
    sliderInput(ns("angle"), "Rotation angle (°)", min = -180, max = 180, value = 0, step = 1),
    fluidRow(
      column(4, actionButton(ns("ccw90"), "↺ 90°")),
      column(4, actionButton(ns("reset"), "Reset")),
      column(4, actionButton(ns("cw90"),  "↻ 90°"))
    ),
    fluidRow(
      column(6, actionButton(ns("flip_h"), "Flip horizontally")),
      column(6, actionButton(ns("flip_v"), "Flip vertically"))
    ),
    helpText("Only rotate or flip the background HE image; spot coordinates remain unchanged.")
  )
}


rotateImageServer <- function(id, image_path_reactive) {
  moduleServer(id, function(input, output, session) {
    requireNamespace("magick")

    # Internal state: angle and flip
    angle <- reactiveVal(0)
    flip_h <- reactiveVal(FALSE)
    flip_v <- reactiveVal(FALSE)

    # Rotation angle controls
    observeEvent(input$angle,  ignoreInit = TRUE, { angle(input$angle) })
    observeEvent(input$ccw90,                     { angle((angle() - 90) %% 360) })
    observeEvent(input$cw90,                      { angle((angle() + 90) %% 360) })
    observeEvent(input$reset, {
      angle(0); flip_h(FALSE); flip_v(FALSE)
    })

    # Flip controls
    observeEvent(input$flip_h, { flip_h(!flip_h()) })
    observeEvent(input$flip_v, { flip_v(!flip_v()) })

    # Output: rotated + flipped raster and width/height
    rotated <- reactive({
      path <- image_path_reactive()
      req(path)

      im <- magick::image_read(path)
      # Rotation
      if (angle() != 0) {
        im <- magick::image_rotate(im, degrees = angle())
      }
      # Flip
      if (isTRUE(flip_h())) im <- magick::image_flop(im)
      if (isTRUE(flip_v())) im <- magick::image_flip(im)

      rs <- as.raster(im)
      inf <- magick::image_info(im)
      list(img = rs, img_w = inf$width[1], img_h = inf$height[1],
           angle = angle(), flip_h = flip_h(), flip_v = flip_v())
    })

    # Expose getter to outside
    list(get = rotated,
         angle = reactive(angle()),
         flip_h = reactive(flip_h()),
         flip_v = reactive(flip_v()))
  })
}
