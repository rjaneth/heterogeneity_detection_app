# app.R
# SAR Heterogeneity Detection
req_pkgs <- c("shiny", "viridisLite", "fields", "png")
for (pkg in req_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}


lapply(req_pkgs, library, character.only = TRUE)
# library(shiny)
# library(viridisLite)
# library(fields)          # for image.plot legend
# library(png)             # for readPNG
# library(shinycssloaders) # for the spinner


source("./Code/al_omari_1_estimator.R")
source("./Code/bootstrap_al_omari_1_estimator.R")
source("./Code/renyi_entropy_estimator_v1.R")
source("./Code/bootstrap_renyi_entropy_estimator_v1.R")
source("./Code/tsallis_estimator_optimized.R")
source("./Code/bootstrap_tsallis_entropy_optimized.R")
source("./Code/read_ENVI_images.R")
source("./Code/imagematrix_visualizer.R")


shannon_theoretical <- function(L, mu) {
  log(mu) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L))
}
renyi_theoretical <- function(L, λ, mu) {
  num <- λ * lgamma(L) - lgamma(λ * (L - 1) + 1) + (λ * (L - 1) + 1) * log(λ)
  (num/(λ-1)) + log(mu) - log(L)
}
tsallis_theoretical <- function(L, λ, mu) {
  (1 - exp((1 - λ)*log(mu) + (λ - 1)*log(L) +
             lgamma(λ*(L-1) + 1) - λ*lgamma(L) -
             (λ*(L-1) + 1)*log(λ))) / (λ - 1)
}
calc_pvals <- function(stat_mat) {
  μ   <- mean(stat_mat, na.rm = TRUE)
  σ   <- sd(stat_mat,   na.rm = TRUE)
  eps <- stat_mat / σ
  2 * pnorm(-abs(eps))
}

# --- UI ---------------------------------------------------------------
ui <- navbarPage("SAR Heterogeneity Detection",
                 
                 # App tab
                 tabPanel("App",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("input_type", "Data source:",
                                          choices = c("Simulated Images Examples" = "sim",
                                                      "ENVI Images Examples"       = "envi",
                                                      "Upload ENVI"       = "upload")
                              ),
                              conditionalPanel(
                                "input.input_type == 'sim'",
                                selectInput("sim_choice", "Choose simulated image:",
                                            choices = list(
                                              "Simulated image (L = 5)" = "Phantom_4_z.Rdata",
                                              "Simulated image (L = 9)" = "Phantom4_L9_1_8.Rdata"
                                            ))
                              ),
                              conditionalPanel(
                                "input.input_type == 'envi'",
                                selectInput("envi_choice", "Choose ENVI sample:",
                                            choices = list(
                                              "Dublin L16 (600×600)"      = "L16_envi_dublin_size_600",
                                              "London L1 (1000×1000)"     = "L1_envi_london_size_1000",
                                              "New Orleans L12 (600×600)" = "L12_envi_New_Orleans_size_600"
                                            )
                                ),
                                textInput("envi_band", "Header basename:", value = "Intensity_HH")
                              ),
                              conditionalPanel(
                                "input.input_type == 'upload'",
                                fileInput("upload_img", "Upload .img file", accept = ".img"),
                                fileInput("upload_hdr", "Upload .hdr file", accept = ".hdr")
                              ),
                              selectInput("entropy", "Entropy type:",
                                          choices = c("Shannon" = "shannon",
                                                      "Rényi"   = "renyi",
                                                      "Tsallis" = "tsallis")
                              ),
                              conditionalPanel(
                                "input.entropy != 'shannon'",
                                sliderInput("lambda", "Order λ:", min = 0.1, max = 3,
                                            value = 0.9, step = 0.05)
                              ),
                              numericInput("looks", "Number of looks L:", value = 5, min = 1),
                              checkboxInput("bootstrap", "Enable bootstrap", value = TRUE),
                              conditionalPanel(
                                "input.bootstrap",
                                numericInput("B", "Bootstrap replicates:", value = 5, min = 1)
                              ),
                              numericInput("window", "Window side (pixels):",
                                           value = 7, min = 3, step = 2),
                              actionButton("go", "Run Detection", class = "btn-primary"),
                              width = 3
                            ),
                            mainPanel(
                              #withSpinner(plotOutput("pvalPlot", height = "600px")),
                              plotOutput("pvalPlot", height = "600px"),
                              verbatimTextOutput("timing"),
                              width = 9
                            )
                          )
                 ),
                 
                 # 
                 tabPanel("About",
                          fluidRow(
                            column(8,
                                   h4("About this app"),
                                   p("This interactive Shiny application allows users to detect heterogeneity in SAR images."),
                                   tags$ul(
                                     tags$li("Load predefined simulated images, ENVI image examples, or upload your own .img/.hdr files."),
                                     tags$li("Select entropy estimator: Shannon, Rényi, or Tsallis."),
                                     tags$li("Optional bootstrap resampling and sliding-window analysis (n×n pixels)."),
                                     tags$li("Visualize resulting p-value map in color with interactive scale bar.")
                                   ),
                                   h5("Recommendation"),
                                   p("⚙️ On the hosted version, heavy bootstrap (B > 10) can be slow. We suggest starting with 5–10 replicates for interactive use."),
                                   h5("Note on performance"),
                                   p("⚠️ This app is hosted on the free tier of Shinyapps.io, which has limited computing resources (25 active hours/month)."),
                                   p("For unlimited use and faster performance, we recommend running the app or the standalone script locally from the GitHub repository:"),
                                   tags$ul(
                                     tags$li(a("GitHub repository with app and script", 
                                               href="https://github.com/rjaneth/heterogeneity_detection_app",target="_blank")),
                                     tags$li("The same repository also contains a R script (`heterogeneity_detection_v1.R`) designed for large images with parallel processing (no Shiny interface).")
                                   ),
                                   p("To run locally, clone the repository, open `app.R` in RStudio, and click “Run App.”"),
                                   p("To process large SAR images, edit the `opt` parameters in `heterogeneity_detection_v1.R` and run it directly.")
                            )
                          )
                 )
)  # 
# --- SERVER 
server <- function(input, output, session) {
  
  observeEvent(input$go, {
    
    #
    output$pvalPlot <- renderPlot({ plot.new() })
    
    # 
    if (input$input_type == "sim") {
      load(file.path("www", input$sim_choice))  # 
      img_mat <- Z
    } else if (input$input_type == "envi") {
      folder   <- file.path("www", input$envi_choice)
      img_file <- file.path(folder, paste0(input$envi_band, ".img"))
      hdr_file <- file.path(folder, paste0(input$envi_band, ".hdr"))
      img_mat  <- myread.ENVI(img_file, headerfile = hdr_file)
    } else {
      req(input$upload_img, input$upload_hdr)
      img_mat <- myread.ENVI(input$upload_img$datapath,
                             headerfile = input$upload_hdr$datapath)
    }
    
    # 2) Select entropy estimator 
    estimator <- switch(input$entropy,
                        shannon = {
                          if (input$bootstrap && input$looks > 1) {
                            function(z) bootstrap_al_omari_1_estimator(z, input$B)
                          } else {
                            function(z) al_omari_1_estimator(z)
                          }
                        },
                        renyi = {
                          if (input$bootstrap && input$looks > 1) {
                            function(z) bootstrap_renyi_entropy_estimator_v1(
                              z, input$B, input$lambda)
                          } else {
                            function(z) renyi_entropy_estimator_v1(z, input$lambda)
                          }
                        },
                        tsallis = {
                          if (input$bootstrap && input$looks > 1) {
                            function(z) bootstrap_tsallis_entropy_optimized(
                              z, input$B, input$lambda)
                          } else {
                            function(z) tsallis_estimator_optimized(z, input$lambda)
                          }
                        }
    )
    
    # 3) Theoretical function
    theoretical_fn <- switch(input$entropy,
                             shannon = function(L, λ, mu) shannon_theoretical(L, mu),
                             renyi   = renyi_theoretical,
                             tsallis = tsallis_theoretical
    )
    
    # 4) Sliding-window processing with progress per row
    start_time <- Sys.time()
    nr <- nrow(img_mat) - input$window + 1
    nc <- ncol(img_mat) - input$window + 1
    stat_mat <- matrix(NA_real_, nr, nc)
    
    withProgress(message = "Running detection...", value = 0, {
      for (i in seq_len(nr)) {
        for (j in seq_len(nc)) {
          w            <- img_mat[i:(i+input$window-1),
                                  j:(j+input$window-1)]
          stat_mat[i,j]<- estimator(w) -
            theoretical_fn(input$looks,
                           input$lambda,
                           mean(w))
        }
        incProgress(1 / nr)
      }
    })
    
    # 5) Compute p-values & format runtime mm:ss
    pvals    <- calc_pvals(stat_mat)
    end_time <- Sys.time()
    elapsed  <- as.numeric(difftime(end_time, start_time, units = "secs"))
    mins     <- floor(elapsed / 60)
    secs     <- round(elapsed %% 60)
    
    output$timing <- renderText({
      sprintf("Run time: %02d:%02d (mm:ss)", mins, secs)
    })
    
    # 6) Render the colored p-value map
    output$pvalPlot <- renderPlot({
      previewImagematrixPanel(
        imagematrix_color(pvals),
        palette_option     = "viridis-H",
        significance_level = NULL
      )
    })
  })
}

# Launch the app
shinyApp(ui, server)
