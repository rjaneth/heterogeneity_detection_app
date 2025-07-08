# ------------------------------------------------------------
# heterogeneity_detection_v1.R  (v1.)
# ------------------------------------------------------------
# Author: Janeth Alpala
# Date: 30-jun-2025
# Description:
#   Script for heterogeneity detection in SAR images (ENVI format)
#   or simulated (RData), using non-parametric estimators
#   of Shannon, Rényi, or Tsallis entropy (with or without bootstrap).
#
# Quick-start:
#   1. Edit the `opt` list below with your own paths, entropy type,
#      window size, etc.
#   2. Run in R: source("heterogeneity_detection_v1.R")
#   3. At the end, a color p-value map with scale bar is displayed.
# Outputs:
#   • On-screen plot of color p-values.
#   • PNG saved in ./PNG/ for p-values.
#   • .RData saved in ./Data/ containing img_mat, entropy_test, p_values, and opt.

# ----------------------------------------------------------------------
# Required packages
# ----------------------------------------------------------------------
req_pkgs <- c("viridisLite", "png", "progress", "fields")
for (p in req_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(viridisLite)
library(png)
library(progress)
library(fields)   # for image.plot legend

# ----------------------------------------------------------------------
# Optional: Parallelization setup 
# ----------------------------------------------------------------------
if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
library(future.apply)
plan(multisession)   # Automatically uses available CPU cores 
parallel <- TRUE     # Set to FALSE to disable parallel execution
set.seed(1234567890, kind = "Mersenne-Twister")  

# ----------------------------------------------------------------------
# External functions
# ----------------------------------------------------------------------
source("./Code/al_omari_1_estimator.R")
source("./Code/bootstrap_al_omari_1_estimator.R")
source("./Code/renyi_entropy_estimator_v1.R")
source("./Code/bootstrap_renyi_entropy_estimator_v1.R")
source("./Code/tsallis_estimator_optimized.R")
source("./Code/bootstrap_tsallis_entropy_optimized.R")
source("./Code/read_ENVI_images.R")

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
shannon_theoretical <- function(L, mu) {
  log(mu) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L))
}

renyi_theoretical <- function(L, lambda, mu) {
  num <- lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
    (lambda * (L - 1) + 1) * log(lambda)
  (num / (lambda - 1)) + log(mu) - log(L)
}

tsallis_theoretical <- function(L, lambda, mu) {
  (1 - exp((1 - lambda)*log(mu) +
             (lambda - 1)*log(L) +
             lgamma(lambda*(L - 1) + 1) -
             lambda*lgamma(L) -
             (lambda*(L - 1) + 1)*log(lambda))) / (lambda - 1)
}

calculate_p_values_matrix <- function(stat_mat) {
  mu  <- mean(stat_mat, na.rm = TRUE)
  sig <- sd(stat_mat,   na.rm = TRUE)
  eps <- (stat_mat - 0) / sig
  2 * pnorm(-abs(eps))
}

# ----------------------------------------------------------------------
# Manual parameters
# ----------------------------------------------------------------------
# examples: SAR data in envi format is located in the Data folder.
# image        = "./Data/SAR/L12_envi_munich_size_1024/Intensity_HH.img",
# header       = "./Data/SAR/L12_envi_munich_size_1024/Intensity_HH.hdr",

# sim_rdata    = "./Data/L5_simulated_image_500.Rdata",

opt <- list(
  # For ENVI inputs: 
  image        = "./Data/SAR/L16_envi_dublin_size_600/Intensity_HH.img",
  header       = "./Data/SAR/L16_envi_dublin_size_600/Intensity_HH.hdr",
  # For simulated inputs:
  sim_rdata    = "./Data/L9_simulated_image_500.Rdata",
  # Input type selector:
  #   - "envi": load SAR images in ENVI format (.img + .hdr)
  #   - "sim" : load simulated images from an .RData file
  input_type   = "envi",         # choose either "envi" or "sim"
  # Entropy settings:
  entropy      = "tsallis",      # "shannon" | "renyi" | "tsallis"
  lambda       = 0.85,           # order parameter λ for Rényi & Tsallis (for L=1, choose Renyi (λ=3) and for Tsallis (λ=1.1)
  # Bootstrap settings:
  bootstrap    = TRUE,           # TRUE to enable bootstrap resampling or FALSE to disable
  B            = 50,             # number of bootstrap replicates
  # SAR image parameters:
  looks        = 16,             # number of looks L
  window       = 5,              # window side length (5x5,7×7,9x9... pixels)
  # Output & display:
  prefix       = "L16_dublin_600_tsallis_w9_b50", # Prefix for saved files, e.g., "image_size_sar.png"
  progress_step= 10,             # progress bar refresh interval
  p_thr        = 0.05            # significance threshold
)

# ----------------------------------------------------------------------
# Basic checks
# ----------------------------------------------------------------------
stopifnot(opt$input_type %in% c("envi", "sim"))
if (opt$window < 1) stop("Window size must be a positive integer")
if (opt$entropy %in% c("renyi", "tsallis")) {
  if (opt$lambda <= 0 || opt$lambda == 1) 
    stop("Parameter λ must satisfy: λ > 0 and λ ≠ 1")
}

# ----------------------------------------------------------------------
# Load image matrix
# ----------------------------------------------------------------------
cat("Loading data (type =", opt$input_type, ")...\n")
if (opt$input_type == "envi") {
  stopifnot(!is.null(opt$image), !is.null(opt$header))
  img_mat <- myread.ENVI(file = opt$image, headerfile = opt$header)
} else {
  stopifnot(!is.null(opt$sim_rdata))
  load(opt$sim_rdata)   # loads object Z
  img_mat <- Z
}
rows <- nrow(img_mat); cols <- ncol(img_mat)
cat("Dimensions:", rows, "x", cols, "pixels\n")

# ----------------------------------------------------------------------
# Select estimator dynamically
# ----------------------------------------------------------------------
get_estimator <- function(type, L, bootstrap) {
  if (type == "shannon") {
    if (bootstrap && L > 1) return(function(z) bootstrap_al_omari_1_estimator(z, opt$B))
    else                     return(al_omari_1_estimator)
  }
  if (type == "renyi") {
    if (bootstrap && L > 1) return(function(z) bootstrap_renyi_entropy_estimator_v1(z, opt$B, opt$lambda))
    else                     return(function(z) renyi_entropy_estimator_v1(z, opt$lambda))
  }
  if (type == "tsallis") {
    if (bootstrap && L > 1) return(function(z) bootstrap_tsallis_entropy_optimized(z, opt$B, opt$lambda))
    else                     return(function(z) tsallis_estimator_optimized(z, opt$lambda))
  }
  stop("Unknown entropy type")
}
entropy_estimator <- get_estimator(opt$entropy, opt$looks, opt$bootstrap)

get_theoretical <- function(type, L, lambda, mu) {
  if (type == "shannon")  return(shannon_theoretical(L, mu))
  if (type == "renyi")    return(renyi_theoretical(L, lambda, mu))
  if (type == "tsallis")  return(tsallis_theoretical(L, lambda, mu))
}

# ----------------------------------------------------------------------
# Sliding-window processing (parallelized with future.apply)
# ----------------------------------------------------------------------
cat("Processing", opt$window, "x", opt$window, "windows...\n")
start_time <- Sys.time()

out_nrow     <- rows - opt$window + 1
out_ncol     <- cols - opt$window + 1
entropy_test <- matrix(NA_real_, nrow = out_nrow, ncol = out_ncol)

# Generate all top-left coordinates for each window
idx <- expand.grid(i = seq_len(out_nrow), j = seq_len(out_ncol))

# Function to process a single window
compute_one <- function(k) {
  i <- idx$i[k]
  j <- idx$j[k]
  slice <- img_mat[i:(i + opt$window - 1), j:(j + opt$window - 1)]
  est   <- entropy_estimator(slice)
  mu    <- mean(slice)
  theo  <- get_theoretical(opt$entropy, opt$looks, opt$lambda, mu)
  est - theo
}

# Process all windows, in parallel if enabled
if (parallel) {
  # Uses all available CPU cores by default
  res <- future_sapply(seq_len(nrow(idx)), compute_one, future.seed = TRUE)
} else {
  res <- sapply(seq_len(nrow(idx)), compute_one)
}
entropy_test[] <- matrix(res, nrow = out_nrow, byrow = FALSE)

# ----------------------------------------------------------------------
# p-value matrix
# ----------------------------------------------------------------------
cat("Computing p-values...\n")
p_values <- calculate_p_values_matrix(entropy_test)

# ----------------------------------------------------------------------
# Quick exploratory plots (commented)
# ----------------------------------------------------------------------
source("./Code/imagematrix_visualizer.R")
# dev.new()
# plot(imagematrix(equalize(img_mat)), main = "SAR image")
# plot(imagematrix(p_values), main = "p-values (greyscale)")
# plot(imagematrix(p_values < opt$p_thr), main = "p-values < threshold")
# previewImagematrix(imagematrix_color(p_values), palette_option = "viridis-H",
#                    main = "p-values (colour)")

# ----------------------------------------------------------------------
# Composite plot function
# ----------------------------------------------------------------------

plot_composite <- function(img_mat, p_values,
                           thr = 0.05,
                           palette_option = "viridis-H") {
  
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
  
  # SAR image
  plot(imagematrix(equalize(img_mat)), main = "SAR image")
  
  # grayscale p-values
  plot(imagematrix(p_values), main = "p-values (grayscale)")
  
  # color p-values 
  
  previewImagematrixPanel(
    imagematrix_color(p_values),
    palette_option     = palette_option,
    significance_level = thr
    #main               = "p-values (color)"
  )
  par(mar = c(1, 1, 3, 1))
  #  threshold map
  plot(imagematrix(p_values < thr), main = paste0("p-values < ", thr))
}

# ----------------------------------------------------------------------
# Show composite in RStudio panel
# ----------------------------------------------------------------------
cat("Plotting 2×2 composite figure...\n")
plot_composite(img_mat, p_values, thr = opt$p_thr)

# ----------------------------------------------------------------------
# Show only color p-values with scale bar
# ----------------------------------------------------------------------
# cat("Plotting color p-values with scale bar...\n")
# previewImagematrixPanel(
#   imagematrix_color(p_values),
#   palette_option     = "viridis-H",
#   significance_level = opt$p_thr
# )

# ----------------------------------------------------------------------
# Save PNG output
# ----------------------------------------------------------------------
cat("Saving PNG output...\n")
if (!dir.exists("./PNG")) dir.create("./PNG")

imagematrixPNG(
  imagematrix(equalize(img_mat)),
  file.path("./PNG", paste0(opt$prefix, "_sar.png"))
)

imagematrixPNG(
  imagematrix(equalize(entropy_test)),
  file.path("./PNG", paste0(opt$prefix, "_statistic.png"))
)

imagematrixPNG(
  imagematrix(p_values),
  file.path("./PNG", paste0(opt$prefix, "_pvals_gray.png"))
)

imagematrixPNG(
  imagematrix(p_values < opt$p_thr),
  file.path("./PNG", paste0(opt$prefix, "_pvals_thr.png"))
)
imagematrix_colorPNG(
  imagematrix_color(p_values),
  name            = file.path("./PNG", paste0(opt$prefix, "_pvals_color.png")),
  palette_option  = "viridis-H",
  legend_width_px = 300,
  scale_factor    = 2
)

# ----------------------------------------------------------------------
# Save workspace objects
# ----------------------------------------------------------------------
if (!dir.exists("./Data")) dir.create("./Data")
base_name <- sprintf("%s_%s_L%d%s_w%d%s%s",
                     opt$prefix,
                     ifelse(opt$entropy == "shannon", "shannon",
                            ifelse(opt$entropy == "renyi",   "renyi",   "tsallis")),
                     opt$looks,
                     ifelse(opt$bootstrap, paste0("_b", opt$B), ""),
                     opt$window,
                     if (opt$entropy %in% c("renyi","tsallis"))
                       paste0("_lambda", gsub("\\.", "", sprintf("%.2f", opt$lambda))) else "",
                     format(Sys.Date(), "_%Y%m%d"))
save(file = file.path("./Data", paste0(base_name, ".RData")),
     list = c("img_mat", "entropy_test", "p_values", "opt"))

# ----------------------------------------------------------------------
# Execution time formatted hh:mm:ss
# ----------------------------------------------------------------------
end_time <- Sys.time()
run_time <- difftime(end_time, start_time, units = "secs")
secs <- as.numeric(run_time)
h <- floor(secs / 3600)
m <- floor((secs %% 3600) / 60)
s <- round(secs %% 60)
cat(sprintf("Total run time: %02d:%02d:%02d (hh:mm:ss)\n", h, m, s))
