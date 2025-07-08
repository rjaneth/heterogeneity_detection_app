# SAR Heterogeneity Detection App

This repository contains:

- A **Shiny application** for interactive detection of heterogeneity in SAR images.
- A **standalone R script** (`heterogeneity_detection_v1.R`) for processing large images using parallel computing (non-interactive).

---

## 🛰️ App Features

- Load simulated images, ENVI examples, or your own `.img` + `.hdr` files.
- Choose entropy estimator: **Shannon**, **Rényi**, or **Tsallis**.
- Optional **bootstrap** resampling and **sliding-window** analysis.
- Visualize interactive p-value map with color scale.

---

## 🧭 How the Interface Works

1. **Data source**: Select one of the three options:
   - `Simulated Images Examples`: two predefined simulated images.
   - `ENVI Images Examples`: three sample image+header pairs for testing.
   - `Upload ENVI`: upload your own `.img` and `.hdr` files.

2. **Choose image**: Based on the data source, select an image from the dropdown.

3. **Entropy type**: Select the estimator:
   - `Shannon`, `Rényi`, or `Tsallis`.

4. **Enable bootstrap (optional)**:
   - If enabled, set:
     - Number of bootstrap replicates (`B`)
     - Sliding window size (in pixels)

5. Click **Run Detection** to start the analysis.

> ⚠️ On the online version, use small values of `B` (e.g., 5–10) to avoid slow performance.

---

## 🚀 Running the App Locally

1. Clone the repository:
   ```bash
   git clone https://github.com/rjaneth/heterogeneity_detection_app.git
   ```
2. Open the project in **RStudio**.
3. Open the file `app.R`.
4. Click **Run App**.

---

## ⚙️ Using the Script (`heterogeneity_detection_v1.R`)

This script is designed for **large SAR images** and allows **faster processing** via **parallel computing**. No Shiny interface.

### Steps:

1. Open the file `heterogeneity_detection_v1.R` in RStudio.
2. Edit the `opt` list to set your:
   - Input `.img` and `.hdr` paths
   - Entropy type, window size, bootstrap options, number of replicates (`B`), etc.
3. Run the script:
   ```r
   source("heterogeneity_detection_v1.R")
   ```

### Output:

- A 2×2 composite plot with:
  - Equalized SAR image
  - p-value maps (color and grayscale)
  - Significance map (p < threshold)
- Saved outputs:
  - High-resolution PNGs in `./PNG/`
  - `.RData` file with results in `./Data/`

---

## 📁 Folder Structure

```
heterogeneity_detection_app/
├── app.R                          # Shiny app
├── heterogeneity_detection_v1.R   # Standalone script
├── Code/                          # Entropy and bootstrap functions
├── Data/                          # ENVI image samples
├── www/                           # Example images only for Shiny app
```

---

## 🌐 Online Demo

Try the hosted Shiny app (free tier — limited usage):

👉 [https://janeth-alpala.shinyapps.io/heterogeneity_detection_app](https://janeth-alpala.shinyapps.io/heterogeneity_detection_app)

> ⚠️ The free plan allows **15 hours/month** of usage. For best performance and full control, use the local version.
