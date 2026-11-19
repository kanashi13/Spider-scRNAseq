suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

dirs <- c(
  "data/raw",
  "data/extracted",
  "data/processed",
  "results/figures",
  "results/tables"
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
