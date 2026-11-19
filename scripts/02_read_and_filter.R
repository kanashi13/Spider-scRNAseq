source("scripts/00_setup.R")

read_rhapsody <- function(path, project) {
  dt <- fread(path)
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- paste0(project, "_", dt$Cell_Index)
  CreateSeuratObject(
    counts = as(t(mat), "dgCMatrix"),
    project = project,
    min.cells = 3,
    min.features = 200
  )
}

files <- c(
  st5l_SC   = "data/extracted/GSM6069859_SC-stage5l-01_RSEC_MolsPerCell.csv.gz",
  st5l_cell = "data/extracted/GSM6069860_st5l_cell01_RSEC_MolsPerCell.csv.gz",
  st5l_nuc  = "data/extracted/GSM6069861_st5l_nuc01_RSEC_MolsPerCell.csv.gz"
)

objs <- lapply(names(files), function(n) read_rhapsody(files[[n]], n))
names(objs) <- names(files)
saveRDS(objs, "data/processed/objs_raw_list.rds")

qc <- do.call(rbind, lapply(names(objs), function(n) {
  data.frame(
    sample = n,
    nFeature = objs[[n]]$nFeature_RNA,
    nCount = objs[[n]]$nCount_RNA
  )
}))
write.csv(qc, "results/tables/qc_per_cell.csv", row.names = FALSE)

p_feature <- ggplot(qc, aes(sample, nFeature, fill = sample)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.size = 0.3) +
  theme_bw() +
  theme(legend.position = "none")

p_count <- ggplot(qc, aes(sample, nCount, fill = sample)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.size = 0.3) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none")

p_scatter <- ggplot(qc, aes(nCount, nFeature, color = sample)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_x_log10() +
  facet_wrap(~sample) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("results/figures/qc_nFeature.png", p_feature, width = 6, height = 4, dpi = 300)
ggsave("results/figures/qc_nCount.png", p_count, width = 6, height = 4, dpi = 300)
ggsave("results/figures/qc_scatter.png", p_scatter, width = 9, height = 3.5, dpi = 300)

thresholds <- list(
  st5l_SC   = c(fmin = 2000, fmax = 5500, cmin = 0,    cmax = 45000),
  st5l_cell = c(fmin = 500,  fmax = 4500, cmin = 1000, cmax = 35000),
  st5l_nuc  = c(fmin = 1000, fmax = 4500, cmin = 1500, cmax = 18000)
)

objs_filtered <- lapply(names(objs), function(n) {
  t <- thresholds[[n]]
  subset(
    objs[[n]],
    nFeature_RNA >= t["fmin"] & nFeature_RNA <= t["fmax"] &
      nCount_RNA >= t["cmin"] & nCount_RNA <= t["cmax"]
  )
})
names(objs_filtered) <- names(objs)

filter_table <- data.frame(
  sample = names(objs),
  before = sapply(objs, ncol),
  after = sapply(objs_filtered, ncol)
)
print(filter_table)
write.csv(filter_table, "results/tables/filter_before_after.csv", row.names = FALSE)

saveRDS(objs_filtered, "data/processed/objs_filtered_list.rds")
