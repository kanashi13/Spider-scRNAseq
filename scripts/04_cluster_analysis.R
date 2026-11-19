source("scripts/00_setup.R")

merged <- readRDS("data/processed/merged_rpca.rds")
Idents(merged) <- "seurat_clusters"

p_umap <- DimPlot(merged, reduction = "umap_rpca", label = TRUE, repel = TRUE) +
  ggtitle(paste0("Spider embryo clusters (n = ", ncol(merged), " cells)"))
ggsave("results/figures/umap_clusters.png", p_umap, width = 7, height = 6, dpi = 300)

p_sample <- DimPlot(merged, reduction = "umap_rpca", group.by = "sample") +
  ggtitle("Samples")
ggsave("results/figures/umap_samples.png", p_sample, width = 7, height = 6, dpi = 300)

p_qc <- VlnPlot(
  merged,
  features = c("nCount_RNA", "nFeature_RNA"),
  group.by = "seurat_clusters",
  pt.size = 0
)
ggsave("results/figures/cluster_qc_violin.png", p_qc, width = 10, height = 4, dpi = 300)

cluster_counts <- as.data.frame(table(cluster = merged$seurat_clusters))
write.csv(cluster_counts, "results/tables/cluster_cell_counts.csv", row.names = FALSE)

sample_table <- table(cluster = merged$seurat_clusters, sample = merged$sample)
write.csv(
  as.data.frame.matrix(sample_table),
  "results/tables/cluster_sample_counts.csv"
)
write.csv(
  round(as.data.frame.matrix(prop.table(sample_table, 1)), 3),
  "results/tables/cluster_sample_proportions.csv"
)

qc_stats <- merged@meta.data |>
  group_by(cluster = seurat_clusters) |>
  summarise(
    cells = n(),
    median_nCount = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA)
  )
write.csv(qc_stats, "results/tables/cluster_qc_statistics.csv", row.names = FALSE)

composition <- as.data.frame(sample_table)
p_comp <- ggplot(composition, aes(cluster, Freq, fill = sample)) +
  geom_col(position = "fill") +
  labs(x = "Cluster", y = "Proportion of cells") +
  theme_bw()
ggsave("results/figures/cluster_sample_composition.png", p_comp, width = 8, height = 4, dpi = 300)

markers <- FindAllMarkers(
  merged,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
)
write.csv(markers, "results/tables/cluster_markers_all.csv", row.names = FALSE)

significant <- markers |> filter(p_val_adj < 0.05)
write.csv(significant, "results/tables/cluster_markers_significant.csv", row.names = FALSE)

top10 <- significant |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) |>
  ungroup()
write.csv(top10, "results/tables/cluster_markers_top10.csv", row.names = FALSE)

top20 <- significant |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 20, with_ties = FALSE) |>
  ungroup()
write.csv(top20, "results/tables/cluster_markers_top20.csv", row.names = FALSE)

saveRDS(merged, "data/processed/merged_rpca_clustered.rds")
