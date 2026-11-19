source("scripts/00_setup.R")

objs <- readRDS("data/processed/objs_filtered_list.rds")

merged <- merge(objs[[1]], y = objs[-1])
merged <- JoinLayers(merged)
merged$sample <- sub("_[0-9]+$", "", colnames(merged))
print(table(merged$sample))

merged[["RNA"]] <- split(merged[["RNA"]], f = merged$sample)

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 30)
merged <- RunUMAP(
  merged,
  reduction = "pca",
  dims = 1:20,
  reduction.name = "umap_before"
)

p_before <- DimPlot(merged, reduction = "umap_before", group.by = "sample") +
  ggtitle("Before integration")
ggsave("results/figures/umap_before_integration.png", p_before, width = 6, height = 5, dpi = 300)

merged <- IntegrateLayers(
  merged,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "rpca",
  k.weight = 50,
  verbose = FALSE
)
merged <- JoinLayers(merged)

merged <- FindNeighbors(merged, reduction = "rpca", dims = 1:20)
merged <- FindClusters(merged, resolution = 0.5)
merged <- RunUMAP(
  merged,
  reduction = "rpca",
  dims = 1:20,
  reduction.name = "umap_rpca"
)

p_after <- DimPlot(merged, reduction = "umap_rpca", group.by = "sample") +
  ggtitle("After RPCA integration")
p_clusters <- DimPlot(
  merged,
  reduction = "umap_rpca",
  group.by = "seurat_clusters",
  label = TRUE
)
ggsave("results/figures/umap_rpca.png", p_after + p_clusters, width = 12, height = 5, dpi = 300)

saveRDS(merged, "data/processed/merged_rpca.rds")
