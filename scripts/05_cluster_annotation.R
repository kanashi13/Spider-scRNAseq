source("scripts/00_setup.R")
library(rentrez)

merged <- readRDS("data/processed/merged_rpca_clustered.rds")
markers <- read.csv("results/tables/cluster_markers_significant.csv")

top <- markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 50, with_ties = FALSE) |>
  ungroup()

fetch_summary <- function(ids, attempts = 3) {
  for (i in seq_len(attempts)) {
    result <- tryCatch(entrez_summary(db = "gene", id = ids), error = function(e) NULL)
    if (!is.null(result)) {
      if (inherits(result, "esummary")) result <- setNames(list(result), ids)
      return(result)
    }
    Sys.sleep(2)
  }
  NULL
}

gene_ids <- unique(sub("^LOC", "", top$gene))
gene_ids <- gene_ids[grepl("^[0-9]+$", gene_ids)]
chunks <- split(gene_ids, ceiling(seq_along(gene_ids) / 50))

annotation <- do.call(rbind, lapply(chunks, function(ch) {
  s <- fetch_summary(ch)
  Sys.sleep(0.4)
  if (is.null(s)) return(NULL)
  data.frame(
    gene = paste0("LOC", names(s)),
    symbol = sapply(s, function(x) x$name),
    description = sapply(s, function(x) x$description),
    stringsAsFactors = FALSE
  )
}))

top_annotated <- left_join(top, annotation, by = "gene")
write.csv(top_annotated, "results/tables/cluster_markers_top50_annotated.csv", row.names = FALSE)

cluster_labels <- c(
  "0"  = "C0 patterning (weak)",
  "1"  = "C1 ventral ectoderm EGFR?",
  "2"  = "C2 neuroectoderm/segmental?",
  "3"  = "C3 anterior head/brain?",
  "4"  = "C4 BMP/dorsal?",
  "5"  = "C5 yolk/gut lipid?",
  "6"  = "C6 mesoderm?",
  "7"  = "C7 yolk/gut lipid?",
  "8"  = "C8 epidermis?",
  "9"  = "C9 proliferating?",
  "10" = "C10 midline/glia?"
)
merged$cluster_label <- unname(cluster_labels[as.character(merged$seurat_clusters)])

p_labeled <- DimPlot(
  merged,
  reduction = "umap_rpca",
  group.by = "cluster_label",
  label = TRUE,
  repel = TRUE,
  label.size = 3
) + NoLegend()
ggsave("results/figures/umap_labeled.png", p_labeled, width = 8, height = 6, dpi = 300)

described <- top_annotated |>
  filter(!is.na(description), !grepl("^uncharacterized", description))

top_features <- described |>
  group_by(cluster) |>
  slice_head(n = 3) |>
  ungroup() |>
  distinct(gene, .keep_all = TRUE)
top_features <- top_features |> filter(gene %in% rownames(merged))

p_top <- DotPlot(merged, features = top_features$gene, group.by = "seurat_clusters") +
  RotatedAxis() +
  scale_x_discrete(labels = setNames(substr(top_features$description, 1, 32), top_features$gene)) +
  labs(x = "Top marker genes", y = "Cluster") +
  theme(axis.text.x = element_text(size = 7))
ggsave("results/figures/dotplot_top_markers.png", p_top, width = 14, height = 5, dpi = 300)

reference_pattern <- "OTX2|Sox-3|twist|Hox|barH|DLX|T-box|GATA|endothelial transcription|cyclin|slit|sonic|noggin|chordin|gooseberry"
reference_features <- described |>
  filter(grepl(reference_pattern, description, ignore.case = TRUE)) |>
  distinct(gene, .keep_all = TRUE) |>
  filter(gene %in% rownames(merged))

if (nrow(reference_features) > 0) {
  p_ref <- DotPlot(merged, features = reference_features$gene, group.by = "seurat_clusters") +
    RotatedAxis() +
    scale_x_discrete(labels = setNames(substr(reference_features$description, 1, 32), reference_features$gene)) +
    labs(x = "Reference genes", y = "Cluster") +
    theme(axis.text.x = element_text(size = 7))
  ggsave("results/figures/dotplot_reference_genes.png", p_ref, width = 14, height = 5, dpi = 300)
}

saveRDS(merged, "data/processed/merged_rpca_labeled.rds")
