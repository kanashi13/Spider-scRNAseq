source("scripts/00_setup.R")

tar_file <- "data/raw/GSE201705_RAW.tar"
geo_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201705/suppl/GSE201705_RAW.tar"

if (!file.exists(tar_file)) {
  options(timeout = 3600)
  download.file(geo_url, destfile = tar_file, mode = "wb")
}

untar(tar_file, exdir = "data/extracted")

extracted <- list.files("data/extracted", full.names = TRUE)
manifest <- data.frame(
  file_name = basename(extracted),
  size_MB = round(file.info(extracted)$size / 1024^2, 3)
)
print(manifest)
write.csv(manifest, "results/tables/file_manifest.csv", row.names = FALSE)
