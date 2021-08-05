library(baydem)
rm(list = ls())

tikal_file <- file.path("outputs","tikal.rds")
if (!file.exists(tikal_file)) {
  stop("Missing tikal.rds")
}
tikal <- readRDS(tikal_file)

all_file <- file.path("outputs","all.rds")
if (!file.exists(all_file)) {
  stop("Missing all.rds")
}
all <- readRDS(all_file)

pdf(file.path("outputs","FigS3_tikal_waic.pdf"),width=8,height=6)
  plot(tikal$density_model$K,tikal$waic_vect,xlab="K",ylab="WAIC",main="Tikal")
dev.off()

pdf(file.path("outputs","FigS4_all_waic.pdf"),width=8,height=6)
  plot(all$density_model$K,all$waic_vect,xlab="K",ylab="WAIC",main="All")
dev.off()