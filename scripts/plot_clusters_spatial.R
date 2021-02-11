# Plot clusters spatially


#' SET UP

#' Load libs
library(Seurat)
library(STutility)
library(ggplot2)


#' Decide data
ANALYSIS <- "baseline"


#' Define project paths
PROJECT_ID <- "visium"
DIR_ROOT <- file.path(getwd())  # 

DIR_DATA <- file.path(DIR_ROOT, "data", PROJECT_ID)
DIR_RES <- file.path(DIR_ROOT, "results" , PROJECT_ID)
DIR_FIG <- file.path(DIR_RES, "figures")


#' Colors
source(file.path(DIR_WD, "colors.R"))


#' READ DATA

metadata <- read.delim(file.path(DIR_DATA, "visium_sample_metadata.tsv"), sep = "\t", stringsAsFactors = F)

if (ANALYSIS == "baseline") {
  se <- readRDS(file.path(DIR_RES, "se-object.visium_baseline.rds"))
  c_anno <- read.csv(file = file.path(DIR_RES, "tables", "visium_baseline.clustering_annotations.csv"), stringsAsFactors = F)
  metadata <- subset(metadata, insulin_stim == 0 | is.na(insulin_stim))
} else if (ANALYSIS == "insulin") {
  se <- readRDS(file.path(DIR_RES, "se-object.visium_insulin.rds"))
  se <- AddMetaData(se, metadata = paste0(se$subject_id, "_", se$insulin_stim, "_", se$sample_id), col.name = "subject_id")
  c_anno <- read.csv(file = file.path(DIR_RES, "tables", "visium_insulin.clustering_annotations.csv"), stringsAsFactors = F)
  metadata <- subset(metadata, insulin_stim == 0 | insulin_stim == 1)
}


#' PLOT

#' All clusters on all samples
p <- ST.FeaturePlot(se,
                    features = "seurat_clusters",
                    cols = c_anno[order(c_anno$seurat_clusters), "cluster_color"],
                    pt.size = 1.2);p

fname <- paste0("plot_clusters_spatial_01.", ANALYSIS)
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 12.5, height = 10, useDingbats = F);p;dev.off()


#' Selected clusters
c_plot <- c(3, 4, 6, 9, 11, 14, 18, 19, 20)
s_plot <- c(2, 3)

fname <- paste0("plot_clusters_spatial_02.", ANALYSIS)
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 6, height = 6, useDingbats = F)
for (s in s_plot) {
  for (c in c_plot) {
    colors_cluster <- rep("grey95", length(unique(c_anno$seurat_clusters)))
    colors_cluster[c] <- c_anno[c_anno$seurat_clusters == c, "cluster_color"]
    p <- ST.FeaturePlot(se, 
                       indices = s, 
                       features = "seurat_clusters",
                       cols = colors_cluster,
                       pt.size = 2)
    print(p)
  }
}; dev.off()


#' Selected cluster-pairs
c_pair_plot <- list(c(3, 4), c(8, 9), c(17, 19), c(6, 12, 15), c(9,11), c(14,15))
s_plot <- c(3)

fname <- paste0("plot_clusters_spatial_03.", ANALYSIS)
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 6, height = 6, useDingbats = F)
for (s in s_plot) {
  for (c_pair in c_pair_plot) {
    colors_cluster <- rep("grey95", length(unique(c_anno$seurat_clusters)))
    colors_cluster[c_pair] <- c_anno[c_anno$seurat_clusters %in% c_pair, "cluster_color"]
    p <- ST.FeaturePlot(se, 
                        indices = s, 
                        features = "seurat_clusters",
                        cols = colors_cluster,
                        pt.size = 2)
    print(p)
  }
}; dev.off()


