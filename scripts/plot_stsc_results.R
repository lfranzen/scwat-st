########################
#' Stereoscope results
#' 
#' Plot output from Stereoscope (A Andersson, 2020). scRNA-seq data from Grundberg publication (Vijay et al, 2020).
#' 
#' L. Franzén, lovisa.franzen@scilifelab.se
#' Feb 2020
########################

# ===================================
#' SET UP

#' Select analysis:
ANALYSIS <- "baseline"
# ANALYSIS <- "insulin"

# STSC_RUN <- "stereoscope_200925"
STSC_RUN <- "stereoscope_210322"
# STSC_RUN <- "stereoscope_210407"

ANALYSIS_ID <- paste0("stsc", strsplit(STSC_RUN, "_")[[1]][2], "_", ANALYSIS)

#' Load libs
library(Seurat)
library(STutility)

library(magrittr)

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scico)


#' Define paths
DIR_ROOT <- getwd()
DIR_WD <- file.path(getwd(), "scripts")
DIR_DATA <- file.path(DIR_ROOT, "data", STSC_RUN)
DIR_RES <- file.path(DIR_ROOT, "results" , "visium")
DIR_FIG <- file.path(DIR_RES, "figures")

source(file.path(DIR_WD, "colors.R"))
  
# ===================================
#' DEF FUNCTIONS
read_st_datasets <- function (data_path, file_pattern = "stdata", transpose_data = TRUE) {
  #' @description Function to read multiple st count matrices and store as a list of matrices
  #' @param data_path Path to where the raw data is located (data format given as spots in rows and genes in columns)
  #' @param file_pattern Regex pattern by which the count .tsv files are recognised
  #' @param transpose_data Whether to transpose the input data matrix in order to obtain gene id in rows and spot id in cols (default is TRUE)
  #' @return List with count matrices with genes in rows and spots in columns
  
  if (!require(data.table)) {install.packages('data.table')}
  library(data.table)
  
  file_names <- list.files(path = data_path, pattern = file_pattern, full.names = F, recursive = T)
  
  exp_data_list <- lapply(1:length(file_names), function(i) {
    # Read input data from list of filenames
    file_path <- file.path(data_path, file_names[i])
    message("Reading file: ", file_path)
    
    if (transpose_data == TRUE) {
      exp_matrix <- t(data.frame(fread(input = file_path, sep = "\t"), row.names = 1))
    } else {
      exp_matrix <- (data.frame(fread(input = file_path, sep = "\t"), row.names = 1))
    }
    
    return(exp_matrix)
  })
  exp_data_list <- setNames(exp_data_list, nm = as.character( seq(1:length(file_names)) ))
  
  return(exp_data_list)
}

merge_st_datasets <- function (exp_data_list, labels = 1:length(exp_data_list), label_spot_sep = "x", label_first = TRUE) {
  #' @description Function used to merge expression datasets and store in sparse matrix (dgCMatrix) format
  #' @param exp_data_list List of expression matrices
  #' @param labels Labels for each dataset, must be same length as exp_data_list (default is 1:n)
  #' @param label_spot_sep Define character to separate sample label to spot identifier
  #' @param label_first Add sample label before spot identifier. If FALSE, add after (default TRUE)
  #' @return Sparse matrix with all given datasets merged
  
  if (!require(Matrix)) {install.packages('Matrix')}
  library(Matrix)
  stopifnot(length(exp_data_list) == length(labels))
  
  genes_union <- unique(unlist(lapply(exp_data_list, rownames)))
  
  #' Start with first dataset
  message("Attaching dataset 1")
  A <- as.data.frame(exp_data_list[[1]])
  if (label_first == TRUE) {
    colnames(A) <- paste(labels[1], colnames(A), sep = label_spot_sep)
  } else {
    colnames(A) <- paste(colnames(A), labels[1], sep = label_spot_sep)
  }
  A <- A[genes_union, ]
  A[is.na(A)] <- 0
  A <- as(object = as.matrix(A), Class = "dgCMatrix")
  
  #' Attach following datasets to the first
  for (i in 2:length(exp_data_list)) {
    message("Attaching dataset ", i)
    B <- as.data.frame(exp_data_list[[i]])
    if (label_first == TRUE) {
      colnames(B) <- paste(labels[i], colnames(B), sep = label_spot_sep)
    } else {
      colnames(B) <- paste(colnames(B), labels[i], sep = label_spot_sep)
    }
    
    B <- B[genes_union, ]
    B[is.na(B)] <- 0
    B <- as(object = as.matrix(B), Class = "dgCMatrix")
    
    A <- cbind(A, B)
  }
  return(A)
}


# ===================================
#' READ DATA

#' Stsc
path_data_files <- paste0(DIR_DATA, "/", list.files(DIR_DATA, recursive = T))
sample_fnames <- gsub(pattern = ".tsv", replacement = "", grep(pattern = "ADI_", x = unlist(strsplit(x = path_data_files, split = "/")), value = T))
sample_names <- grep(pattern = "ADI_", x = unlist(strsplit(x = sample_fnames, split = "-")), value = T)

stsc_data_raw <- read_st_datasets(data_path = DIR_DATA, file_pattern = "W.20")
names(stsc_data_raw) <- sample_names


#' se obj
if (ANALYSIS=="baseline") {
  #' 'Baseline' Visium data
  message("Read seurat object from analysis ", ANALYSIS)
  c_anno <- read.csv(file.path(DIR_RES, "tables", "visium_baseline.clustering_annotations.csv"), stringsAsFactors = F)
  se <- readRDS(file.path(DIR_RES, "se-object.visium_baseline.rds"))
  # se <- LoadImages(se, time.resolve = T, verbose = T, xdim = 100)
} else if (ANALYSIS=="insulin") {
  #' 'Insulin' Visium data
  message("Read seurat object from analysis ", ANALYSIS)
  c_anno <- read.csv(file.path(DIR_RES, "tables", "visium_insulin.clustering_annotations.csv"), stringsAsFactors = F)
  se <- readRDS(file.path(DIR_RES, "se-object.visium_insulin.rds"))
  # se <- LoadImages(se, time.resolve = T, verbose = T, xdim = 100)
  se <- AddMetaData(se, metadata = paste0(se$subject_id, "_", se$insulin_stim), col.name = "subject_id")
}


#' Subset and merge stsc data
sample_ids <- unique(se$sample_name)
sample_rep <- unique(se$seu_n)
stsc_data_raw_subset <- stsc_data_raw[sample_ids]

stsc_data_merged <- merge_st_datasets(exp_data_list = stsc_data_raw_subset, 
                                      labels = sample_rep,
                                      label_spot_sep = "_", 
                                      label_first = FALSE)

stsc_data_merged_t <- t(stsc_data_merged)
rownames(stsc_data_merged_t) <- gsub(pattern = "\\.", replacement = "-", rownames(stsc_data_merged_t))
colnames(stsc_data_merged_t) <- paste0("cell.", colnames(stsc_data_merged_t))


#' Subset se data and add stsc data to metadata
spots_intersect <- intersect(x = colnames(se), y =  rownames(stsc_data_merged_t)); print(length(spots_intersect))
se_subset <- SubsetSTData(se, spots = spots_intersect)
se_subset <- AddMetaData(object = se_subset, metadata = stsc_data_merged_t[spots_intersect, ])


#' ==================================================================
#' PLOT

#' Prep plot data
fname <- paste0(ANALYSIS_ID, ".overlap_seu_clusters")

stsc_cell_feats <- colnames(stsc_data_merged_t)
stsc_cell_feats2 <- gsub("cell.", "", stsc_cell_feats)

d_plot <- se_subset@meta.data[, c("seurat_clusters", "cluster_anno", "cluster_group", "sample_name", stsc_cell_feats)]
colnames(d_plot) <- gsub("cell.", "", colnames(d_plot))
d_plot_l <- reshape2::melt(d_plot, measure.vars=stsc_cell_feats2, variable.name="cell")

d_plot_l$order_group <- as.numeric(as.factor(d_plot_l$cluster_group))


#' Summarised stats by cluster
d_plot_stats <- d_plot_l %>%
  dplyr::group_by(seurat_clusters, cell) %>%
  dplyr::summarise(value_median = median(value, na.rm = TRUE),
                   value_mean = mean(value, na.rm = TRUE))


#' Export plot data/tables
write.table(d_plot, file.path(DIR_RES, "tables", paste0(ANALYSIS_ID, ".plot_data_wide.tsv")), quote = F, sep = "\t", row.names = T, col.names = T)
write.table(d_plot_l, file.path(DIR_RES, "tables", paste0(ANALYSIS_ID, ".plot_data_long.tsv")), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(d_plot_stats, file.path(DIR_RES, "tables", paste0(ANALYSIS_ID, ".plot_data_medians.tsv")), quote = F, sep = "\t", row.names = F, col.names = T)


#' VlnPlot + Box (Seurat)
feat_plot <- grep("immune", stsc_cell_feats, value = T)
pdf(file = file.path(DIR_FIG, paste0(fname, "_violinbox.pdf")), width = 16, height = 28, useDingbats = F)
VlnPlot(object = se_subset,
        features = sort(stsc_cell_feats),
        group.by = "cluster_anno",
        pt.size = 0,
        cols = c_anno$cluster_color,
        ncol = 3
        ) &
  geom_boxplot(width=0.6, outlier.size = .5) &
  theme(legend.position = 'none') &
  labs(x="")
dev.off()


#' Boxplot 1 (ggplot2)
plots <- list()
for(c in sort(stsc_cell_feats2)) {
  plots[[c]] <- ggplot(subset(d_plot_l, cell==c), aes(x=reorder(cluster_anno, order_group), y=value, fill=cluster_anno)) +
    geom_boxplot(width=.8, outlier.size = .2, outlier.colour = "grey80", color = "black", lwd=0.3) +
    labs(x="", y=c) +
    facet_grid(~cell) + 
    scale_fill_manual(values = c_anno$cluster_color) +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
          # strip.text = element_text(size = 10), 
          strip.background = element_blank(),
          panel.grid = element_blank()) +
    NoLegend()
}
pdf(file = file.path(DIR_FIG, paste0(fname, "_box.pdf")), width = 14, height = 26, useDingbats = F)
patchwork::wrap_plots(plots, ncol=3)
dev.off()


#' Boxplot 2 (ggplot2)
plots <- list()
for(c in sort(stsc_cell_feats2)) {
  plots[[c]] <- ggplot(subset(d_plot_l, cell==c), aes(x=reorder(seurat_clusters, order_group), y=value, fill=cluster_anno)) +
    geom_boxplot(width=.8, outlier.size = .2, outlier.colour = "grey80", color = "black", lwd=0.3) +
    labs(x="", y=c) +
    facet_grid(~cell) + 
    scale_fill_manual(values = c_anno$cluster_color) +
    theme_linedraw() +
    theme(
      # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
      # strip.text = element_text(size = 10), 
      strip.background = element_blank(),
      panel.grid = element_blank()) +
    NoLegend()
}
pdf(file = file.path(DIR_FIG, paste0(fname, "_box2.pdf")), width = 14, height = 20, useDingbats = F)
patchwork::wrap_plots(plots, ncol=3)
dev.off()


#' no outliers
plots <- list()
for(c in sort(stsc_cell_feats2)) {
  plots[[c]] <- ggplot(subset(d_plot_l, cell==c), aes(x=reorder(seurat_clusters, order_group), y=value, fill=cluster_anno)) +
    geom_boxplot(width=.8, outlier.shape = NA, color = "black", lwd=0.3) +
    labs(x="", y=c) +
    facet_grid(~cell) + 
    scale_fill_manual(values = c_anno$cluster_color) +
    coord_cartesian(ylim = c(0, .25)) +
    theme_linedraw() +
    theme(
      strip.background = element_blank(),
      panel.grid = element_blank()) +
    NoLegend()
}
pdf(file = file.path(DIR_FIG, paste0(fname, "_box2-nooutliers.pdf")), width = 14, height = 20, useDingbats = F)
patchwork::wrap_plots(plots, ncol=3)
dev.off()


#' Boxplot 3 (ggplot2): selected samples
feat_plot <- c("E1", "E2", "IS4", "IS5", "IS9", "IS10", "P4", "P6")
txt_size <- 8
plots <- list()
for(c in feat_plot) {
  plots[[c]] <- ggplot(subset(d_plot_l, cell==c), aes(x=reorder(seurat_clusters, order_group), y=value, fill=cluster_anno)) +
    geom_boxplot(width=.8, outlier.size = .2, outlier.colour = "grey80", color = "black", lwd=0.3) +
    labs(x="", y="") +
    facet_grid(~cell) + 
    coord_flip() +
    scale_fill_manual(values = c_anno$cluster_color) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = txt_size),
      axis.text.y = element_text(size = txt_size),
      strip.background = element_rect(fill = "white", colour = "white"),
      strip.text = element_text(colour = "black", size = txt_size),
      # strip.background = element_blank(),
      panel.grid = element_blank()) +
    NoLegend()
}
pdf(file = file.path(DIR_FIG, paste0(fname, "_box3b.pdf")), width = 4, height = 12, useDingbats = F)
patchwork::wrap_plots(plots, ncol=2)
dev.off()


#' Bar plot of medians
txt_size <- 6
p <- ggplot(d_plot_stats, aes(x=seurat_clusters, y=value_median, fill=seurat_clusters)) +
  geom_col() +
  facet_wrap(~cell, scales = "free_y", ncol = 3) +
  labs(x="", y="Spot proportion median value") +
  scale_fill_manual(values = c_anno[order(c_anno$seurat_clusters), "cluster_color"]) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = txt_size),
    axis.text.y = element_text(size = txt_size),
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(colour = "black", size = txt_size),
    panel.grid = element_blank()) + NoLegend(); p
pdf(file = file.path(DIR_FIG, paste0(fname, "_bars-median.pdf")), width = 8, height = 12, useDingbats = F);p;dev.off()

#' selected clusters
c_rm <- c(1,2,3,4,6)
col_fill_d <- c_anno[!c_anno$seurat_clusters %in% c_rm, ] %>% dplyr::arrange(seurat_clusters)
p2 <- ggplot(d_plot_stats[!d_plot_stats$seurat_clusters %in% c_rm, ], aes(x=seurat_clusters, y=value_median, fill=seurat_clusters)) +
  geom_col() +
  geom_hline(yintercept = 0, color="black", lwd=.3) +
  facet_wrap(~cell, scales = "free", ncol = 5) +
  labs(x="", y="Spot proportion median value") +
  scale_fill_manual(values = col_fill_d$cluster_color) +
  coord_flip() +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = txt_size),
    axis.text.y = element_text(size = txt_size),
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(colour = "black", size = txt_size),
    panel.grid = element_blank()) + NoLegend(); p2
pdf(file = file.path(DIR_FIG, paste0(fname, "_bars-median_subset.pdf")), width = 8, height = 12, useDingbats = F);p2;dev.off()


#' Correlation plot
M <- cor(d_plot[, stsc_cell_feats2])
pdf(file = file.path(DIR_FIG, paste0(fname, "_corrmatrix.pdf")), width = 8, height = 7, useDingbats = F)
corrplot(M, method="pie", 
         type="upper",
         # order="hclust", 
         col=brewer.pal(n=8, name="RdYlBu"), 
         tl.col="black")
dev.off()


#' UMAPs
feat_plot <- paste0("cell.", c("E1", "E2", "IS4", "IS5", "IS9", "IS10", "P4", "P6"))

p <- FeaturePlot(se_subset, reduction = "umap", features = feat_plot, 
                ncol=2,
                # min.cutoff = 0, max.cutoff = .2,
                cols = c("grey95", color_high2) 
                ) & theme(axis.title = element_blank(), 
                          axis.line = element_blank(), 
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

pdf(file = file.path(DIR_FIG, paste0(fname, "_umap.pdf")), width = 10, height = 18, useDingbats = F);p;dev.off()


feat_plot2 <- paste0("cell.", c("E2", "IS9", "P4"))
p2 <- FeaturePlot(se_subset, reduction = "umap", features = feat_plot2, 
                 ncol=1,
                 min.cutoff = 0, max.cutoff = .4,
                 cols = c("grey95", color_high2) 
                 ) & theme(axis.title = element_blank(), 
                            axis.line = element_blank(), 
                            axis.text = element_blank(),
                            axis.ticks = element_blank());p2
pdf(file = file.path(DIR_FIG, paste0(fname, "_umap2.pdf")), width = 5, height = 12, useDingbats = F);p2;dev.off()


#' Spatial plots
# feat_plot <- paste0("cell.", c("E1", "E2", "IS4", "IS4", "IS9", "P6"))
# ST.FeaturePlot(object = se_subset,
#                features = feat_plot, 
#                indices = 4)


#' ==================================================================



