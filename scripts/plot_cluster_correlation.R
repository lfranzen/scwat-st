#' CLUSTER CORRELATION PLOTS
#' Hypergeometric test of jaccard index based on similarities in cluster marker genes.

library(ggplot2)
library(RColorBrewer)
library(magrittr)
library(readxl)
library(pheatmap)
library(GeneOverlap)

PROJECT_ID <- "visium"
DIR_ROOT <- getwd()
DIR_WD <- file.path(getwd(), "scripts")  # 
DIR_DATA <- file.path(DIR_ROOT, "data", PROJECT_ID)
DIR_RES <- file.path(DIR_ROOT, "results" , PROJECT_ID)
DIR_FIG <- file.path(DIR_RES, "figures")

source(file.path(DIR_WD, "colors.R"))


# ===================================
#' READ DATA
#' 'Baseline' Visium data
se_base <- readRDS(file.path(DIR_RES, "se-object.visium_baseline.rds"))
canno_base <- read.csv(file.path(DIR_RES, "tables", "visium_baseline.clustering_annotations.csv"), stringsAsFactors = F)
colnames(canno_base)[colnames(canno_base)=="seurat_clusters"] <- "cluster"

markers_list_base <- list()
for(c in sort(canno_base$cluster)){
  markers_list_base[[c]] <- as.data.frame(read_excel(file.path(DIR_RES, "tables", "visium_baseline.markers_clusterdea.xlsx"), sheet = paste0("cluster_", c)))
  markers_list_base[[c]] <- subset(markers_list_base[[c]], avg_logFC>0 & p_val_adj<0.05)
  markers_list_base[[c]]$cluster <- c
}
markers_base <- do.call(rbind.data.frame, markers_list_base)
markers_base <- merge(x=markers_base, y=canno_base, by="cluster") %>% dplyr::arrange(as.numeric(cluster))

#' 'Insulin' Visium data
se_ins <- readRDS(file.path(DIR_RES, "se-object.visium_insulin.rds"))
canno_ins <- read.csv(file.path(DIR_RES, "tables", "visium_insulin.clustering_annotations.csv"), stringsAsFactors = F)
colnames(canno_ins)[colnames(canno_ins)=="seurat_clusters"] <- "cluster"

markers_list_ins <- list()
for(c in sort(canno_ins$cluster)){
  markers_list_ins[[c]] <- as.data.frame(read_excel(file.path(DIR_RES, "tables", "visium_insulin.markers_clusterdea.xlsx"), sheet = paste0("cluster_", c)))
  markers_list_ins[[c]] <- subset(markers_list_ins[[c]], avg_logFC>0 & p_val_adj<0.05)
  markers_list_ins[[c]]$cluster <- c
}
markers_ins <- do.call(rbind.data.frame, markers_list_ins)
markers_ins <- merge(x=markers_ins, y=canno_ins, by="cluster") %>% dplyr::arrange(as.numeric(cluster))


# ===================================
#' JACCARD IDXs
computeJaccardMatrix <- function (cluster.marker.table, cluster.column.name = "cluster", padj.cutoff = 0.05, genome.size) {
  c_markers_filt <- cluster.marker.table[cluster.marker.table$p_val_adj < padj.cutoff, ]
  c_markers_filt <- c_markers_filt[order(c_markers_filt[,cluster.column.name]),]
  gene_sets <- split(c_markers_filt$gene, c_markers_filt[,cluster.column.name])
  
  jaccards <- c()
  p_values <- c()
  for (i in names(gene_sets)) {
    for (j in names(gene_sets)) {
      set1 <- gene_sets[[i]]
      set2 <- gene_sets[[j]]
      go.obj <- GeneOverlap::newGeneOverlap(listA = set1, listB = set2, genome.size = genome.size)
      go.obj <- GeneOverlap::testGeneOverlap(go.obj)
      jaccards <- c(jaccards, go.obj@Jaccard)
      p_values <- c(p_values, go.obj@pval)
    }
  }
  jaccards <- matrix(data = jaccards, nrow = length(gene_sets), ncol = length(gene_sets))
  colnames(jaccards) <- rownames(jaccards) <- names(gene_sets)
  
  return(jaccards)
}

jacc_base <- computeJaccardMatrix(markers_base, cluster.column.name = "cluster_anno", genome.size=nrow(se_base))
jacc_ins <- computeJaccardMatrix(markers_ins, cluster.column.name = "cluster_anno", genome.size=nrow(se_ins))


# ===================================
#' PLOT

#' Pheatmap plot function
plot_jaccard_pheatmap <- function (jaccards, 
                                   jaccard.max = 0.5,
                                   cluster.annotation.df, 
                                   column.name.use="cluster",
                                   hm.category = "cluster_group",
                                   hm.cluster.colors = "cluster_color",
                                   colors.hm = colorRampPalette(c("white", 'grey30'))(40)) {
  
  jaccards2 <- jaccards
  jaccards2[jaccards2 > jaccard.max] <- jaccard.max
  diag(jaccards2) <- NA
  
  groups_clus <- data.frame(cluster.annotation.df[, c(column.name.use, hm.category, hm.cluster.colors)], check.names = F)
  colnames(groups_clus) <- c("Cluster", "Category", "Color")
  rownames(groups_clus) <- groups_clus$Cluster
  groups_clus <- groups_clus[colnames(jaccards2), ]
  
  cats <- as.character(unique(groups_clus$Category))
  groups_colors <- c()
  for(c in cats){
    # print(c)
    color_add <- cluster.annotation.df[cluster.annotation.df[,hm.category] == c, hm.cluster.colors][1]
    groups_colors <- c(groups_colors, color_add)
  }
  names(groups_colors) <- cats
  
  cluster_colors <- groups_clus$Color
  names(cluster_colors) <- groups_clus$Cluster
  
  legend_colors <- list(Category = groups_colors,
                        Cluster = cluster_colors)
  
  p <- pheatmap::pheatmap(jaccards2,
                          annotation_row = groups_clus[,1:2],
                          annotation_col = groups_clus[,1:2],
                          annotation_colors = legend_colors,
                          color = colors.hm,
                          border_color = NA,
                          annotation_legend = F,
                          legend = T,
                          fontsize = 14)
  return(p)
}

#' plot
p1 <- plot_jaccard_pheatmap(jaccards = jacc_base, 
                            jaccard.max = 0.3,
                            cluster.annotation.df = canno_base, 
                            column.name.use = "cluster_anno", 
                            colors.hm = viridis::viridis(20))
 
p2 <- plot_jaccard_pheatmap(jaccards = jacc_ins, 
                            jaccard.max = 0.3,
                            cluster.annotation.df = canno_ins, 
                            column.name.use = "cluster_anno", 
                            colors.hm = viridis::viridis(20))



fname <- "plot_cluster_corr_jaccard.visium_baseline"
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 9.8, height = 9);p1;dev.off()
tiff(filename = file.path(DIR_FIG, paste0(fname, ".tiff")), units = "cm", width = 9.8*2.5, height = 9*2.5, res = 600); p1; dev.off()

fname <- "plot_cluster_corr_jaccard.visium_insulin"
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 9.8, height = 9);p2;dev.off()
tiff(filename = file.path(DIR_FIG, paste0(fname, ".tiff")), units = "cm", width = 9.8*2.5, height = 9*2.5, res = 600); p2; dev.off()


# ===================================
#' 'Baseline' vs 'Insulin' clustering

se_markers_filt_base <- markers_base[markers_base$p_val_adj < 0.01, ]
se_markers_filt_ins <- markers_ins[markers_ins$p_val_adj < 0.01, ]

GeneSets_ins <- split(se_markers_filt_ins$gene, se_markers_filt_ins$cluster_anno)
GeneSets_base <- split(se_markers_filt_base$gene, se_markers_filt_base$cluster_anno)

# Calculate hypergeometric test 
jaccards <- c()
p_values <- c()
for (i in names(GeneSets_ins)) {
  for (j in names(GeneSets_base)) {
    set1 <- GeneSets_ins[[i]]
    set2 <- GeneSets_base[[j]]
    go.obj <- GeneOverlap::newGeneOverlap(listA = set1, listB = set2, genome.size = nrow(se_base))
    go.obj <- GeneOverlap::testGeneOverlap(go.obj)
    jaccards <- c(jaccards, go.obj@Jaccard) #c(comparisons, go.obj@Jaccard)
    p_values <- c(p_values, go.obj@pval) # c(p_values, go.obj@pval)
  }
}

# Create matrix from jaccard distances
jaccards <- matrix(data = jaccards, nrow = length(GeneSets_base), ncol = length(GeneSets_ins))

rownames(jaccards) <- names(GeneSets_base)
colnames(jaccards) <- names(GeneSets_ins)

#' Prep plot
groups_clus_ins <- data.frame(row.names = canno_ins$cluster_anno,
                              Cluster_ins = canno_ins$cluster_anno,
                              Category = canno_ins$cluster_group,
                              stringsAsFactors = F
                              )
groups_clus_base <- data.frame(row.names = canno_base$cluster_anno,
                               Cluster_base = canno_base$cluster_anno,
                               Category = canno_base$cluster_group,
                               stringsAsFactors = F
                               )

legend_colors_ins <- canno_ins$cluster_color
names(legend_colors_ins) <- groups_clus_ins$Cluster_ins

legend_colors_base <- canno_base$cluster_color
names(legend_colors_base) <- groups_clus_base$Cluster_base

cats <- as.character(unique(groups_clus_ins$Category))
groups_colors <- c()
for(c in cats){
  # print(c)
  color_add <- canno_ins[canno_ins$cluster_group == c, "cluster_color"][1]
  groups_colors <- c(groups_colors, color_add)
}
names(groups_colors) <- cats

legend_colors <- list(Cluster_ins = legend_colors_ins,
                      Cluster_base = legend_colors_base,
                      Category = groups_colors)


#' plot
max_jacc <- 0.4
jaccards2 <- jaccards
jaccards2[jaccards2>max_jacc] <- max_jacc

p <- pheatmap::pheatmap(jaccards2, 
                        annotation_row = groups_clus_base,
                        annotation_col = groups_clus_ins,
                        annotation_colors = legend_colors,
                        color = viridis::viridis(20), 
                        border_color = NA, 
                        annotation_legend = F,
                        legend = T,
                        fontsize = 14)

fname <- "plot_cluster_corr_jaccard.comparison"
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 9.8, height = 9);p;dev.off()
tiff(filename = file.path(DIR_FIG, paste0(fname, ".tiff")), units = "cm", width = 9.8*2.5, height = 9*2.5, res = 600); p; dev.off()


#==========================================================================
#==========================================================================


