########################
#' Betweeen-cluster analysis: "Heterotypic score algorithm"
#' 
#' Identify and quantify cluster-cluster neighbours. 
#' Visualize adjacency matrices as graphs and heatmaps. 
#' 
#' L. Franzén, lovisa.franzen@scilifelab.se
#' Jan 2021
########################

# ===================================
#' SET UP

#' Select analysis:
ANALYSIS <- "baseline"
# ANALYSIS <- "insulin"

#' Load libs
library(ggplot2)
library(magrittr)
library(STutility)

library(igraph)
library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggraph)

library(pheatmap)
library(scico)


#' Define paths
PROJECT_ID <- "visium"
DIR_ROOT <- getwd()
DIR_WD <- file.path(getwd(), "scripts")  # 
DIR_DATA <- file.path(DIR_ROOT, "data", PROJECT_ID)
DIR_RES <- file.path(DIR_ROOT, "results" , PROJECT_ID)
DIR_FIG <- file.path(DIR_RES, "figures")

source(file.path(DIR_WD, "colors.R"))


# ===================================
#' READ DATA & LOAD IMGS

if (ANALYSIS=="baseline") {
  #' 'Baseline' Visium data
  canno_base <- read.csv(file.path(DIR_RES, "tables", "visium_baseline.clustering_annotations.csv"), stringsAsFactors = F)
  se_base <- readRDS(file.path(DIR_RES, "se-object.visium_baseline.rds"))
  se_base <- LoadImages(se_base, time.resolve = T, verbose = T, xdim = 100)
  
  se <- se_base
  c_anno <- canno_base
  c_include <- seq(3, max(as.numeric(se[[]]$seurat_clusters)), 1)
  
  rm(se_base)
  
} else if (ANALYSIS=="insulin") {
  #' 'Insulin' Visium data
  canno_ins <- read.csv(file.path(DIR_RES, "tables", "visium_insulin.clustering_annotations.csv"), stringsAsFactors = F)
  se_ins <- readRDS(file.path(DIR_RES, "se-object.visium_insulin.rds"))
  se_ins <- LoadImages(se_ins, time.resolve = T, verbose = T, xdim = 100)
  
  se <- se_ins
  se <- AddMetaData(se, metadata = paste0(se$subject_id, "_", se$insulin_stim), col.name = "subject_id")
  c_anno <- canno_ins
  c_include <- seq(4, max(as.numeric(se[[]]$seurat_clusters)), 1)
  
  rm(se_ins)
}


# ===================================
#' Define functions

#' Create Adjacency Matrix
#' 
#' @param nbs.df Data frame with output from STUtility::GetSpatNet() data for all clusters of interest. Rows should correspond to spot ID and columns should include cluster IDs as well as colums starting with "nbs_".
#' @param column.clusters.id Column name in nbs.df corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param cluster.include Vector of cluster IDs to include in your analysis.
#' @return Adjacency matrix with number of neighbours present between each cluster pair
#' @export
CreateAdjMatrix <- function(
  nbs.df,
  column.clusters.id = "seurat_clusters",
  cluster.include
  ){
  nbs_adjmat <- matrix(0L, nrow = length(cluster.include), ncol = length(cluster.include))
  for (i in seq_along(cluster.include)) {
    c <- cluster.include[i]
    c_nbs_df <- nbs.df[nbs.df[,column.clusters.id]==c, ]
    for (j in seq_along(cluster.include)) {
      nbs <- paste0("nbs_", cluster.include[j])
      if (j==i) {
        n_nbs <- sum(!is.na(c_nbs_df[, nbs]==c))
      } else {
        n_nbs <- sum(!is.na(c_nbs_df[, nbs]==nbs))
      }
      nbs_adjmat[i,j] <- n_nbs
      if (nbs_adjmat[j,i] > 0) {
        nbs_adjmat[i,j] <- nbs_adjmat[j,i] <- max(c(nbs_adjmat[i,j], nbs_adjmat[j,i]))
      }
    }
    nbs_adjmat[i,i] <- sum(nbs_adjmat[i,][-i])
  }
  return(nbs_adjmat)
}

#' Create Adjacency Matrix per Subject
#' 
#' @param nbs.df Data frame with output from STUtility::GetSpatNet() data for all clusters of interest. Rows should correspond to spot ID and columns should include cluster IDs as well as colums starting with "nbs_".
#' @param se.metadata Metadata dataframe from seurat object containing Sample IDs for each spot.
#' @param column.clusters.id Column name in nbs.df corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param cluster.include Vector of cluster IDs to include in your analysis.
#' @param column.subject.id Column name in se.metadata corresponding to Sample ID of the spots.
#' @return List of adjacency matrices with number of neighbours present between each cluster pair. Each matrix in the list corresponds to data from one sample.
#' @export
CreateAdjMatrixPerSubject <- function(
  nbs.df,
  se.metadata,
  column.clusters.id = "seurat_clusters",
  cluster.include,
  column.subject.id
){
  nbs_adjmat_list <- list()
  subjects_include <- unique(as.character(se.metadata[, column.subject.id]))
  
  for(subject_id in subjects_include){
    rows_subject <- rownames(se.metadata[se.metadata[,column.subject.id] %in% subject_id, ])
    nbs.df_subject <- nbs.df[rows_subject, ]
    
    nbs_adjmat <- CreateAdjMatrix(nbs.df = nbs.df_subject, 
                                  cluster.include = cluster.include, 
                                  column.clusters.id = column.clusters.id)
    
    nbs_adjmat_list[[subject_id]] <- nbs_adjmat
  }
  
  return(nbs_adjmat_list)
}

#' Randomise Cluster IDs within Subject data
#' 
#' @param se.object Seurat (STUtility) object containing cluster and sample identities for each spot in the metadata.
#' @param column.clusters.id Column name in metadata corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param column.subject.id Column name in metadata corresponding to Sample ID of the spots.
#' @return New Seurat object with shuffled cluster identities per sample
#' @export
RandomiseClusteringIDs <- function (
  se.object,
  column.cluster.id,
  column.sample.id,
  random.seed = NA
) {
  if(!is.na(random.seed)){
    message(paste("Setting random seed to", random.seed))
    set.seed(random.seed)
  }
  
  #' Shuffle cluster ids for each sample
  se_metadata <- se.object@meta.data[, c(column.cluster.id, column.sample.id)]
  se_metadata$clusters_original <- se_metadata[, column.cluster.id]
  se_metadata$sample_id <- se_metadata[, column.sample.id]
  
  se_metadata_perm <- se_metadata %>% 
    dplyr::group_by(sample_id) %>% 
    dplyr::mutate(clusters_perm = clusters_original[sample(dplyr::row_number())])
  
  #' Add shuffled clusters to se object metadata
  se.object <- AddMetaData(se.object, as.character(se_metadata_perm$clusters_perm), col.name = "clusters_perm")
  
  return(se.object)
}

minmax_norm <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


# ===================================
#' RUN NBS ANALYSIS
#' 
#' "EXPECTED" VALUES: RegionNeighbours() with permuted cluster IDs for each sample
n_perm <- 50
perm_adj_mat_list <- list()
for(i in 1:n_perm){
  se <- RandomiseClusteringIDs(se.object = se, column.cluster.id = "seurat_clusters", column.sample.id = "sample_id", random.seed = i)
  se <- SetIdent(se, value = "clusters_perm")
  for(column_rm in grep(pattern = "nbs_", colnames(se[[]]), value = T)){
    se[[column_rm]] <- NULL
  }
  for(c in c_include){
    se <- RegionNeighbours(se, id = c, keep.within.id = T, verbose = TRUE)
  }
  perm_nbs_df <- se[[]][, c("clusters_perm", grep(pattern = "nbs_", colnames(se[[]]), value = T))]
  perm_adj_mat <- CreateAdjMatrix(nbs.df = perm_nbs_df, cluster.include = c_include, column.clusters.id = "clusters_perm")
  perm_adj_mat_list[[i]] <- perm_adj_mat
}

n_cols <- dim(perm_adj_mat_list[[1]])[1]
perm_adj_mat_avg <- perm_adj_mat_sd <- matrix(0L, nrow = n_cols, ncol = n_cols)
for(i in 1:n_cols){
  for(j in 1:n_cols){
    list_ij <- c()
    for(list_n in 1:length(perm_adj_mat_list)){
      list_ij <- c(list_ij, perm_adj_mat_list[[list_n]][i,j])
    }
    perm_adj_mat_avg[i,j] <- mean(list_ij)
    perm_adj_mat_sd[i,j] <- sd(list_ij)
  }
}


#' Check if random data is normally dist.
perm_adj_mat_df <- data.frame(matrix(unlist(perm_adj_mat_list), nrow=length(perm_adj_mat_list), byrow=T))  # dims: n_permX(18*18) = 50x324
perm_adj_mat_df$perm_n <- rownames(perm_adj_mat_df)
perm_adj_mat_df_long <- reshape2::melt(perm_adj_mat_df, id.vars = c("perm_n"))

png(file.path(DIR_OUT, "../figures", paste0("nbs_analysis.permutation_distribution_check.", ANALYSIS, ".png")), width = 10*300, height = 6*300, res=300)
ggplot(subset(perm_adj_mat_df_long, variable %in% paste0("X", 1:18)), aes(x=variable, y=value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free", ncol = 6)
dev.off()


#####
#' OBSERVED VALUES: Run RegionNeighbours()
se <- SetIdent(se, value = "seurat_clusters")
for(column_rm in grep(pattern = "nbs", colnames(se[[]]), value = T)){
  se[[column_rm]] <- NULL
}
for(c in c_include){
  se <- RegionNeighbours(se, id = c, keep.within.id = T, verbose = TRUE)
}

#' Create df with nbs info
nbs_df <- se[[]][, c("seurat_clusters", grep(pattern = "nbs", colnames(se[[]]), value = T))]

# write.table(nbs_df, file.path(DIR_RES, "tables", paste0("nbs_data.", ANALYSIS, ".tsv")), quote = F, sep = "\t", row.names = T)
# nbs_df <- read.table(file = file.path(DIR_RES, "tables", paste0("nbs_data.", ANALYSIS, ".tsv")), header = T, sep = "\t")

c_count <- nbs_df %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::count()
colnames(c_count) <- c("cluster", "n")
c_count$id <- paste0("C", c_count$cluster)

c_nodes <- cbind(as.data.frame(c_count[c_include, ]), c_anno[order(c_anno$seurat_clusters), ][c_include, c("cluster_anno", "cluster_group", "cluster_color")])
rownames(c_nodes) <- c_nodes$id
c_nodes <- c_nodes[, c("id", "cluster", "cluster_anno", "cluster_group", "cluster_color", "n")]


#' Adjacency matrix
nbs_adjmat <- CreateAdjMatrix(nbs.df = nbs_df, 
                              cluster.include = c_include, 
                              column.clusters.id = "seurat_clusters"
                              )


#' Normalize method 1
nbs_adjmat_norm <- nbs_adjmat
for (i in seq(1:dim(nbs_adjmat_norm)[1])) {
  for (j in seq(1:dim(nbs_adjmat_norm)[1])) {
    e_sum <- sum(nbs_adjmat_norm[i,i], nbs_adjmat_norm[j,j])
    if (i!=j) {
      nbs_adjmat_norm[i,j] <- round(nbs_adjmat_norm[i,j]/e_sum, 4)*100
    }
  }
}


#' "Normalize"/Standardise method 2: Using permuted values
nbs_adjmat_permscore <- round( ((nbs_adjmat - perm_adj_mat_avg) / perm_adj_mat_sd) , digits = 3)
diag(nbs_adjmat_permscore) <- diag(nbs_adjmat)


#' Set row/column names
rownames(nbs_adjmat) <- rownames(nbs_adjmat_norm) <- rownames(nbs_adjmat_permscore) <- paste0("C", c_include)
colnames(nbs_adjmat) <- colnames(nbs_adjmat_norm) <- colnames(nbs_adjmat_permscore) <- paste0("C", c_include)

rownames(perm_adj_mat_avg) <- rownames(perm_adj_mat_sd) <- colnames(perm_adj_mat_avg) <- colnames(perm_adj_mat_sd) <- paste0("C", c_include)


#' Save tables
write.table(nbs_adjmat, file.path(DIR_RES, "tables", paste0("nbs_adj_matrix.", ANALYSIS, ".tsv")), quote = F, sep = "\t", row.names = T, col.names = T)
write.table(nbs_adjmat_norm, file.path(DIR_RES, "tables", paste0("nbs_adj_matrix_norm.", ANALYSIS, ".tsv")), quote = F, sep = "\t", row.names = T, col.names = T)
write.table(nbs_adjmat_permscore, file.path(DIR_RES, "tables", paste0("nbs_adj_matrix_permscore.", ANALYSIS, ".tsv")), quote = F, sep = "\t", row.names = T, col.names = T)

write.table(perm_adj_mat_avg, file.path(DIR_RES, "tables", paste0("nbs_adj_matrix_perm_avg.", ANALYSIS, ".tsv")), quote = F, sep = "\t", row.names = T, col.names = T)
write.table(perm_adj_mat_sd, file.path(DIR_RES, "tables", paste0("nbs_adj_matrix_perm_sd.", ANALYSIS, ".tsv")), quote = F, sep = "\t", row.names = T, col.names = T)



#####
#' PER SUBJECT, OBSERVED VALUES

#' Adjacency matrix
nbs_adjmat_list_subject <- CreateAdjMatrixPerSubject(nbs.df = nbs_df, 
                                                     se.metadata = se@meta.data,
                                                     column.subject.id = "subject_id",
                                                     cluster.include = c_include, 
                                                     column.clusters.id = "seurat_clusters"
                                                     )

#' PER SUBJECT, EXPECTED VALUES
n_perm <- 50
perm_adj_mat_list_subject <- list()
for(i in 1:n_perm){
  se <- RandomiseClusteringIDs(se.object = se, column.cluster.id = "seurat_clusters", column.sample.id = "sample_id", random.seed = i)
  se <- SetIdent(se, value = "clusters_perm")
  for(column_rm in grep(pattern = "nbs_", colnames(se[[]]), value = T)){
    se[[column_rm]] <- NULL
  }
  for(c in c_include){
    se <- RegionNeighbours(se, id = c, keep.within.id = T, verbose = TRUE)
  }
  perm_nbs_df <- se[[]][, c("clusters_perm", grep(pattern = "nbs_", colnames(se[[]]), value = T))]
  perm_adj_mat_list_subject_id <- CreateAdjMatrixPerSubject(nbs.df = perm_nbs_df, 
                                                           se.metadata = se@meta.data, 
                                                           cluster.include = c_include, 
                                                           column.clusters.id = "clusters_perm", 
                                                           column.subject.id = "subject_id")
  perm_adj_mat_list_subject[[i]] <- perm_adj_mat_list_subject_id
}


n_cols <- dim(perm_adj_mat_list_subject[[1]][[1]])[1]
perm_adj_mat_avg_list_subject <- perm_adj_mat_sd_list_subject <- list()
perm_adj_mat_avg_subject <- perm_adj_mat_sd_subject <- matrix(0L, nrow = n_cols, ncol = n_cols)

for (subject_id in names(perm_adj_mat_list_subject[[1]])) {
  perm_adj_mat_list_subject_id <- list()
  for (i in 1:length(perm_adj_mat_list_subject)) {
    perm_adj_mat_list_subject_id[[i]] <- perm_adj_mat_list_subject[[i]][[subject_id]]
  }
  for (i in 1:n_cols) {
    for (j in 1:n_cols) {
      list_ij <- c()
      for (list_n in 1:length(perm_adj_mat_list_subject_id)) {
        list_ij <- c(list_ij, perm_adj_mat_list_subject_id[[list_n]][i,j])
      }
      perm_adj_mat_avg_subject[i,j] <- mean(list_ij)
      perm_adj_mat_sd_subject[i,j] <- sd(list_ij)
    }
    perm_adj_mat_avg_list_subject[[subject_id]] <- perm_adj_mat_avg_subject
    perm_adj_mat_sd_list_subject[[subject_id]] <- perm_adj_mat_sd_subject
  }
}


#' PER SUBJECT, OBS - EXP_AVG
nbs_adjmat_permscore_subject <- list()
for (subject_id in names(nbs_adjmat_list_subject)) {
  nbs_adjmat_permscore_subject[[subject_id]] <- round( (nbs_adjmat_list_subject[[subject_id]] - perm_adj_mat_avg_list_subject[[subject_id]]) / perm_adj_mat_sd_list_subject[[subject_id]], digits = 3)
  diag(nbs_adjmat_permscore_subject[[subject_id]]) <- diag(nbs_adjmat_list_subject[[subject_id]])

  #' Set row/column names
  rownames(nbs_adjmat_permscore_subject[[subject_id]]) <- colnames(nbs_adjmat_permscore_subject[[subject_id]]) <- paste0("C", c_include)
}


#' Save tables
for (subject_id in names(nbs_adjmat_permscore_subject)) {
  write.table(nbs_adjmat_permscore_subject[[subject_id]], 
              file.path(DIR_RES, "tables/nbs_kavg_subjects", paste0("nbs_adj_matrix_permscore_subject-", subject_id, ".", ANALYSIS, ".tsv")), 
              quote = F, sep = "\t", row.names = T, col.names = T)
  write.table(nbs_adjmat_list_subject[[subject_id]], 
              file.path(DIR_RES, "tables/nbs_kavg_subjects", paste0("nbs_adj_matrix_subject-", subject_id, ".", ANALYSIS, ".tsv")), 
              quote = F, sep = "\t", row.names = T, col.names = T)
  write.table(perm_adj_mat_avg_list_subject[[subject_id]], 
              file.path(DIR_RES, "tables/nbs_kavg_subjects", paste0("nbs_adj_matrix_perm_avg_subject-", subject_id, ".", ANALYSIS, ".tsv")), 
              quote = F, sep = "\t", row.names = T, col.names = T)
  write.table(perm_adj_mat_sd_list_subject[[subject_id]], 
              file.path(DIR_RES, "tables/nbs_kavg_subjects", paste0("nbs_adj_matrix_perm_sd_subject-", subject_id, ".", ANALYSIS, ".tsv")), 
              quote = F, sep = "\t", row.names = T, col.names = T)
}

# ===================================
#' PLOTS
color_low2 <- "#62376e"

#' PLOT 1
#' Make graph
nbs_adjmat_df <- nbs_adjmat_norm
g <- graph.adjacency(nbs_adjmat_df, mode = "undirected", weighted = TRUE, diag = F)

e_sum <- diag(nbs_adjmat_df)
df <- data.frame(name = names(e_sum),
                 size = e_sum,
                 norm_size = e_sum/max(e_sum)*40, 
                 group = subset(c_nodes, id %in% names(e_sum))$cluster_group,
                 color = subset(c_nodes, id %in% names(e_sum))$cluster_color,
                 c_size = subset(c_nodes, id %in% names(e_sum))$n,
                 stringsAsFactors = F)
links <- data.frame(
  source=as_edgelist(g, names = TRUE)[,1],
  target=as_edgelist(g, names = TRUE)[,2])

g_df_norm <- data.frame(links,
                        weight = E(g)$weight,
                        weight_minmax = minmax_norm(abs(E(g)$weight)))

g2 <- graph_from_data_frame(d = links, vertices = df, directed = F)

g2 <- set_edge_attr(g2, "weight", value=E(g)$weight)
E(g2)$width <- E(g2)$weight*1.5

#' plot 1.1
# fname <- paste0("nbs_analysis.graph.", ANALYSIS)
# pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 5.5, height = 5.5, useDingbats = F)
# plot(g2, 
#      layout=layout_in_circle,
#      vertex.label.family = "Helvetica", vertex.label.color = "black", vertex.label.cex = 1.5,
#      vertex.size = V(g2)$norm_size,
#      vertex.color = adjustcolor(V(g2)$color, alpha.f=.9),
#      vertex.frame.color="white",
#      edge.color = adjustcolor("grey70", alpha.f = .5),
#      edge.curved=0.1)
# 
# for(cluster_plot in names(e_sum)){
#   e_pairs <- c()
#   for(c in names(e_sum)[!names(e_sum)==cluster_plot]){
#     e_pairs <- c(e_pairs, cluster_plot, c)
#   }
#   ecol <- rep("gray95", ecount(g2))
#   ecol[get.edge.ids(g, vp = e_pairs)] <- "gray20"
#   
#   plot(g2, 
#        layout=layout_in_circle,
#        vertex.label.family = "Helvetica", vertex.label.color = "black", vertex.label.cex = 1.5,
#        vertex.size = V(g2)$norm_size,
#        vertex.color = adjustcolor(V(g2)$color, alpha.f=.9),
#        vertex.frame.color="white",
#        edge.color = adjustcolor(ecol, alpha.f = .25),
#        edge.curved=0.1)
# }
# dev.off()

#' plot 1.2
# fname <- paste0("nbs_analysis.graph2.", ANALYSIS)
# pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 5.5, height = 5.5, useDingbats = F)
# for(cluster_plot in names(e_sum)){
#   e_pairs <- c()
#   for(c in names(e_sum)[!names(e_sum)==cluster_plot]){
#     e_pairs <- c(e_pairs, cluster_plot, c)
#   }
#   ecol <- rep(NA, ecount(g2))
#   ecol[get.edge.ids(g, vp = e_pairs)] <- "gray50"
#   
#   node_ids <- names(V(g3))
#   node_sizes <- data.frame(node = node_ids,
#                            size = V(g3)$norm_size, 
#                            stringsAsFactors = F)
#   LS <- layout_as_star(g3, center = cluster_plot, order = order(node_sizes$size, decreasing = T))
#   
#   plot(g2, 
#        layout=LS,
#        vertex.label.family = "Helvetica", vertex.label.color = "black", vertex.label.cex = 1.5,
#        vertex.size = V(g2)$norm_size,
#        vertex.color = adjustcolor(V(g2)$color, alpha.f=.9),
#        vertex.frame.color="white",
#        edge.color = adjustcolor(ecol, alpha.f = .5),
#        edge.curved=0)
# }
# dev.off()


#' plot 1.3
# fname <- paste0("nbs_analysis.graph3.", ANALYSIS)
# pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 5.5, height = 5.5, useDingbats = F)
# g3 <- g2
# for(cluster_plot in names(e_sum)){
#   e_pairs <- c()
#   for(c in names(e_sum)[!names(e_sum)==cluster_plot]){
#     e_pairs <- c(e_pairs, cluster_plot, c)
#   }
#   ecol <- rep(NA, ecount(g3))
#   ecol[get.edge.ids(g, vp = e_pairs)] <- "gray50"
#   
#   node_ids <- names(V(g3))
#   node_widths <- E(g3)$weight[get.edge.ids(g3, vp = e_pairs)]
#   node_sizes <- data.frame(node = node_ids,
#                            size = V(g3)$norm_size,
#                            stringsAsFactors = F)
#   node_sizes$width <- 0
#   node_sizes[node_sizes$node!=cluster_plot, "width"] <- node_widths
#   LS <- layout_as_star(g3, center = cluster_plot, order = order(node_sizes$width, decreasing = T))
#   
#   x <- E(g3)$weight[get.edge.ids(g3, vp = e_pairs)]
#   norm_x = (x-min(x))/(max(x)-min(x))
#   norm_x <- (norm_x+.1)*20
#   E(g3)$width[get.edge.ids(g3, vp = e_pairs)] <- norm_x
#   
#   plot(g3, 
#        layout=LS,
#        vertex.label.family = "Helvetica", vertex.label.color = "black", vertex.label.cex = 1.5,
#        vertex.size = V(g3)$norm_size,
#        vertex.color = adjustcolor(V(g3)$color, alpha.f=.9),
#        vertex.frame.color="white",
#        edge.color = adjustcolor(ecol, alpha.f = .5),
#        edge.curved=0)
# }
# dev.off()


#' PLOT 2: cluster in center
g3 <- g2
E(g3)$width <- E(g2)$weight*4

for(cluster_plot in c("C3", "C4", "C6")){
  e_pairs <- c()
  for(c in names(e_sum)[!names(e_sum)==clusters_plot]){
    e_pairs <- c(e_pairs, cluster_plot, c)
  }
  ecol <- rep(NA, ecount(g3))
  ecol[get.edge.ids(g, vp = e_pairs)] <- "gray50"
  
  node_ids <- names(V(g3))
  node_sizes <- data.frame(node = node_ids,
                           size = V(g3)$norm_size, 
                           stringsAsFactors = F)
  LS <- layout_as_star(g3, center = cluster_plot, order = order(node_sizes$size, decreasing = T))

  fname <- paste0("nbs_analysis.graph_adipocytes_", cluster_plot, ".", ANALYSIS)
  cairo_ps(filename = file.path(DIR_FIG, paste0(fname, ".eps")), width = 5.5, height = 5.5)
  plot(g3, 
       layout=LS,
       vertex.label.family = "Helvetica", 
       vertex.label.color = "black", 
       vertex.label.cex = 1.5,
       vertex.size = V(g3)$norm_size*1.25,
       vertex.color = adjustcolor(V(g3)$color, alpha.f=1),
       vertex.frame.color="white",
       edge.color = adjustcolor(ecol, alpha.f = .5),
       edge.curved=0)
  dev.off()
}


#' PLOT 3: Permuted score values
#' Make graph
nbs_adjmat_df <- nbs_adjmat_permscore
g <- graph.adjacency(nbs_adjmat_df, mode = "undirected", weighted = TRUE, diag = F)

e_sum <- diag(nbs_adjmat_df)
df <- data.frame(name = names(e_sum),
                 size = e_sum,
                 norm_size = e_sum/max(e_sum)*40, 
                 group = subset(c_nodes, id %in% names(e_sum))$cluster_group,
                 color = subset(c_nodes, id %in% names(e_sum))$cluster_color,
                 c_size = subset(c_nodes, id %in% names(e_sum))$n,
                 stringsAsFactors = F)
links <- data.frame(
  source=as_edgelist(g, names = TRUE)[,1],
  target=as_edgelist(g, names = TRUE)[,2])

g_df_permscore <- data.frame(links,
                             weight = E(g)$weight,
                             weight_minmax = minmax_norm(abs(E(g)$weight)))


g2 <- graph_from_data_frame(d = links, vertices = df, directed = F)

g2 <- set_edge_attr(g2, "weight", value = minmax_norm(abs(E(g)$weight)))
E(g2)$width <- (E(g2)$weight+0.1)*14


fname <- paste0("nbs_analysis.permscore.", ANALYSIS)
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 5.5, height = 5.5, useDingbats = F)
for(cluster_plot in names(e_sum)){
  e_pairs <- c()
  for(c in names(e_sum)[!names(e_sum)==cluster_plot]){
    e_pairs <- c(e_pairs, cluster_plot, c)
  }
  e_weights <- E(g)$weight
  ecol <- rep("grey90", ecount(g2))
  ecol[get.edge.ids(g2, vp = e_pairs)] <- "grey20"
  ecol[e_weights>0 & ecol=="grey20"] <- color_high
  ecol[e_weights<0 & ecol=="grey20"] <- color_low2
  ecol[ecol=="grey90"] <- NA
  
  node_ids <- names(V(g2))
  node_e_widths <- E(g2)$width[get.edge.ids(g2, vp = e_pairs)]
  node_sizes <- data.frame(node = node_ids,
                           size = V(g2)$norm_size, 
                           # width = c(node_e_widths[1:match(cluster_plot, node_ids)-1], 0, node_e_widths[match(cluster_plot, node_ids):length(node_e_widths)]),
                           stringsAsFactors = F)
  LS <- layout_as_star(g2, center = cluster_plot, order = order(node_sizes$size, decreasing = T))
  plot(g2, 
       layout=LS,
       vertex.label.family = "Helvetica", vertex.label.color = "black", vertex.label.cex = 1.5,
       vertex.size = V(g2)$norm_size,
       vertex.color = adjustcolor(V(g2)$color, alpha.f=1),
       vertex.frame.color="white",
       edge.color = adjustcolor(ecol, alpha.f = .8),
       edge.curved=0)
}
dev.off()


#' Export tables for plots
write.csv(g_df_norm, file.path(DIR_RES, "tables", paste0("nbs_graph_df_norm.", ANALYSIS, ".csv")), quote = F, row.names = F)
write.csv(g_df_permscore, file.path(DIR_RES, "tables", paste0("nbs_graph_df_permscore.", ANALYSIS, ".csv")), quote = F, row.names = F)


#' PLOT HEATMAP
hm_df <- nbs_adjmat_permscore
diag(hm_df) <- NA
hm_df[lower.tri(hm_df)] <- NA
hm_df <- hm_df[,rev(colnames(hm_df))]

breaksList <- seq(-6, 6, by = .5)
p1 <- pheatmap(hm_df,
              cluster_rows = F, cluster_cols = F,
              breaks = breaksList,
              color = colorRampPalette(c(color_low3, "white", color_high2))(length(breaksList)-1),
              border_color = "white",
              na_col = "white")

fname <- paste0("nbs_analysis.permscore_hm.", ANALYSIS)
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 5, height = 4.8, useDingbats = F);p1;dev.off()


#' PER SUBEJECT, PLOT HEATMAP
fname <- paste0("nbs_analysis.permscore_subject_hm.", ANALYSIS)
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 5, height = 5, useDingbats = F)
for (subject_id in names(nbs_adjmat_permscore_subject)) {
  hm_df <- nbs_adjmat_permscore_subject[[subject_id]]
  
  diag(hm_df) <- NA
  hm_df[lower.tri(hm_df)] <- NA
  hm_df[is.infinite(hm_df)] <- NA
  hm_df <- hm_df[,rev(colnames(hm_df))]
  breaksList <- seq(-6, 6, by = .5)
  pheatmap(hm_df,
           cluster_rows = F, 
           cluster_cols = F,
           breaks = breaksList,
           color = colorRampPalette(c(color_low3, "white", color_high2))(length(breaksList)-1),
           border_color = "white",
           na_col = "grey80",
           main  = subject_id)
}
dev.off()


#===============
sessionInfo()
