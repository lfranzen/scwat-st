#' Cluster neighbourhoo analysis
#' Idnetify cluster neighbours and plot as graph 

library(ggplot2)
library(magrittr)
library(STutility)

library(igraph)
library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggraph)


PROJECT_ID <- "visium"
DIR_ROOT <- getwd()
DIR_WD <- file.path(getwd(), "scripts")  # 
DIR_DATA <- file.path(DIR_ROOT, "data", PROJECT_ID)
DIR_RES <- file.path(DIR_ROOT, "results" , PROJECT_ID)
DIR_FIG <- file.path(DIR_RES, "figures")

source(file.path(DIR_WD, "colors.R"))


# ===================================
#' READ DATA & LOAD IMGS
#' 'Baseline' Visium data
canno_base <- read.csv(file.path(DIR_RES, "tables", "visium_baseline.clustering_annotations.csv"), stringsAsFactors = F)

se_base <- readRDS(file.path(DIR_RES, "se-object.visium_baseline.rds"))
se_base <- LoadImages(se_base, time.resolve = T, verbose = T, xdim = 100)


#' 'Insulin' Visium data
canno_ins <- read.csv(file.path(DIR_RES, "tables", "visium_insulin.clustering_annotations.csv"), stringsAsFactors = F)

se_ins <- readRDS(file.path(DIR_RES, "se-object.visium_insulin.rds"))
se_ins <- LoadImages(se_ins, time.resolve = T, verbose = T, xdim = 100)


# ===================================
#' RUN NBS ANALYSIS

# ANALYSIS <- "baseline"
ANALYSIS <- "insulin"


if(ANALYSIS=="baseline"){
  se <- se_base
  c_anno <- canno_base
  c_include <- seq(3, max(as.numeric(se[[]]$seurat_clusters)), 1)
  } else if (ANALYSIS=="insulin") {
  se <- se_ins
  c_anno <- canno_ins
  c_include <- seq(4, max(as.numeric(se[[]]$seurat_clusters)), 1)
}

#' Run RegionNeighbours()
se <- SetIdent(se, value = "seurat_clusters")
for(column_rm in grep(pattern = "nbs", colnames(se[[]]), value = T)){
  se[[column_rm]] <- NULL
}
for(c in c_include){
  se <- RegionNeighbours(se, id = c, keep.within.id = T, verbose = TRUE)
}

#' Create df with nbs info
nbs_df <- se[[]][, c("seurat_clusters", grep(pattern = "nbs", colnames(se[[]]), value = T))]
write.table(nbs_df, file.path(DIR_RES, "tables", paste0("nbs_data.", ANALYSIS, ".tsv")), quote = F, sep = "\t", row.names = T)

# nbs_df <- read.table(file = file.path(DIR_RES, "tables", paste0("nbs_data.", ANALYSIS, ".tsv")), header = T, sep = "\t")


c_count <- nbs_df %>%
  group_by(seurat_clusters) %>%
  dplyr::count()
colnames(c_count) <- c("cluster", "n")
c_count$id <- paste0("C", c_count$cluster)

c_nodes <- cbind(as.data.frame(c_count[c_include, ]), c_anno[order(c_anno$seurat_clusters), ][c_include, c("cluster_anno", "cluster_group", "cluster_color")])
rownames(c_nodes) <- c_nodes$id
c_nodes <- c_nodes[, c("id", "cluster", "cluster_anno", "cluster_group", "cluster_color", "n")]


#' Adjacency matrix
nbs_adjmat <- matrix(0L, nrow = length(c_include), ncol = length(c_include))
for (i in seq_along(c_include)) {
  c <- c_include[i]
  c_nbs_df <- nbs_df[nbs_df$seurat_clusters==c, ]
  for (j in seq_along(c_include)) {
    nbs <- paste0("nbs_", c_include[j])
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

nbs_adjmat_norm <- nbs_adjmat
for (i in seq(1:dim(nbs_adjmat_norm)[1])) {
  for (j in seq(1:dim(nbs_adjmat_norm)[1])) {
    e_sum <- sum(nbs_adjmat_norm[i,i], nbs_adjmat_norm[j,j])
    if (i!=j) {
      nbs_adjmat_norm[i,j] <- round(nbs_adjmat_norm[i,j]/e_sum, 4)*100
    }
  }
}


rownames(nbs_adjmat) <- rownames(nbs_adjmat_norm) <- paste0("C", c_include)
colnames(nbs_adjmat) <- colnames(nbs_adjmat_norm) <- paste0("C", c_include)

write.table(nbs_adjmat, file.path(DIR_RES, "tables", paste0("nbs_adj_matrix.", ANALYSIS, ".tsv")), quote = F, sep = "\t", row.names = T, col.names = T)
write.table(nbs_adjmat_norm, file.path(DIR_RES, "tables", paste0("nbs_adj_matrix_norm.", ANALYSIS, ".tsv")), quote = F, sep = "\t", row.names = T, col.names = T)

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

g2 <- graph_from_data_frame(d = links, vertices = df, directed = F)

g2 <- set_edge_attr(g2, "weight", value=E(g)$weight)
E(g2)$width <- E(g2)$weight*1.5


# ===================================
#' PLOT 
fname <- paste0("nbs_analysis.graph.", ANALYSIS)
pdf(file = file.path(DIR_FIG, paste0(fname, ".pdf")), width = 5.5, height = 5.5)
plot(g2, 
     layout=layout_in_circle,
     vertex.label.family = "Helvetica", vertex.label.color = "black", vertex.label.cex = 1.5,
     vertex.size = V(g2)$norm_size,
     vertex.color = adjustcolor(V(g2)$color, alpha.f=.9),
     vertex.frame.color="white",
     edge.color = adjustcolor("grey70", alpha.f = .5),
     edge.curved=0.1)

for(cluster_plot in names(e_sum)){
  e_pairs <- c()
  for(c in names(e_sum)[!names(e_sum)==cluster_plot]){
    e_pairs <- c(e_pairs, cluster_plot, c)
  }
  ecol <- rep("gray95", ecount(g2))
  ecol[get.edge.ids(g, vp = e_pairs)] <- "gray20"
  
  plot(g2, 
       layout=layout_in_circle,
       vertex.label.family = "Helvetica", vertex.label.color = "black", vertex.label.cex = 1.5,
       vertex.size = V(g2)$norm_size,
       vertex.color = adjustcolor(V(g2)$color, alpha.f=.9),
       vertex.frame.color="white",
       edge.color = adjustcolor(ecol, alpha.f = .25),
       edge.curved=0.1)
}
dev.off()




