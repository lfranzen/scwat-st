#'====================================================
#' Run neighborhood analyses on public 10x Visium example data
#'
#'====================================================

ANALYSIS_ID_1 <- "10x_HuLymphNode"
ANALYSIS_ID_2 <- "10x_HuCerebellum"
ANALYSIS_ID <- ANALYSIS_ID_2

DIR_ROOT <- file.path(getwd())  # , ".."
DIR_RES <- file.path(DIR_ROOT, "results", "10x_visium_data")

#'====================================================
#' PREP
#' Read functions
source(file = file.path(DIR_ROOT, "scripts", "nbs_analysis_functions.R"))

#' Read se object
se <- readRDS(file = file.path(DIR_RES, paste0("se_obj.", ANALYSIS_ID, ".rds")))
se <- LoadImages(se, xdim = 100)


#'====================================================
#' CALCULATE HOMOTYPIC SCORES
obs_k_avg_df <- CalculateAvgDegree(se.object = se, column.cluster.id = "Cluster", column.sample.id = "orig.ident")
exp_k_avg_df_list <- RandomClustersAvgDegree(se.object = se, column.cluster.id = "Cluster", column.sample.id = "orig.ident", n.perm = 50)
homoscore_df <- CalculateHomotypicScore(avk.k.observed = obs_k_avg_df, avk.k.expected = exp_k_avg_df_list);homoscore_df

#' Export table
write.csv(homoscore_df, file = file.path(DIR_RES, paste0(ANALYSIS_ID, "_homotypic_score.csv")), row.names = T)

#' Plot
d_plot <- as.data.frame(t(homoscore_df))
d_plot$cluster <- paste0("Cluster_", as.character(sub(pattern = "kavg_cluster_", "", rownames(d_plot))))

p <- ggplot(d_plot, aes(x=reorder(cluster, S1), y=S1)) +
  geom_col(fill="black") +
  labs(x="", y="Homotypic score", title = paste0(ANALYSIS_ID, ": Average degree obs-exp")) +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(color="black", size=12, face = "bold"),
        axis.text = element_text(color="black", size=12),
        axis.title = element_text(color="black", size=12));p

res_use <- 300
png(filename = file.path(DIR_RES, paste0(ANALYSIS_ID, "_homotypic_score.png")), width = 5*res_use, height = 3*res_use, res = res_use);p;dev.off()
pdf(file = file.path(DIR_RES, paste0(ANALYSIS_ID, "_homotypic_score.pdf")), width = 5, height = 3, useDingbats = F);p;dev.off()


#'====================================================
#' CALCULATE HETEROTYPIC SCORES

#' Observed: adjacency matrix
se <- SetIdent(se, value = "Cluster")
for(column_rm in grep(pattern = "nbs", colnames(se[[]]), value = T)){se[[column_rm]] <- NULL}

for(c in 1:as.numeric(max(se$Cluster))){
  se <- STutility::RegionNeighbours(se, id = c, keep.within.id = T, verbose = TRUE)
  }
obs_nbs_df <- se@meta.data[, c("Cluster", grep(pattern = "nbs", colnames(se@meta.data), value = T))]

obs_adjmat <- CreateAdjMatrix(nbs.df = obs_nbs_df, column.cluster.id = "Cluster")


#' Expected (50 permutations): adjacency matrix
perm_adj_mat_list <- RandomClustersAdjMat(se.object = se, 
                                          column.cluster.id = "Cluster", 
                                          column.sample.id = "orig.ident", 
                                          n.perm = 50)

perm_adj_mat_obj <- CalulateMeanSdAdjMatList(adj.mat.list = perm_adj_mat_list)


#' Compute score
adj_mat_score <- CalculateHeterotypicScore(adj.mat.observed = obs_adjmat, 
                                           adj.mat.expected.mean = perm_adj_mat_obj@mean, 
                                           adj.mat.expected.sd = perm_adj_mat_obj@sd)

#' Export table
write.csv(adj_mat_score, file = file.path(DIR_RES, paste0(ANALYSIS_ID, "_heterotypic_score.csv")), row.names = T)

#' Plot - heatmap
library(pheatmap)

hm_df <- adj_mat_score
hm_df <- read.csv(file = file.path(DIR_RES, paste0(ANALYSIS_ID, "_heterotypic_score.csv")), header = T, row.names = 1)
diag(hm_df) <- NA
# hm_df[lower.tri(hm_df)] <- NA
hm_df <- hm_df[,rev(colnames(hm_df))]

# breaksList <- seq(-6, 6, by = .5)
# breaksList <- seq(-max(abs(hm_df), na.rm = T),
#                   max(abs(hm_df), na.rm = T),
#                   by = 1)
breaksList <- seq(-20, 20, by = 1)

p <- pheatmap(hm_df,
             cluster_rows = F, 
             cluster_cols = F,
             main = paste0(ANALYSIS_ID, ": Heterotypic score"),
             breaks = breaksList, 
             cellheight = 30, cellwidth = 30,
             color = colorRampPalette(c("#772A82", "white", "#1D793D"))(length(breaksList)-1),
             border_color = "white",
             na_col = "grey50")

res_use <- 300
png(filename = file.path(DIR_RES, paste0(ANALYSIS_ID, "_heterotypic_score_hm.png")), width = 4.5*res_use, height = 4.5*res_use, res = res_use);p;dev.off()
pdf(file = file.path(DIR_RES, paste0(ANALYSIS_ID, "_heterotypic_score_hm.pdf")), width = 4.5, height = 4.5, useDingbats = F);p;dev.off()

#'====================================================
