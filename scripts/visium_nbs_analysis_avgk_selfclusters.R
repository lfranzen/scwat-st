########################
#' Visium neighbourhood analysis of "self-clustering"
########################

#' Set up
library(magrittr)
library(STutility)

DIR_IN <- file.path(getwd(), "data", "visium")
DIR_OUT <- file.path(getwd(), "results", "visium", "tables")


#####
#' Read data
metadata <- read.delim(file.path(DIR_IN, "visium_sample_metadata.tsv"), sep = "\t", stringsAsFactors = F)

se_base <- readRDS(file.path(DIR_OUT, "..", "se-object.visium_baseline.rds"))
se <- se_base
se <- LoadImages(se, time.resolve = T, verbose = T, xdim = 100)

se_stats1 <- se@meta.data %>%
  dplyr::group_by(subject_alias, subject_id, tissue_id, seu_n, novaseq_id) %>% 
  dplyr::count(name = "spots")
se_stats1$seu_n <- paste0("S", se_stats1$seu_n)
colnames(se_stats1)[colnames(se_stats1)=="seu_n"] <- "sample"


#####
#' Create "spatial" network for each sample
spatnet_init <- GetSpatNet(se)

spatnet <- do.call(rbind, lapply(seq_along(spatnet_init), function(i) {
  spnet <- spatnet_init[[i]]
  spnet$cluster_from <- se[[]][spnet$from, "seurat_clusters"]
  spnet$cluster_to <- se[[]][spnet$to, "seurat_clusters"]
  spnet$sample <- paste0(i)
  return(spnet)
}))

spatnet$cluster_from <- as.factor(as.numeric(spatnet$cluster_from))
spatnet$cluster_to <- as.factor(as.numeric(spatnet$cluster_to))


#####
#' Plot
colors_clusters <- se@meta.data[, c("seurat_clusters", "cluster_anno", "cluster_group", "cluster_color")] %>%
  dplyr::distinct() %>%
  dplyr::arrange(seurat_clusters)

pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_k.baseline2.pdf")), width = 12, height = 9)
for(cluster_view in c("1", "3", "4", "6", "8", "16", "20")){
  spatnet.subset <- subset(spatnet, cluster_from %in% cluster_view)
  spatnet.subset_conly <- subset(spatnet.subset, cluster_to %in% cluster_view)

  p <- ggplot() +
    geom_point(data = spatnet.subset, aes(start_x, -start_y, color = cluster_from), size = .5) +
    geom_segment(data = spatnet.subset_conly, aes(x = start_x, xend = end_x, y = -start_y, yend = -end_y), size=0.3) +
    labs(color="") +
    scale_color_manual(values = colors_clusters[colors_clusters$seurat_clusters==cluster_view,"cluster_color"]) +
    facet_wrap(~sample) +
    theme_classic();
  print(p)
  # png(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_k_",cluster_view,".png")), width = 6000, height = 4500, res = 500);p;dev.off()
}
dev.off()


#####
#' Count avgerage degree (k_avg)
N <- length(unique(c(spatnet.subset$from)))
ki <- spatnet.subset_conly %>%
  dplyr::group_by(from) %>%
  dplyr::count(name = "ki")
L <- sum(ki$ki)/2
k_avg <- (2*L)/N


#' Calculate avg degree for all samplse and clusters
k_avg_list <- list()
k_avg_df <- data.frame(row.names = paste0("S", seq(1,10)))

for(cluster_view in as.character(seq(1,20))){
  for(s in seq(1,10)){
    spatnet.subset <- subset(spatnet, cluster_from %in% cluster_view)
    spatnet.subset <- subset(spatnet.subset, sample == s)
    spatnet.subset_conly <- subset(spatnet.subset, cluster_to %in% cluster_view)
    
    N <- length(unique(c(spatnet.subset$from)))
    ki <- spatnet.subset_conly %>%
      dplyr::group_by(from) %>%
      dplyr::count(name = "ki")
    L <- sum(ki$ki)/2
    k_avg <- (2*L)/N
    k_avg_list[[paste0("S", s)]] <- k_avg
  }
  k_avg_df_add <- as.data.frame(unlist(k_avg_list))
  colnames(k_avg_df_add) <- paste0("kavg_cluster_", cluster_view)
  k_avg_df <- cbind(k_avg_df, k_avg_df_add)
}

k_avg_df$sample <- rownames(k_avg_df)


#' Reorder sample ids to match correct tissue
sample_conv <- t(data.frame(S1 = "S42",
                            S2 = "S55",
                            S3 = "S44",
                            S4 = "S46",
                            S5 = "S48",
                            S6 = "S49",
                            S7 = "S50",
                            S8 = "S51",
                            S9 = "S52",
                            S10 = "S54", 
                            stringsAsFactors = F, row.names = "novaseq_id"))
sample_conv <- as.data.frame(sample_conv)
sample_conv$sample <- rownames(sample_conv)
se_stats_neworder <- se_stats1
se_stats_neworder$sample <- NULL
se_stats_neworder <- merge(se_stats_neworder, sample_conv, by="novaseq_id")
k_avg_out <- merge(se_stats_neworder[,c(1:4,6)], k_avg_df, by = "sample")


write.csv(k_avg_out, file = file.path(DIR_OUT, "nbs_cluster_kavg.baseline.csv"), row.names = F)
# write.csv(k_avg_out[,1:5], file = file.path(DIR_OUT, "nbs_cluster_kavg.sample_id_conversion.csv"), row.names = F)


#####
#' Randomise spot cluster locations - permute multiple times to obtain avg k for random state

RandomClusteringSpatNet <- function (
  se.object,
  column.cluster.id,
  column.sample.id,
  random.seed = NA,
  se.SpatNet = NA
) {
  if(!is.na(random.seed)){
    message(paste("Setting random seed to", random.seed))
    set.seed(random.seed)
  } else {
    rm(.Random.seed, envir=globalenv())
  }
  
  #' Shuffle cluster ids for each sample
  se_metadata <- se.object@meta.data[, c(column.cluster.id, column.sample.id)]
  se_metadata$clusters_original <- se_metadata[, column.cluster.id]
  se_metadata$sample_id <- se_metadata[, column.sample.id]
  
  se_metadata_perm <- se_metadata %>% 
    dplyr::group_by(sample_id) %>% 
    dplyr::mutate(clusters_perm = seurat_clusters[sample(dplyr::row_number())])

  #' Add shuffled clusters to se object metadata
  se.object <- AddMetaData(se.object, as.character(se_metadata_perm$clusters_perm), col.name = "clusters_perm")
  
  #' Create SpatNet and make network for permuted clusters
  if(is.list(se.SpatNet)){
    spatnet_init <- se.SpatNet
  } else {
    message("Creating network...")
    spatnet_init <- GetSpatNet(se.object)
  }
  
  spatnet <- do.call(rbind, lapply(seq_along(spatnet_init), function(i) {
    spnet <- spatnet_init[[i]]
    spnet$cluster_from <- se.object[[]][spnet$from, "clusters_perm"]
    spnet$cluster_to <- se.object[[]][spnet$to, "clusters_perm"]
    spnet$sample <- paste0(i)
    return(spnet)
  }))
  spatnet$cluster_from <- as.character(spatnet$cluster_from)
  spatnet$cluster_to <- as.character(spatnet$cluster_to)
  
  #' Calculate avgerage degree for all samplse and clusters
  message("Calculating average degree for each cluster and sample...")
  sample_ids <- seq(1, length(unique(se_metadata_perm$sample_id)))
  cluster_ids <- sort(as.character(unique(se_metadata_perm$clusters_perm)))
  k_avg_list <- list()
  k_avg_df <- data.frame(row.names = paste0("S",sample_ids))
  
  for(cluster_view in as.character(cluster_ids)){
    for(s in sample_ids){
      spatnet.subset <- subset(spatnet, cluster_from %in% cluster_view)
      spatnet.subset <- subset(spatnet.subset, sample == s)
      spatnet.subset_conly <- subset(spatnet.subset, cluster_to %in% cluster_view)

      N <- length(unique(c(spatnet.subset$from)))
      ki <- spatnet.subset_conly %>%
        dplyr::group_by(from) %>%
        dplyr::count(name = "ki")
      L <- sum(ki$ki)/2
      k_avg <- (2*L)/N
      k_avg_list[[paste0("S", s)]] <- k_avg
    }
    k_avg_df_add <- as.data.frame(unlist(k_avg_list))
    colnames(k_avg_df_add) <- paste0("kavg_cluster_", cluster_view)
    k_avg_df <- cbind(k_avg_df, k_avg_df_add)
  }
  
  message("Done!")
  return(k_avg_df)
}


#' Run for multiple random seeds
se_spat_net <- GetSpatNet(se)
n_perm <- 50
avgk_df_perm_list <- list()
for(i in seq(1,n_perm)){
  avgk_df_perm_list[[i]] <- RandomClusteringSpatNet(se, column.cluster.id = "seurat_clusters", column.sample.id = "sample_id", 
                                                    se.SpatNet = se_spat_net,
                                                    random.seed = i
                                                    )
}

avgk_df_perm_avg <- data.frame(row.names = row.names(avgk_df_perm_list[[1]]))
avgk_df_perm_sd <- data.frame(row.names = row.names(avgk_df_perm_list[[1]]))
avgk_df_perm_me <- data.frame(row.names = row.names(avgk_df_perm_list[[1]]))
n_cols <- dim(avgk_df_perm_list[[1]])[2]
for(c in seq(1:n_cols)){
  cluster_avgks <- do.call(rbind, lapply(avgk_df_perm_list, `[[`, c))
  avgk_df_perm_avg$c_avg <- colMeans(cluster_avgks)
  avgk_df_perm_sd$c_sd <- colSds(cluster_avgks)
  avgk_df_perm_me$c_me <- qnorm(.95) * ( avgk_df_perm_avg$c_avg / sqrt(length(avgk_df_perm_list)) )  # margin of error at 90% CI
  colnames(avgk_df_perm_avg)[colnames(avgk_df_perm_avg)=="c_avg"] <- colnames(avgk_df_perm_sd)[colnames(avgk_df_perm_sd)=="c_sd"] <- colnames(avgk_df_perm_me)[colnames(avgk_df_perm_me)=="c_me"] <- colnames(avgk_df_perm_list[[1]])[c]
}


#' Compare random averages with real results
rownames(k_avg_out) <- k_avg_out$sample
k_avg_perm_diff <- k_avg_out[rownames(avgk_df_perm_avg), colnames(avgk_df_perm_avg)] - avgk_df_perm_avg
k_avg_perm_zscore <- (k_avg_out[rownames(avgk_df_perm_avg), colnames(avgk_df_perm_avg)] - avgk_df_perm_avg) / avgk_df_perm_sd



write.csv(avgk_df_perm_avg, file = file.path(DIR_OUT, "nbs_cluster_kavg-perm_avg.baseline.csv"), row.names = T)
write.csv(avgk_df_perm_sd, file = file.path(DIR_OUT, "nbs_cluster_kavg-perm_sd.baseline.csv"), row.names = T)
write.csv(avgk_df_perm_me, file = file.path(DIR_OUT, "nbs_cluster_kavg-perm_me.baseline.csv"), row.names = T)
write.csv(k_avg_perm_diff, file = file.path(DIR_OUT, "nbs_cluster_kavg-perm_diff_score.baseline.csv"), row.names = T)
write.csv(k_avg_perm_zscore, file = file.path(DIR_OUT, "nbs_cluster_kavg-perm_zscore.baseline.csv"), row.names = T)


#' Plot summary stats
#' kavg diff column
summary_df_kavg_diff <- as.data.frame(t(do.call(cbind, lapply(k_avg_perm_diff, summary))))
summary_df_kavg_diff$cluster <- paste0("C", as.character(sub(pattern = "kavg_cluster_", "", rownames(summary_df_kavg_diff))))

p1 <- ggplot(summary_df_kavg_diff, aes(x=reorder(cluster, Max.), y=Max.)) +
  geom_col(fill=color_low) +
  labs(x="") +
  coord_flip() +
  theme_classic()

p2 <- ggplot(summary_df_kavg_diff, aes(x=reorder(cluster, Mean), y=Mean)) +
  geom_col(fill=color_low) +
  labs(x="") +
  coord_flip() +
  theme_classic()

p <- p1 + p2 + plot_annotation(title = '<k> obs-exp difference') & theme(plot.title = element_text(hjust=0.5));p
pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_diff_stats.baseline.pdf")), width = 4, height = 3.5);p;dev.off()
png(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_diff_stats.baseline.png")), width = 4*300, height = 3.5*300, res=300);p;dev.off()


#' kavg diff box plot
df_kavg_diff <- as.data.frame(t(k_avg_perm_diff))
df_kavg_diff$cluster <- paste0("C", as.character(sub(pattern = "kavg_cluster_", "", rownames(df_kavg_diff))))
df_kavg_diff_long <- reshape2::melt(df_kavg_diff, id.vars = c("cluster"))

p3 <- ggplot(df_kavg_diff_long, aes(x=reorder(cluster, value), y=value)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_point(color=color_low2, size=.5) +
  labs(x="", y="<k>", title="<k> obs-exp diff.") +
  coord_flip() +
  theme_classic();p3

pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_diff_boxplot.baseline.pdf")), width = 3, height = 3.5);p3;dev.off()
png(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_diff_boxplot.baseline.png")), width = 3*300, height = 3.5*300, res=300);p3;dev.off()



#' kavg zscore column plot
summary_df_kavg_zscore <- as.data.frame(t(do.call(cbind, lapply(k_avg_perm_zscore, summary))))
summary_df_kavg_zscore$cluster <- paste0("C", as.character(sub(pattern = "kavg_cluster_", "", rownames(summary_df_kavg_zscore))))

p1 <- ggplot(summary_df_kavg_zscore, aes(x=reorder(cluster, Max.), y=Max.)) +
  geom_col(fill=color_low) +
  labs(x="") +
  coord_flip() +
  theme_classic();p1

p2 <- ggplot(summary_df_kavg_zscore, aes(x=reorder(cluster, Mean), y=Mean)) +
  geom_col(fill=color_low) +
  labs(x="") +
  coord_flip() +
  theme_classic()

p <- p1+p2 + plot_annotation(title = '<k> obs-exp z-score') & theme(plot.title = element_text(hjust=0.5));p
pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_zscore_stats.baseline.pdf")), width = 4, height = 3.5);p;dev.off()
png(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_zscore_stats.baseline.png")), width = 4*300, height = 3.5*300, res=300);p;dev.off()



#' kavg zscore box plot
df_kavg_zscore <- as.data.frame(t(k_avg_perm_zscore))
df_kavg_zscore$cluster <- paste0("C", as.character(sub(pattern = "kavg_cluster_", "", rownames(df_kavg_zscore))))
df_kavg_zscore_long <- reshape2::melt(df_kavg_zscore, id.vars = c("cluster"))

p1 <- ggplot(df_kavg_zscore_long, aes(x=reorder(cluster, value), y=value)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_point(color=color_low2, size=.5) +
  labs(x="", y="z-score", title="<k> obs-exp diff.") +
  coord_flip() +
  theme_classic();p1


pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_zscore_boxplot.baseline.pdf")), width = 3, height = 3.5);p1;dev.off()
png(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_zscore_boxplot.baseline.png")), width = 3*300, height = 3.5*300, res=300);p1;dev.off()


#####
sessionInfo()







