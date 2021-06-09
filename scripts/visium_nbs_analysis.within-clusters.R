########################
#' Within-cluster analysis: "Homotypic score algorithm"
#' 
#' Visium neighbourhood analysis of cluster community formation or aggregation  
#' by computing the average degree for each cluster and sample
#' 
#' L. Franzén, lovisa.franzen@scilifelab.se
#' Jan 2021
########################

#' SET UP

#' Select analysis:

ANALYSIS <- "baseline"
# ANALYSIS <- "insulin"


#' Load libs
library(magrittr)
library(STutility)
library(patchwork)

#' Define paths
DIR_IN <- file.path(getwd(), "data", "visium")
DIR_OUT <- file.path(getwd(), "results", "visium", "tables")
DIR_WD <- file.path(getwd(), "scripts")

source(file.path(DIR_WD, "colors.R"))

# ===================================
#' READ DATA & LOAD IMGS

metadata <- read.delim(file.path(DIR_IN, "visium_sample_metadata.tsv"), sep = "\t", stringsAsFactors = F)

if (ANALYSIS == "baseline") {
  se_base <- readRDS(file.path(DIR_OUT, "..", "se-object.visium_baseline.rds"))
  se <- se_base
} else if (ANALYSIS == "insulin") {
  se_ins <- readRDS(file.path(DIR_OUT, "..", "se-object.visium_insulin.rds"))
  se <- se_ins
  se <- AddMetaData(se, metadata = paste0(se$subject_id, "_", se$insulin_stim, "_", se$sample_id), col.name = "subject_id")
}

se <- LoadImages(se, time.resolve = T, verbose = T, xdim = 100)

#' Summary data
se_stats <- se@meta.data %>%
  dplyr::group_by(subject_id, tissue_id, seu_n, novaseq_id) %>% 
  dplyr::count(name = "spots")
se_stats$seu_n <- paste0("S", se_stats$seu_n)
colnames(se_stats)[colnames(se_stats)=="seu_n"] <- "sample"

n_samples <- length(se_stats$sample)
n_clusters <- length(unique(se$seurat_clusters))


# ===================================
#' Define functions

#' Randomise Cluster IDs within Subject data and calculate average degree
#' 
#' @param se.object Seurat (STUtility) object containing cluster and sample identities for each spot in the metadata.
#' @param column.clusters.id Column name in metadata corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param column.sample.id Column name in metadata corresponding to Sample ID of the spots.
#' @param random.seed (Optional) Random seed for the random sampling. Default is NA.
#' @param se.SpatNet (Optional) Provide the output from GetSpatNet(se.object). Default is NA (i.e computed within the function).
#' @return Data frame containing average degree for network.
#' @export
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


# ===================================
#' RUN NBS ANALYSIS: "OBSERVED"
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

pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_k.spatial2.", ANALYSIS, ".pdf")), width = 12, height = 9, useDingbats = F)
for(cluster_view in as.character(sort(unique(se$seurat_clusters))) ){  #c("1", "4", "5", "6", "7", "12", "14", "17")
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
}
dev.off()


#####
#' Count avgerage degree (k_avg)
# N <- length(unique(c(spatnet.subset$from)))
# ki <- spatnet.subset_conly %>%
#   dplyr::group_by(from) %>%
#   dplyr::count(name = "ki")
# L <- sum(ki$ki)/2
# k_avg <- (2*L)/N


#' Calculate avg degree for all samples and clusters
k_avg_list <- list()
k_avg_df <- data.frame(row.names = paste0("S", seq(1, n_samples)))

for(cluster_view in as.character(seq(1,n_clusters))){
  for(s in seq(1, n_samples)){
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
if (ANALYSIS == "baseline") {
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
} else if ( ANALYSIS == "insulin") {
  sample_conv <- t(data.frame(
                              S1 = "S41",
                              S2 = "S42",
                              S3 = "S43",
                              S4 = "S44",
                              S5 = "S45",
                              S6 = "S46",
                              S7 = "S47",
                              S8 = "S52",
                              S9 = "S53", 
                              stringsAsFactors = F, row.names = "novaseq_id"))
}
sample_conv <- as.data.frame(sample_conv)
sample_conv$sample <- rownames(sample_conv)
se_stats_neworder <- se_stats
se_stats_neworder$sample <- NULL
se_stats_neworder <- merge(se_stats_neworder, sample_conv, by="novaseq_id")
k_avg_out <- merge(se_stats_neworder[,c(1:3,5)], k_avg_df, by = "sample")

#' Export table
write.csv(k_avg_out, file = file.path(DIR_OUT, paste0("nbs_cluster_kavg.", ANALYSIS, ".csv")), row.names = F)


# ===================================
#' EXPECTED: RANDOMISED DATA
#' 
#' Run for multiple random seeds
se_spat_net <- GetSpatNet(se)
n_perm <- 50
avgk_df_perm_list <- list()
for(i in seq(1, n_perm)){
  avgk_df_perm_list[[i]] <- RandomClusteringSpatNet(se, 
                                                    column.cluster.id = "seurat_clusters", 
                                                    column.sample.id = "sample_id", 
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
  colnames(avgk_df_perm_avg)[colnames(avgk_df_perm_avg)=="c_avg"] <- 
    colnames(avgk_df_perm_sd)[colnames(avgk_df_perm_sd)=="c_sd"] <- 
    colnames(avgk_df_perm_me)[colnames(avgk_df_perm_me)=="c_me"] <- 
    colnames(avgk_df_perm_list[[1]])[c]
}


#' Compare random averages with real results
rownames(k_avg_out) <- k_avg_out$sample
k_avg_perm_diff <- k_avg_out[rownames(avgk_df_perm_avg), colnames(avgk_df_perm_avg)] - avgk_df_perm_avg
k_avg_perm_zscore <- (k_avg_out[rownames(avgk_df_perm_avg), colnames(avgk_df_perm_avg)] - avgk_df_perm_avg) / avgk_df_perm_sd

#' Export tables
write.csv(avgk_df_perm_avg, file = file.path(DIR_OUT, paste0("nbs_cluster_kavg-perm_avg.", ANALYSIS, ".csv")), row.names = T)
write.csv(avgk_df_perm_sd, file = file.path(DIR_OUT, paste0("nbs_cluster_kavg-perm_sd.", ANALYSIS, ".csv")), row.names = T)
write.csv(avgk_df_perm_me, file = file.path(DIR_OUT, paste0("nbs_cluster_kavg-perm_me.", ANALYSIS, ".csv")), row.names = T)
write.csv(k_avg_perm_diff, file = file.path(DIR_OUT, paste0("nbs_cluster_kavg-perm_diff_score.", ANALYSIS, ".csv")), row.names = T)
write.csv(k_avg_perm_zscore, file = file.path(DIR_OUT, paste0("nbs_cluster_kavg-perm_zscore.", ANALYSIS, ".csv")), row.names = T)

# k_avg_perm_diff <- read.csv(file = file.path(DIR_OUT, paste0("nbs_cluster_kavg-perm_diff_score.", ANALYSIS, ".csv")), row.names = 1)
# k_avg_perm_zscore <- read.csv(file = file.path(DIR_OUT, paste0("nbs_cluster_kavg-perm_zscore.", ANALYSIS, ".csv")), row.names = 1)



# ===================================
#' PLOT SUMMARY STATS
#' 
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
pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_diff_stats.", ANALYSIS, ".pdf")), width = 4, height = 3.5, useDingbats = F);p;dev.off()
# png(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_diff_stats.", ANALYSIS, ".png")), width = 4*300, height = 3.5*300, res=300);p;dev.off()


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

p4 <- ggplot(df_kavg_diff_long, aes(x=reorder(cluster, value), y=value)) +
  geom_col(data = summary_df_kavg_diff, aes(x=reorder(cluster, Mean), y=Mean), fill="white", color = "black", width = .8) +
  geom_hline(yintercept = 0) +
  # geom_boxplot() +
  geom_point(color=color_low2, size=1) +
  labs(x="", y="<k>", title="<k> obs-exp diff.") +
  coord_flip() +
  theme_classic();p4

pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_diff_boxplot.", ANALYSIS, ".pdf")), width = 2.5, height = 3.5, useDingbats = F);p3;dev.off()
pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_diff_barplot.", ANALYSIS, ".pdf")), width = 2.5, height = 3.5, useDingbats = F);p4;dev.off()


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
pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_zscore_stats.", ANALYSIS, ".pdf")), width = 4, height = 3.5, useDingbats = F);p;dev.off()
# png(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_zscore_stats.", ANALYSIS, ".png")), width = 4*300, height = 3.5*300, res=300);p;dev.off()



#' Plot kavg zscore box plot
df_kavg_zscore <- as.data.frame(t(k_avg_perm_zscore))
df_kavg_zscore$cluster <- paste0("C", as.character(sub(pattern = "kavg_cluster_", "", rownames(df_kavg_zscore))))
df_kavg_zscore_long <- reshape2::melt(df_kavg_zscore, id.vars = c("cluster"))

p1 <- ggplot(df_kavg_zscore_long, aes(x=reorder(cluster, value), y=value)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  geom_point(color=color_low2, size=.5) +
  labs(x="", y="z-score", title="<k> obs-exp diff.") +
  coord_flip() +
  theme_classic()
pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_kavg-perm_zscore_boxplot.", ANALYSIS, ".pdf")), width = 3, height = 3.5, useDingbats = F);p1;dev.off()

  


# ===================================
#' EXTRA PLOTS
#' Plot networks for manuscript figures:

# OBSERVED

subjects_include <- c("S44", "S46", "S49", "S50", "S51", "S55")  # sample_conv[sample_conv$novaseq_id %in% c("S44", "S46", "S49", "S50", "S51", "S55"), "sample"])
clusters_include <- c(6, 8, 9, 10, 11, 14, 18, 19, 20)

for(sample_novaseq in subjects_include) {  # 6
  pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_k.spatial_observed_example_", sample_novaseq, "." , ANALYSIS, ".pdf")), width = 8, height = 7.5, useDingbats = F)
  
  for(cluster_view in clusters_include){  # 4
    sample_view <- gsub("S", "", sample_conv[sample_conv$novaseq_id == sample_novaseq, "sample"])
    
    spatnet.subset <- subset(spatnet, cluster_from %in% cluster_view & sample %in% sample_view)
    spatnet.subset_conly <- subset(spatnet.subset, cluster_to %in% cluster_view & sample %in% sample_view)
    
    N <- length(unique(c(spatnet.subset$from)))
    ki <- spatnet.subset_conly %>%
      dplyr::group_by(from) %>%
      dplyr::count(name = "ki")
    
    spatnet.subset_ki <- spatnet.subset
    spatnet.subset_ki <- merge(spatnet.subset_ki, as.data.frame(ki), by="from", all=T)
    spatnet.subset_ki[is.na(spatnet.subset_ki$ki), "ki"] <- 0
    
    p_ob <- ggplot() +
      geom_segment(data = spatnet.subset_conly, aes(x = start_x, xend = end_x, y = -start_y, yend = -end_y), color = "white", size=0.3) +
      geom_point(data = spatnet.subset_ki, aes(start_x, -start_y, color = as.factor(ki)), size = 1) +
      labs(color="", title = paste("sample", sample_novaseq, " - cluster", cluster_view)) +
      # scale_color_manual(values = colors_clusters[colors_clusters$seurat_clusters==cluster_view, "cluster_color"]) +
      scale_color_manual(values = viridis::viridis(7)) +
      xlim(c(0,2000)) +
      ylim(c(0,-2000)) +
      theme_void() +
      theme(panel.background = element_rect(fill = "black"))
      # NoLegend()
    print(p_ob)
    
    write.table(spatnet.subset_ki, file.path(DIR_OUT, "../tables/nbs_kavg_plot_tables", paste0("nbs_cluster_k.spatial_observed_example_", sample_novaseq, "-C", cluster_view, "_spatnetSubsetKi." , ANALYSIS, ".tsv")),
                sep = "\t", row.names = T, col.names = T, quote = F)
    write.table(spatnet.subset_conly, file.path(DIR_OUT, "../tables/nbs_kavg_plot_tables", paste0("nbs_cluster_k.spatial_observed_example_", sample_novaseq, "-C", cluster_view, "_spatnetSubsetCOnly." , ANALYSIS, ".tsv")),
                sep = "\t", row.names = T, col.names = T, quote = F)

    # pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_k.spatial_observed_example_S", sample_novaseq, "-C", cluster_view, "." , ANALYSIS, ".pdf")), width = 4, height = 3.5, useDingbats = F);p_ob;dev.off()
  }
  dev.off()
}

# pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_k.spatial_observed_example_S", sample_view, "-C", cluster_view, "." , ANALYSIS, ".pdf")), width = 4, height = 3.5, useDingbats = F);p_ob;dev.off()


ggplot() +
  geom_segment(data = spatnet.subset_conly, aes(x = start_x, xend = end_x, y = -start_y, yend = -end_y), color = "white", size=0.3) +
  geom_point(data = spatnet.subset_ki, aes(start_x, -start_y, color = as.factor(ki)), size = 1.5) +
  labs(color="", title = paste("sample", sample_novaseq, " - cluster", cluster_view)) +
  scale_color_manual(values = viridis::viridis(7)) +
  # scale_x_continuous(limits = c(0,2000)) +
  # scale_y_continuous(limits = c(0,-2000)) +
  xlim(c(0,2000)) +
  ylim(c(0,-2000)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "black"))



# ===================================
# PERMUTED
se.object <- se
se_metadata <- se.object@meta.data[, c("seurat_clusters", "sample_id")]
se_metadata$clusters_original <- se_metadata[, "seurat_clusters"]
se_metadata$sample_id <- se_metadata[, "sample_id"]

for(rand_seed in 1:3){
  set.seed(seed = rand_seed)
  se_metadata_perm <- se_metadata %>% 
    dplyr::group_by(sample_id) %>% 
    dplyr::mutate(clusters_perm = clusters_original[sample(dplyr::row_number())])
  se.object <- AddMetaData(se.object, as.character(se_metadata_perm$clusters_perm), col.name = "clusters_perm")
  
  spatnet_perm <- do.call(rbind, lapply(seq_along(spatnet_init), function(i) {
    spatnet_perm <- spatnet_init[[i]]
    spatnet_perm$cluster_from <- se.object[[]][spatnet_perm$from, "clusters_perm"]
    spatnet_perm$cluster_to <- se.object[[]][spatnet_perm$to, "clusters_perm"]
    spatnet_perm$sample <- paste0(i)
    return(spatnet_perm)
  }))
  spatnet_perm$cluster_from <- as.factor(as.numeric(spatnet_perm$cluster_from))
  spatnet_perm$cluster_to <- as.factor(as.numeric(spatnet_perm$cluster_to))
  
  for(cluster_view in c("6")){
    for(sample_view in "4"){
      spatnet.subset <- subset(spatnet_perm, cluster_from %in% cluster_view & sample %in% sample_view)
      spatnet.subset_conly <- subset(spatnet.subset, cluster_to %in% cluster_view & sample %in% sample_view)
      
      p_perm <- ggplot() +
        geom_point(data = spatnet.subset, aes(start_x, -start_y, color = cluster_from), size = .5) +
        geom_segment(data = spatnet.subset_conly, aes(x = start_x, xend = end_x, y = -start_y, yend = -end_y), size=0.3) +
        labs(color="", title = paste("sample", sample_view, " - cluster", cluster_view)) +
        scale_color_manual(values = colors_clusters[colors_clusters$seurat_clusters==cluster_view,"cluster_color"]) +
        theme_void() +
        NoLegend()
      print(p_perm)
    }
  }
  pdf(file.path(DIR_OUT, "../figures", paste0("nbs_cluster_k.spatial_permuted_example_seed", rand_seed, ".", ANALYSIS, ".pdf")), width = 3, height = 3.5, useDingbats = F);print(p_perm);dev.off()
}


#===============
sessionInfo()







