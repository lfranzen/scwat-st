#' Functions for performing the Visium neighborhood analyses 
#' to compute homotypic and heterotypic clustering scores.

library(STutility)
library(magrittr)


#'====================================================================================
# HOMOTYPIC SCORE CALCULATION - WITHIN-CLUSTER ANALYSIS
#'====================================================================================

#' Calculate average degree for spatial data
#' 
#' @param se.object Seurat (STUtility) object containing cluster and sample identities for each spot in the metadata.
#' @param column.clusters.id Column name in metadata corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param column.sample.id Column name in metadata corresponding to Sample ID of the spots.
#' @param se.SpatNet (Optional) Provide the output from GetSpatNet(se.object). Default is NA (i.e computed within the function).
#' @return Data frame containing average degree for network.
#' @export
CalculateAvgDegree <- function (
  se.object,
  column.cluster.id,
  column.sample.id,
  se.SpatNet = NA
) {
  n_samples <- length(unique(se.object@meta.data[, column.sample.id]))
  n_clusters <- length(unique(se.object@meta.data[, column.cluster.id]))
  
  #' Create SpatNet and make network
  if(is.list(se.SpatNet)){
    spatnet_init <- se.SpatNet
  } else {
    message("Creating network...")
    spatnet_init <- GetSpatNet(se.object)
  }
  
  spatnet <- do.call(rbind, lapply(seq_along(spatnet_init), function(i) {
    spnet <- spatnet_init[[i]]
    spnet$cluster_from <- se.object[[]][spnet$from, column.cluster.id]
    spnet$cluster_to <- se.object[[]][spnet$to, column.cluster.id]
    spnet$sample <- paste0(i)
    return(spnet)
  }))
  spatnet$cluster_from <- as.factor(as.numeric(spatnet$cluster_from))
  spatnet$cluster_to <- as.factor(as.numeric(spatnet$cluster_to))
  
  
  #' Calculate avgerage degree for all samplse and clusters
  message("Calculating average degree for each cluster and sample...")
  
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
  
  message("Done!")
  return(k_avg_df)
}



#' Randomise Cluster IDs within Subject data and calculate average degree
#' 
#' @param se.object Seurat (STUtility) object containing cluster and sample identities for each spot in the metadata.
#' @param column.clusters.id Column name in metadata corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param column.sample.id Column name in metadata corresponding to Sample ID of the spots.
#' @param random.seed (Optional) Random seed for the random sampling. Default is NA.
#' @return Data frame containing average degree for network.
#' @export
RandomClusteringAvgDegree <- function (
  se.object,
  column.cluster.id,
  column.sample.id,
  random.seed = NA,
  se.SpatNet = NA,
  n.perm = 1
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
    dplyr::mutate(clusters_perm = clusters_original[sample(dplyr::row_number())])
  
  #' Add shuffled clusters to se object metadata
  se.object <- AddMetaData(se.object, as.character(se_metadata_perm$clusters_perm), col.name = "clusters_perm")
  
  #' Calculate average degree for randomised cluster network. Run for all permutations if n.perm > 1
  if (n.perm == 1) {
    k_avg_df <- CalculateAvgDegree(se.object = se.object, 
                                   column.sample.id = column.sample.id,
                                   column.cluster.id = "clusters_perm"
    )
    return(k_avg_df)
  } else if (n.perm > 1) {
    k_avg_perm_list <- list()
    for(i in seq(1, n.perm)){
      message(paste0("Running for multiple permutations: ", i, "/", n.perm))
      se_spat_net <- GetSpatNet(se.object)
      k_avg_perm_list[[i]] <- RandomClusteringAvgDegree(se.object = se.object, 
                                                        column.sample.id = column.sample.id, 
                                                        column.cluster.id = column.cluster.id, 
                                                        se.SpatNet = se_spat_net,
                                                        n.perm=1)
      }
    return(k_avg_perm_list)
  } else if (n.perm < 0) {warning("Number of random permutations must be a positive integer")}
}


#' Compare observed vs expected average degree scores by computing the "homotypic score" (avg k exp-obs)
#' 
#' @param avk.k.observed data frame with the observed average degree values for each cluster and sample (output from CalculateAvgDegree function)
#' @param avk.k.expected List of data frames with the randomised average degree values for each cluster and sample (output from RandomClusteringAvgDegree function where n.perm>1)
#' @param return.score Which format the homotypic score output should be provided as. Defult is "difference". Other option is "z-score".
#' @return Data frame with homotypic scores for each cluster 
#' @export
CalculateHomotypicScore <- function (
  avk.k.observed,
  avk.k.expected,
  return.score = "difference"
) {
  if (!is.list(avk.k.expected)) {
    stop("The provided avk.k.expected must be a list")
  }
  
  #' Calculate average scores and sd across the random iterations ("expected")
  avgk_df_perm_avg <- data.frame(row.names = row.names(avk.k.expected[[1]]))
  avgk_df_perm_sd <- data.frame(row.names = row.names(avk.k.expected[[1]]))
  n_cols <- dim(avk.k.expected[[1]])[2]-1
  
  for(c in seq(1:n_cols)){
    cluster_avgks <- do.call(rbind, lapply(avk.k.expected, `[[`, c))
    avgk_df_perm_avg$c_avg <- colMeans(cluster_avgks)
    avgk_df_perm_sd$c_sd <- matrixStats::colSds(cluster_avgks)
    colnames(avgk_df_perm_avg)[colnames(avgk_df_perm_avg)=="c_avg"] <- 
      colnames(avgk_df_perm_sd)[colnames(avgk_df_perm_sd)=="c_sd"] <- 
      colnames(avk.k.expected[[1]])[c]
  }
  
  #' Calculate homotypic score by comparing observed vs expected
  rownames(avk.k.observed) <- avk.k.observed$sample
  
  if (return.score == "difference") {
    k_avg_perm_diff <- avk.k.observed[rownames(avgk_df_perm_avg), colnames(avgk_df_perm_avg)] - avgk_df_perm_avg
    return(k_avg_perm_diff)
  } else if (return.score == "z.score") {
    k_avg_perm_zscore <- (avk.k.observed[rownames(avgk_df_perm_avg), colnames(avgk_df_perm_avg)] - avgk_df_perm_avg) / avgk_df_perm_sd
    return(k_avg_perm_zscore)
  }
}


#'====================================================================================
# HETEROTYPIC SCORE CALCULATION - BETWEEN-CLUSTER ANALYSIS - TO DO...
#'====================================================================================






