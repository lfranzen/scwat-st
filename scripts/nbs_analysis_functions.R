#'====================================================
#' Functions for performing the Visium neighborhood analyses 
#' to compute homotypic and heterotypic clustering scores.
#' 
#' L. Franzén, lovisa.franzen@scilifelab.se
#' June 2021
#'====================================================

library(STutility)
library(magrittr)


#'====================================================================================
# HOMOTYPIC SCORE CALCULATION - WITHIN-CLUSTER ANALYSIS
#'====================================================================================

#' Calculate average degree for spatial data
#' 
#' @param se.object Seurat (STUtility) object containing cluster and sample identities for each spot in the metadata.
#' @param column.cluster.id Column name in metadata corresponding to cluster ID of the spots. Default "seurat_clusters".
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
#' @param column.cluster.id Column name in metadata corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param column.sample.id Column name in metadata corresponding to Sample ID of the spots.
#' @param n.perm Number of random permutation to run (positive integer).
#' @param random.seed (Optional) Random seed for the random sampling. Default is NA.
#' @return Data frame, or list of data frames, containing average degree for network.
#' @export
RandomClustersAvgDegree <- function (
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
      k_avg_perm_list[[i]] <- RandomClustersAvgDegree(se.object = se.object, 
                                                      column.sample.id = column.sample.id, 
                                                      column.cluster.id = column.cluster.id, 
                                                      se.SpatNet = se_spat_net,
                                                      n.perm=1)
      }
    return(k_avg_perm_list)
  } else if (n.perm < 0) {stop("Number of random permutations must be a positive integer")}
}


#' Compare observed vs expected average degree scores by computing the "homotypic score" (avg k exp-obs)
#' 
#' @param avk.k.observed data frame with the observed average degree values for each cluster and sample (output from CalculateAvgDegree function)
#' @param avk.k.expected List of data frames with the randomised average degree values for each cluster and sample (output from RandomClustersAvgDegree function where n.perm>1)
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
# HETEROTYPIC SCORE CALCULATION - BETWEEN-CLUSTER ANALYSIS
#'====================================================================================

#' Create Adjacency Matrix
#' 
#' @param nbs.df Data frame with output from STUtility::GetSpatNet() data for all clusters of interest. Rows should correspond to spot ID and columns should include cluster IDs as well as colums starting with "nbs_".
#' @param column.cluster.id Column name in nbs.df corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param cluster.include Vector of cluster IDs to include in your analysis. Default is all (NA).
#' @return Adjacency matrix with number of neighbours present between each cluster pair
#' @export
CreateAdjMatrix <- function(
  nbs.df,
  column.cluster.id = "seurat_clusters",
  cluster.include = NA
){
  if (length(cluster.include)==1 & is.na(cluster.include)) {
    cluster.include <- seq(1, max(as.numeric(nbs.df[,1])), 1)
  }
  
  nbs_adjmat <- matrix(0L, nrow = length(cluster.include), ncol = length(cluster.include))
  
  for (i in seq_along(cluster.include)) {
    c <- cluster.include[i]
    c_nbs_df <- nbs.df[nbs.df[,column.cluster.id]==c, ]
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
#' @param column.cluster.id Column name in nbs.df corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param cluster.include Vector of cluster IDs to include in your analysis.
#' @param column.subject.id Column name in se.metadata corresponding to Sample ID of the spots.
#' @return List of adjacency matrices with number of neighbours present between each cluster pair. Each matrix in the list corresponds to data from one sample.
#' @export
CreateAdjMatrixPerSubject <- function(
  nbs.df,
  se.metadata,
  column.cluster.id = "seurat_clusters",
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
                                  column.cluster.id = column.cluster.id)
    
    nbs_adjmat_list[[subject_id]] <- nbs_adjmat
  }
  
  return(nbs_adjmat_list)
}


#' Randomise Cluster IDs within Subject data
#' 
#' @param se.object Seurat (STUtility) object containing cluster and sample identities for each spot in the metadata.
#' @param column.cluster.id Column name in metadata corresponding to cluster ID of the spots. Default "seurat_clusters".
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


#' Randomise Cluster IDs within Subject data and calculate average degree
#' 
#' @param se.object Seurat (STUtility) object containing cluster and sample identities for each spot in the metadata.
#' @param column.cluster.id Column name in metadata corresponding to cluster ID of the spots. Default "seurat_clusters".
#' @param column.sample.id Column name in metadata corresponding to Sample ID of the spots.
#' @param cluster.include Vector of cluster IDs to include in your analysis. Default is all (NA). 
#' @param n.perm Number of random permutation to run (positive integer).
#' @param random.seed (Optional) TRUE/FALSE. Random seed for the random sampling. Default is FALSE
#' @return List of data frames containing adjacency matrices.
#' @export
RandomClustersAdjMat <- function (
  se.object,
  column.cluster.id,
  column.sample.id,
  cluster.include = NA,
  n.perm = 1,
  random.seed = FALSE
) {
  
  perm_adj_mat_list <- list()
  
  for (i in 1:n.perm) {
    message(paste0("Running function for multiple permutations: ", i, "/", n.perm))
    if (random.seed) {
      se.object <- RandomiseClusteringIDs(se.object = se.object, 
                                          column.cluster.id = column.cluster.id, 
                                          column.sample.id = column.sample.id, 
                                          random.seed = i)
    } else {
      se.object <- RandomiseClusteringIDs(se.object = se.object, 
                                          column.cluster.id = column.cluster.id, 
                                          column.sample.id = column.sample.id)
    }
    
    se.object <- SetIdent(se.object, value = "clusters_perm")
    
    for (column_rm in grep(pattern = "nbs_", colnames(se.object@meta.data), value = T)) {
      se.object[[column_rm]] <- NULL
    }
    
    if (length(cluster.include)==1 & is.na(cluster.include)){
      message("Returning values for all clusters")
      cluster.include <- seq(1, max(as.numeric(se.object@meta.data[, column.cluster.id])), 1)
    }
    
    for (c in cluster.include){
      se.object <- STutility::RegionNeighbours(se.object, 
                                               id = c, 
                                               keep.within.id = T, 
                                               verbose = TRUE)
    }
    
    perm_nbs_df <- se.object@meta.data[, c("clusters_perm", grep(pattern = "nbs_", colnames(se.object@meta.data), value = T))]
    
    perm_adj_mat <- CreateAdjMatrix(nbs.df = perm_nbs_df, 
                                    cluster.include = cluster.include, 
                                    column.cluster.id = "clusters_perm")
    perm_adj_mat_list[[i]] <- perm_adj_mat
  }
  return(perm_adj_mat_list)
}


#' Calculate position-wise mean and standard deviation from list of adjacency matrices
#' 
#' @param adj.mat.list List of ajacency matrices to calculate position-wise mean and standard deviation
#' @return S4 object with mean and standard deviation (sd) matrices
#' @export
CalulateMeanSdAdjMatList <- function (
  adj.mat.list
  ) {
  n_cols <- dim(adj.mat.list[[1]])[1]
  perm_adj_mat_avg <- perm_adj_mat_sd <- matrix(0L, nrow = n_cols, ncol = n_cols)
  for(i in 1:n_cols){
    for(j in 1:n_cols){
      list_ij <- c()
      for(list_n in 1:length(adj.mat.list)){
        list_ij <- c(list_ij, adj.mat.list[[list_n]][i,j])
      }
      perm_adj_mat_avg[i,j] <- mean(list_ij)
      perm_adj_mat_sd[i,j] <- sd(list_ij)
    }
  }
  setClass("AdjMat", representation(mean = "matrix", sd = "matrix"))
  output_obj <- new("AdjMat", 
                    mean=perm_adj_mat_avg, 
                    sd=perm_adj_mat_sd)
  return(output_obj)
}


#' Compare observed vs expected average degree scores by computing the "heterotypic score"
#' 
#' @param adj.mat.observed Adjacency matrix with observed values
#' @param adj.mat.expected.mean Adjacency matrix with position-wise mean of expected values from multiple permutations
#' @param adj.mat.expected.sd Adjacency matrix with position-wise standard deviation (sd) of expected values from multiple permutations
#' @param cluster.ids (Optional) Provide vector with names of clusters
#' @return Adjacency matrix with position-wise heterotypic scores (z-scores)
#' @export
CalculateHeterotypicScore <- function (
  adj.mat.observed,
  adj.mat.expected.mean,
  adj.mat.expected.sd,
  cluster.ids = NA
) {
  nbs_adjmat_zsore <- round( ((adj.mat.observed - adj.mat.expected.mean) / adj.mat.expected.sd) , digits = 3)
  diag(nbs_adjmat_zsore) <- diag(adj.mat.observed)
  if(!is.na(cluster.ids)){
    rownames(nbs_adjmat_zsore) <- colnames(nbs_adjmat_zsore) <- cluster.ids
  } else (
    rownames(nbs_adjmat_zsore) <- colnames(nbs_adjmat_zsore) <- paste0("C", as.character(1:ncol(adj.mat.observed)))
  )
  return(nbs_adjmat_zsore)
}

