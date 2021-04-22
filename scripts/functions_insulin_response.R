#' Functions for looking at insulin response in each cluster and sample/subject
#' L.F, Feb 2021


#' Compute mean expession across spots in each subject and cluster
#'
#' @param se.object Seurat object
#' @param se.data.slot Slot to fetch data from. Default is "data" from active assay (SCT).
#' @param genes.include Vector of gene names to fetch data for. Default is NULL (use all genes in se.object)
#' @param se.assay.data (Optional) Provide Assay Data for seurat object. Analysis is designed for using nomalized data (i.e. slot="data"). Default is NULL.
#' @param clusters.include Vector of cluster IDs to be included from "seurat_clusters" column of se.object metadata
#' @param clusters.include Vector of cluster IDs to be included from "seurat_clusters" column of se.object metadata
#' @param sample.names Vector of sample IDs to be included from "sample_name" column of se.object metadata
#' @param sample.metadata Metadata corresponding to each sample included in analysis (set as row names)
#' @param spotn.cutoff Number of spots at least required in sample and cluster to be included (else set to 0). Default 10.
#'
#' @return List with 1) Output data, and 2) Output metadata
ComputeMeanSpotExpr <- function(
  se.object,
  se.assay.data = NULL, # se_ins_allc
  se.data.slot = "data",
  genes.include = NULL,
  sample.metadata,
  clusters.include,
  sample.names,
  spotn.cutoff = 10
) {
  
  # Check genes.include
  if (is.null(genes.include)) {
    genes.include <- rownames(se.object)
  }
  
  # Fetch assay data
  if (is.null(se.assay.data)) {
    message("\nFetching Assay Data from seurat object...")
    se_assay_data <- GetAssayData(se.object[genes.include, ], slot = se.data.slot)
  } else if ( nrow(se.assay.data) != length(genes.include) ) {
    warning("\nProvided se_assay_data does not match length of genes provided in genes.include")
  } else {
    se_assay_data <- se.assay.data
  }
  
  # Join cluster and sample names
  sc_names <- c()
  for(c in clusters.include){
    sc_names <- c(sc_names, paste(c, sample.names, sep=".") )
  }
  
  # Create frame for output metadata and output data
  output_metadata <- data.frame(row.names = sc_names,
                                cluster = rep("", length(sc_names)),
                                stringsAsFactors = F)
  
  output_data <- data.frame(row.names = rownames(se_assay_data))
  
  for (s in sample.names) {
    message(paste0("\nComputing spot mean values per cluster for sample ", s, "...\n             "))
    
    for (c in clusters.include) {
      mess <- sprintf("Cluster:%s", c)
      message(paste(rep("\b", nchar(mess)), collapse = ""), appendLF = FALSE)
      message(mess, appendLF = c == 9)
      
      # Select spots for sample and cluster
      s_spots <- colnames(subset(se.object, sample_name == s))
      c_spots <- colnames(subset(se.object, seurat_clusters == c))
      sc_spots <- intersect(c_spots, s_spots)
      
      # Calculate row means if number of spots >= spotn.cutoff
      
      if (length(sc_spots) < spotn.cutoff) {
        message(paste0("\nNumber of spots in cluster ", c, " in sample ", s, " is below cut off, setting to NA.\n             "))
        sc_data <- data.frame(value = rep(NA, nrow(se_assay_data)), row.names = rownames(se_assay_data))
        
      } else if (length(sc_spots) == 1) {
        sc_data <- data.frame(se_assay_data[, sc_spots])
        
      }  else if (length(sc_spots) == 0) {
        sc_data <- data.frame(value = rep(0, nrow(se_assay_data)), row.names = rownames(se_assay_data))
        
      } else if (length(sc_spots) >= spotn.cutoff) {
        sc_data <- data.frame(Matrix::rowMeans(se_assay_data[, sc_spots]))
        
      } 
      
      sc_name <- paste(c, s, sep=".")
      colnames(sc_data) <- sc_name
      
      output_data <- cbind(output_data, sc_data)
      output_metadata[sc_name, colnames(sample.metadata)] <- sample.metadata[s, colnames(sample.metadata)]  # c("insulin_stim", "bmi", "condition", "subject_id", "m_value")
      output_metadata[sc_name, "cluster"] <- paste0("C", c)
    }
  }
  output_metadata <- output_metadata[colnames(output_data), ]
  output_list <- list(output_data, output_metadata)
  
  message(paste0("\nReturning list with (1) Output data and (2) Output metadata"))
  return(output_list)
}


#' Compute ratio between before and after insulin stimulation in each cluster and sample
#'
#' @param bulk.data Data output from ComputeMeanSpotExpr().
#' @param bulk.metadata Metadata output from ComputeMeanSpotExpr().
#' @param clusters.include Vector of cluster IDs to be included.
#' @param subject.id.include Vector of sample IDs to be included.
#' @param column.subject.id Column name of column in bulk.metadata corresponding to subject.id.include. Default "subject_id".
#' @param column.insulin Column name of column in bulk.metadata corresponding to insulin stimulation. Must be 0 for no insulin and 1 for insulin. Default "insulin_stim".
#' @param pseudo.count Pseudo-count to be added when calculating ratio. Default 0.01.
#'
#' @return DataFrame with gene ratios for all genes in bulk.data, for each cluster and subject.
CalculateGeneInsulinRatio <- function (
  bulk.data,
  bulk.metadata,
  clusters.include,
  subject.id.include,
  column.subject.id = "subject_id",
  column.insulin = "insulin_stim",
  pseudo.count = 0.01
) {
  pc <- pseudo.count
  output_data_ratios <- data.frame(row.names = rownames(bulk.data),
                                   "gene" = rownames(bulk.data))
  for(s in subject.id.include) {
    message(paste0("\nComputing gene ratio per cluster for subject ", s, "...             "))
    
    for(c in clusters.include){
      c_id <- paste0("C",c)
      metadata_sc <- subset(bulk.metadata, cluster == c_id)
      metadata_sc <- metadata_sc[metadata_sc[, column.subject.id]==s, ]
      
      sample_c_0 <- rownames(metadata_sc[metadata_sc[, column.insulin]==0,])
      sample_c_1 <- rownames(metadata_sc[metadata_sc[, column.insulin]==1,])
      
      #' Calculate ratios
      if (length(sample_c_0) > 1 & length(sample_c_1) > 1) {
        #' If muliple samples per subject: Take gene average for samples and compute ratio
        data_ratios_sc <- data.frame(
          "ratio" = ( rowMeans(bulk.data[, sample_c_1] + pc) / rowMeans(bulk.data[, sample_c_0] + pc) ),
          row.names = rownames(bulk.data), 
          stringsAsFactors = F
        )
        
      } else if (!length(sample_c_0) > 1 & length(sample_c_1) > 1) {
        data_ratios_sc <- data.frame(
          "ratio" = ( rowMeans(bulk.data[, sample_c_1] + pc) / (bulk.data[, sample_c_0] + pc) ),
          row.names = rownames(bulk.data), 
          stringsAsFactors = F
        )
        
      } else if (length(sample_c_0) > 1 & !length(sample_c_1) > 1) {
        data_ratios_sc <- data.frame(
          "ratio" = ( (bulk.data[, sample_c_1] + pc) / rowMeans(bulk.data[, sample_c_0] + pc) ),
          row.names = rownames(bulk.data), 
          stringsAsFactors = F
        )
        
      } else if (length(sample_c_0) == 1 & length(sample_c_1) == 1) {
        data_ratios_sc <- data.frame(
          "ratio" = ( (bulk.data[, sample_c_1] + pc) / (bulk.data[, sample_c_0] + pc) ),
          row.names = rownames(bulk.data), 
          stringsAsFactors = F
        )
        
      } else {
        stop(paste0("No samples detected for subject ", s, " in the provided data"))
      }
      
      #' Prepare output
      colnames(data_ratios_sc) <- paste0(c, ".", s)
      data_ratios_sc$gene <- rownames(data_ratios_sc)
      output_data_ratios <- cbind(output_data_ratios, data_ratios_sc)
      
      if (ncol(output_data_ratios) > 1) {
        output_data_ratios$gene <- NULL
      }
      
    }
    
  }
  output_data_ratios$gene <- NULL
  return(output_data_ratios)
}



#' Calculate mean ratio for up-regulated and down-regulated genes separately
#'
#' @param ratio.data Data output from CalculateGeneInsulinRatio().
#' @param genes.up Vector of upregulated genes to compute score for.
#' @param genes.down Vector of downregulated genes to compute score for.
#' @param gene.cutoff.up Minimum number of genes from provided genes.up vector needed to be present to caluclate mean ratio score. Default 2.
#' @param gene.cutoff.down Minimum number of genes from provided genes.down vector needed to be present to caluclate mean ratio score. Default 2.
#' @param return.log2 Set to TRUE if output ratios should be returned as log2 values of the ratio means. Default FALSE. 
#'
#' @return DataFrame with mean up & down gene ratios for each cluster and subject.
CalculateGeneInsulinRatioMeans <- function (
  ratio.data,
  genes.up,
  genes.down,
  gene.cutoff.up = 2,
  gene.cutoff.down = 2,
  return.log2 = FALSE
) {
  output_mean_ratio <- data.frame(up = rep(1, ncol(ratio.data)),
                                  down = rep(1, ncol(ratio.data)),
                                  cluster = rep(NA, ncol(ratio.data)),
                                  subject = rep(NA, ncol(ratio.data)),
                                  stringsAsFactors = F)
  rownames(output_mean_ratio) <- colnames(ratio.data)
  
  for ( i in seq_along( colnames(ratio.data) ) ) {
    sc <- colnames(ratio.data)[i]
    c <- strsplit(sc, "\\.")[[1]][1]
    s <- strsplit(sc, "\\.")[[1]][2]
    
    output_mean_ratio[sc, "cluster"] <- paste0("C", c)
    output_mean_ratio[sc, "subject"] <- s
    
    ngenes_up <- sum(ratio.data[genes.up, sc] > 1)
    ngenes_down <- sum(ratio.data[genes.down, sc] > 1)
    
    if (is.na(ngenes_up)) {
      output_mean_ratio[sc, "up"] <- NA
      
    } else if (ngenes_up >= gene.cutoff.up) {
      output_mean_ratio[sc, "up"] <- mean(ratio.data[genes.up, sc], na.rm=T)
      
    }
    
    if (is.na(ngenes_down)) {
      output_mean_ratio[sc, "down"] <- NA
      
    } else if (ngenes_down >= gene.cutoff.down) {
      output_mean_ratio[sc, "down"] <- mean(ratio.data[genes.down, sc], na.rm=T)
      
    }

  }
  
  output_mean_log2ratio <- output_mean_ratio
  
  if (return.log2 == TRUE) {
    output_mean_log2ratio[, c("up", "down")] <- log2(output_mean_log2ratio[, c("up", "down")])
    return(output_mean_log2ratio)
    
  } else {
    return(output_mean_log2ratio)
  }
}

