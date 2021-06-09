#'===========================
#' Pathway analysis of adipocytes
#' 
#' 
#' Decription:
#' Markers for adipocytes identified using Seurat::FindMarker() function
#' where the three adipocyte subtypes were compared against all other
#' annotated clusters (i.e not including the "unspecific" clusters).
#' Genes with a avg_logFC>0.1 and adjusted p-value <0.05 were selected and
#' converted to entrez IDs. Two genes lacked entrez IDs (FAM213A, PLA2G16) 
#' and were excluded. Pathway ennrichment analysis was thereafter performed
#' using the clusterProfiler R package to enrich for GO-terms using enrichGO(),
#' Reactome pathways using ReactomePA::enrichPathway(), and KEGG pathways 
#' using enrichKEGG() with default settings applied (BH correction, p-value < 0.05).
#' 
#' 
#' L. Franzén, 2021 June 
#'===========================

ANALYSIS_ID <- "visium_baseline"
ANALYSIS <- "baseline"

#' Load libs
library(magrittr)
library(patchwork)
library(STutility)
library(ReactomePA)
require(clusterProfiler)
library(enrichplot)
library(writexl)


#' Define paths
DIR_IN <- file.path(getwd(), "data", "visium")
DIR_OUT <- file.path(getwd(), "results", "visium", "tables")
DIR_WD <- file.path(getwd(), "scripts")

#' Read data
se <- readRDS(file.path(DIR_OUT, "..", "se-object.visium_baseline.rds"))

#============================================
#' Identify markers
# VlnPlot(se, features = c("PLIN4", "LEP"), pt.size = 0, group.by = "cluster_group")

se <- SetIdent(se, value = "cluster_group")
logf_thr <- 0.1
markers_adi <- FindMarkers(se, 
                           ident.1 = "Adipocyte", 
                           ident.2 = c("Immune", "Vascular", "Fibroblast/Pread"),  # "Unspecific"
                           only.pos = TRUE, 
                           min.pct = 0.1, 
                           logfc.threshold = logf_thr)

marker_genes <- rownames(markers_adi)


#' Pathway enrichment
d <- markers_adi
d <- subset(d, avg_logFC>0 & p_val_adj<0.05)
d$gene_name <- rownames(d)
eg = clusterProfiler::bitr(d$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
rownames(eg) <- eg$SYMBOL
colnames(eg) <- c("gene_name", "gene_entrez")

d <- merge(d, eg, by = "gene_name", all = T)s
rownames(d) <- d$gene_name

pw_genes <- d$gene_entrez[!is.na(d$gene_entrez)]
pw <- clusterProfiler::enrichGO(pw_genes, OrgDb="org.Hs.eg.db")
pw_ra <- enrichPathway(pw_genes)
pw_kegg <- enrichKEGG(pw_genes)


#' Export excel
d$gene_entrez <- as.character(d$gene_entrez)

fname <- paste0(ANALYSIS_ID, ".adipocyte_markers_pathwayenrichment")
sheets <- list()
sheets["Markers"] <- list(d)
sheets["GO_enrichemnt"] <- list(pw@result)
sheets["RA_enrichemnt"] <- list(pw_ra@result)
sheets["KEGG_enrichemnt"] <- list(pw_kegg@result)

write_xlsx(
  x = sheets,
  path = file.path(DIR_OUT, paste0(fname, ".xlsx")),
  col_names = TRUE,
  format_headers = TRUE
)

#' Plots
color_low <- "#441153"
color_low2 <- "#62376e"
p1 <- dotplot(pw, showCategory=10) + ggtitle("GO") + scale_color_gradient(low = color_low2, high = "grey80")
p2 <- dotplot(pw_ra, showCategory=10) + ggtitle("Reactome") + scale_color_gradient(low = color_low2, high = "grey80")
p3 <- dotplot(pw_kegg, showCategory=10) + ggtitle("KEGG") + scale_color_gradient(low = color_low2, high = "grey80")

p <- p1/p2/p3
fname <- paste0(ANALYSIS_ID, ".adipocyte_markers_pathwayenrichment")
pdf(useDingbats = F, file = file.path(DIR_OUT, "../figures", paste0(fname, ".pdf")), width = 11, height = 11);p;dev.off()

#'===========================