#' Analyse impact of bmi/age on "heterotypic clustering score"
#' 
#' LF, 2021 Apr

library(magrittr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

#' Define paths
DIR_ROOT <- getwd()
DIR_WD <- file.path(getwd(), "scripts")
DIR_RES <- file.path(DIR_ROOT, "results" , "visium")
DIR_FIG <- file.path(DIR_RES, "figures")
DIR_OBJ <- file.path(DIR_RES, "tables")


#' Read data and tables
metadata <- read.table(file.path(DIR_ROOT, "data", "visium", "visium_sample_metadata.tsv"), header = T)

file_pattern <- "nbs_adj_matrix_permscore_subject"  # ".baseline.tsv"
fnames <- grep(pattern = "baseline", x = list.files(path = DIR_OBJ, pattern = file_pattern, recursive = T), value = T)

adjmat_list <- list()
for(f in fnames){
  subject_id <- strsplit(strsplit(f, split = "\\.")[[1]][1], split = "-")[[1]][2]
  if(subject_id=="AK2019"){
    subject_id <- "AK2019-27"
  }
  mat <- read.table(file.path(DIR_OBJ, f), header = T)
  adjmat_list[[subject_id]] <- mat
}


#' Format new data
cluster_ids <- colnames(adjmat_list[[1]])

mat2 <- mat
mat2[lower.tri(mat2)] <- NA
diag(mat2) <- NA

n_combos <- sum(!is.na(mat2))
i_combos <- which(!is.na(mat2), arr.ind = TRUE)

dat <- as.data.frame(i_combos)
dat$c1 <- rownames(mat)[dat$row]
dat$c2 <- colnames(mat)[dat$col]
dat$c_c <- paste0(dat$c1, "_", dat$c2)
rownames(dat) <- dat$c_c

for(cc in dat$c_c){
  for(s in names(adjmat_list)){
    mat_s <- adjmat_list[[s]]
    val <- mat_s[dat[cc, "row"], dat[cc, "col"]]
    dat[cc, s] <- val
  }
}

dat_long <- reshape2::melt(dat, measure.vars = names(adjmat_list), variable.name = "subject")
dat_long[, c("row", "col")] <- NULL

#' Add metadata 
metadata_filt <- subset(metadata, insulin_stim==0 | is.na(insulin_stim))
mdat_add <- metadata_filt[, c("subject_id", "novaseq_id", "gender", "condition", "m_value", "bmi", "age")]

dat_long <- merge(dat_long, mdat_add, by.x="subject", by.y="subject_id")
dat_long$subject <- as.character(dat_long$subject)

#' export data
write.table(dat_long, file = file.path(DIR_OBJ, "nbs_adjmat_subjects", "visium_nbs_analysis.between-clusters.summarised_subjects.baseline.tsv"), sep = "\t", row.names = F)


#' Update with age for all
dat_long <- read.table(file.path(DIR_OBJ, "nbs_adjmat_subjects", "visium_nbs_analysis.between-clusters.summarised_subjects.baseline.tsv"), sep = "\t", header = T, stringsAsFactors = F)
df_subj_age_bmi <- read.csv(file.path(DIR_OBJ, "nbs_adjmat_subjects", "subject_age_bmi.csv"), header = T, stringsAsFactors = F)
df_subj_age_bmi$subject <- gsub(pattern = " ", replacement = "", df_subj_age_bmi$subject)

dat_long$age <- NULL 
dat_long <- merge(dat_long, df_subj_age_bmi[,c("subject", "age")], by="subject") 


#' PLOT 
#' BMI
p <- ggplot(dat_long, aes(x=bmi, y=value)) +
  geom_point(size=0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  # facet_grid(c1~c2) +
  facet_wrap(~c_c) +
  theme_linedraw() +
  theme(panel.grid = element_blank())

pdf(file = file.path(DIR_FIG, "visium_nbs_analysis.between-clusters.summarised_subjects_bmi.baseline.pdf"), width = 20, height = 20, useDingbats = F);p;dev.off()

#' AGE
p <- ggplot(dat_long, aes(x=age, y=value)) +
  geom_point(size=0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  # facet_grid(c1~c2) +
  facet_wrap(~c_c) +
  theme_linedraw() +
  theme(panel.grid = element_blank())

pdf(file = file.path(DIR_FIG, "visium_nbs_analysis.between-clusters.summarised_subjects_age.baseline.pdf"), width = 20, height = 20, useDingbats = F);p;dev.off()


#' STATS

#' vs BMI
ttest_res_list <- list()
lm_res_list <- list()
wilcox_res_list <- list()
for(cc in unique(dat_long$c_c)){
  test_dat <- subset(dat_long, c_c==cc)
  a <- subset(test_dat, bmi>=30)$value
  b <- subset(test_dat, bmi<30)$value
  res_ttest <- t.test(a, b, paired = F, alternative = "two.sided")
  res_lm <- lm(formula = value~bmi, data = test_dat)
  res_wilcox <- wilcox.test(a,b,alternative = "two.sided")
  
  ttest_res_list[[cc]] <- res_ttest
  lm_res_list[[cc]] <- res_lm
  wilcox_res_list[[cc]] <- res_wilcox
}

v1 <- sapply(ttest_res_list, function(x) x$p.value)
v2 <- sapply(lm_res_list, function(x) summary(x)$adj.r.squared)
v3 <- sapply(wilcox_res_list, function(x) x$p.value)

names(v1) <- names(v2) <- names(v3) <- names(ttest_res_list)
res_df <- as.data.frame(cbind(v1,v2,v3))
colnames(res_df) <- c("t.test.pval", "lm.adj.r2", "twogroup.mann.whitney.pval")


res_df$t.test.pval.fdr <- p.adjust(p = res_df$t.test.pval, n = length(res_df$t.test.pval), method = "fdr")
res_df$mannwhitney.pval.fdr <- p.adjust(p = res_df$twogroup.mann.whitney.pval, n = length(res_df$twogroup.mann.whitney.pval), method = "fdr")

head(dplyr::arrange(res_df, -abs(lm.adj.r2)), 10)
head(dplyr::arrange(res_df, t.test.pval), 10)
head(dplyr::arrange(res_df, t.test.pval.fdr), 10)
head(dplyr::arrange(res_df, twogroup.mann.whitney.pval), 10)
head(dplyr::arrange(res_df, mannwhitney.pval.fdr), 10)

res_df$cc <- rownames(res_df)
write.table(res_df, file = file.path(DIR_OBJ, "nbs_adjmat_subjects", "visium_nbs_analysis.between-clusters.stats_results_vs_bmi.baseline.tsv"), sep = "\t", row.names = F)


# vs AGE
age_top <- sort(df_subj_age_bmi$age)[1:5]
age_bottom <- sort(df_subj_age_bmi$age)[6:10]

ttest_res_list <- list()
lm_res_list <- list()
wilcox_res_list <- list()
for(cc in unique(dat_long$c_c)){
  test_dat <- subset(dat_long, c_c==cc)
  a <- subset(test_dat, age %in% age_bottom)$value
  b <- subset(test_dat, age %in% age_top)$value
  res_ttest <- t.test(a, b, paired = F, alternative = "two.sided")
  res_lm <- lm(formula = value~age, data = test_dat)
  res_wilcox <- wilcox.test(a,b,alternative = "two.sided")
  
  ttest_res_list[[cc]] <- res_ttest
  lm_res_list[[cc]] <- res_lm
  wilcox_res_list[[cc]] <- res_wilcox
}

v1 <- sapply(ttest_res_list, function(x) x$p.value)
v2 <- sapply(lm_res_list, function(x) summary(x)$adj.r.squared)
v3 <- sapply(wilcox_res_list, function(x) x$p.value)

names(v1) <- names(v2) <- names(v3) <- names(ttest_res_list)
res_df_age <- as.data.frame(cbind(v1,v2,v3))
colnames(res_df_age) <- c("t.test.pval", "lm.adj.r2", "twogroup.mann.whitney.pval")

res_df_age$t.test.pval.fdr <- p.adjust(p = res_df_age$t.test.pval, n = length(res_df_age$t.test.pval), method = "fdr")
res_df_age$mannwhitney.pval.fdr <- p.adjust(p = res_df_age$twogroup.mann.whitney.pval, n = length(res_df_age$twogroup.mann.whitney.pval), method = "fdr")

res_df_age$cc <- rownames(res_df_age)
write.table(res_df_age, file = file.path(DIR_OBJ, "nbs_adjmat_subjects", "visium_nbs_analysis.between-clusters.stats_results_vs_age.baseline.tsv"), sep = "\t", row.names = F)


