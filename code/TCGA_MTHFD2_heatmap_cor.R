# libraries
library(tidyverse)
library(car)
library(biomaRt)
library(data.table)
library(pheatmap)

# 1. Read input ---------------------------------------------------------------------

# files path
pos_cor_genes_path <- list.files(path = "../MTHFD2_cor_genes", pattern = "poscor_genes.tsv$", recursive = F, full.names = T)
names(pos_cor_genes_path) <- projects
neg_cor_genes_path <- list.files(path = "../MTHFD2_cor_genes", pattern = "negcor_genes.tsv$", recursive = F, full.names = T)   
names(neg_cor_genes_path) <- projects

# read all files
read <- function(path) {
    file <- read_tsv(path)
    file                 
    
}
pos_cor_genes_files <- map(pos_cor_genes_path, read)
neg_cor_genes_files <- map(neg_cor_genes_path, read)

str(pos_cor_genes_files) # list of 31 tibbles, n obs x 4 variables

# 2. Heatmap of MTHFD2 core positively correlated genes --------------------------------------------------------------------

# put all data together for a heatmap
pos_cor_genes_files_tumor <- list()
for (i in 1:length(pos_cor_genes_files)) {
    pos_cor_genes_files_tumor[[i]] <- pos_cor_genes_files[[i]] %>%
        mutate(type = rep(projects[[i]], nrow(pos_cor_genes_files[[i]])))
}
names(pos_cor_genes_files_tumor) <- names(pos_cor_genes_files)
df_all_pos_cor_tumor <- rbindlist(pos_cor_genes_files_tumor)
str(df_all_pos_cor_tumor) # 17678 obs 5 var

# change column types
df_all_pos_cor_tumor$estimates <- as.numeric(df_all_pos_cor_tumor$estimates)
df_all_pos_cor_tumor$pvalues <- as.numeric(df_all_pos_cor_tumor$pvalues)
df_all_pos_cor_tumor$pvalues_adjusted <- as.numeric(df_all_pos_cor_tumor$pvalues_adjusted)
str(df_all_pos_cor_tumor)

# remove TCGA-
df_all_pos_cor_tumor$type <- str_remove(df_all_pos_cor_tumor$type, "TCGA-")

# check if they are all significant
table(df_all_pos_cor_tumor$pvalues_adjusted < 0.05) # TRUE

# get matrix of data
matrix_all_pos_cor_tumor <- df_all_pos_cor_tumor %>%
    dplyr::select(Ensembl_Gene_ID, estimates, type) %>%
    arrange(desc(estimates)) %>%
    pivot_wider(names_from = type, values_from = estimates) %>%
    column_to_rownames("Ensembl_Gene_ID") %>%
    as.matrix()

# transpose - genes to columns and types to rows
matrix_all_pos_cor_tumor <- t(matrix_all_pos_cor_tumor)

# get those genes with are present in at least 10 different cancer types
df_all_pos_cor_tumor_top <- matrix_all_pos_cor_tumor[,colSums(!is.na(matrix_all_pos_cor_tumor))>=10]
dim(df_all_pos_cor_tumor_top) # dim 29 and 81

# transform NA values to 0 for the clustering
df_all_pos_cor_tumor_top[is.na(df_all_pos_cor_tumor_top)] <- 0

# clustering - get the order
set.seed(1234)
order <- pheatmap(df_all_pos_cor_tumor_top, show_colnames =  F, cluster_rows = T, cluster_cols = T)
order_genes <- colnames(df_all_pos_cor_tumor_top)[order$tree_col$order]
order_genes

# reorder the matrix, and change 0 to NAs
matrix_reorder <-  df_all_pos_cor_tumor_top[,order_genes]
matrix_reorder[matrix_reorder == 0] <- NA
pheatmap(matrix_reorder,  show_colnames =  F, cluster_rows = F, cluster_cols = F)

# plot heatmap
pdf(file= "../plots/heatmap_top_pos_cor_genes.pdf", width = 7, height = 2.4)
pheatmap(matrix_reorder, show_colnames =  F, cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(c("#9FD7E4", "#000066"))(1000),
         na_col = "whitesmoke", border_color = NA, fontsize = 6)
dev.off()
