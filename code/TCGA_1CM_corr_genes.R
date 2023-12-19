# libraries
library(tidyverse)

# setwd
setwd("/users/ssdelci/npardo/TCGA_RNAseq_analysis")

# 1. Read input ---------------------------------------------------------------------

# obtain sample sheet of tumor data
sample_sheet <- read_tsv("sample_sheet_tumor.tsv") %>%
    column_to_rownames("CompleteName")
projects <- sort(unique(sample_sheet$ProjectID))
print("Head of sample sheet of tumor data")
head(sample_sheet)
print("Projects")
print(projects)

# obtain normalised and filtered expression data 
expr_data <- read_tsv("all_expression_data_TPM_filt_log2.tsv")
print("Dimension of normalised and filtered expression data")
dim(expr_data) 

# 2. Obtain sample sheet and expression data per project--------------------------------------------------

expr_data_m <- expr_data %>%
    column_to_rownames(var = "Ensembl_Gene_ID") %>%
    as.matrix()

# subset expression data matrix to keep only tumor data
tumor_data <- rownames(sample_sheet)
expr_data_m_tumor <- expr_data_m[, tumor_data]
print("Dimension of expression data of tumours")
dim(expr_data_m_tumor)

sample_sheet_per_project <- split(sample_sheet, sample_sheet$ProjectID)
names_per_project <- map(sample_sheet_per_project, rownames)

print("Structure of sample sheet per project")
str(sample_sheet_per_project)

expr_data_per_project <- vector("list", length(names_per_project))
for (i in seq_along(names_per_project)) {
    expr_data_per_project[[i]] <- expr_data_m_tumor[, names_per_project[[i]]]
}

print("Structure of tumor expression data per project")
str(expr_data_per_project)

# eliminate genes lowly expressed in each tumor type, so that log(TPM+0.01) < 0 in all samples of each project
expr_data_per_project_filt <- vector("list", length(expr_data_per_project))
for (i in seq_along(expr_data_per_project)) {
    expr_data_per_project_filt[[i]] <- expr_data_per_project[[i]][apply(expr_data_per_project[[i]], 1, function(x) !all(x < 0)), ]
}
print("Structure of filtered tumor expression data per project")
str(expr_data_per_project_filt)

# 3. Correlations between 1CM genes and other genes  ------------------------------------------------------------------------

# correlation of expression of gene of interest (GOI) compared to all genes with cor.test

# NOTE: change script for each 1CM gene

# MTHFD2 ENSG00000065911.10
# MTHFD2L ENSG00000163738.17
# MTHFD1 ENSG00000100714.14
# MTHFD1L ENSG00000120254.14
# MTR ENSG00000116984.11
# MTHFR ENSG00000177000.9
# TYMS ENSG00000176890.14
# ATIC ENSG00000138363.13
# GART ENSG00000159131.15
# SHMT1 ENSG00000176974.16
# SHMT2 ENSG00000182199.9
# DHFR ENSG00000228716.5

# transpose so that genes in cols and patients in rows
# we want cor between genes and cor is made between cols
get_cor_genes <- function(m, gene) {
    cor_test_results <- apply(t(m), 2, cor.test, t(m)[ , "ENSG00000065911.10"], method = "pearson")
}

cor_genes_per_project <- map(expr_data_per_project_filt, get_cor_genes)
print("Structure of matrix before and after transpose")
dim(expr_data_per_project_filt[[1]])
dim(t(expr_data_per_project_filt[[1]]))

# obtain estimates
estimates_per_project <- map_depth(cor_genes_per_project, 2, "estimate")
estimates_per_project <- map(estimates_per_project, unlist)
print("Estimates")
str(estimates_per_project)

# obtain pvalues
pvalues_per_project <- map_depth(cor_genes_per_project, 2, "p.value")
pvalues_per_project <- map(pvalues_per_project, unlist)
print("pvalues")
str(pvalues_per_project)

# obtain adjusted pvalues
adj_pvalues_per_project <- map(pvalues_per_project, p.adjust, "BH")
print("pvalues adjusted")
str(adj_pvalues_per_project)

# transform results to df
cor_genes_df_per_project <- vector("list", length(estimates_per_project))
for (i in seq_along(estimates_per_project)) {
    cor_genes_df_per_project[[i]] <- tibble(Ensembl_Gene_ID = names(pvalues_per_project[[i]]), estimates = estimates_per_project[[i]], 
                                            pvalues = pvalues_per_project[[i]], pvalues_adjusted = adj_pvalues_per_project[[i]])
}

print("Structure of df of correlated genes")
str(cor_genes_df_per_project)

# filter those with positive corr >=0.6 and padj lower than 0.05 and eliminate GOI (cor = 1)
get_poscor_genes <- function(df, n) {
    df %>%
        filter(estimates >= n & pvalues_adjusted < 0.05) %>%
        filter(Ensembl_Gene_ID != "ENSG00000065911.10")
}

poscor_genes_per_project <- map(cor_genes_df_per_project, get_poscor_genes, 0.6)
print("Structure of positive correlated genes")
str(poscor_genes_per_project)

# filter those with negative corr <=-0.4 and padj lower than 0.05
get_negcor_genes <- function(df, n) {
    df %>%
        filter(estimates <= n & pvalues_adjusted < 0.05)
}

negcor_genes_per_project <- map(cor_genes_df_per_project, get_negcor_genes, -0.4)
print("Structure of negative correlated genes")
str(negcor_genes_per_project)

# 4. Find core correlated genes -----------------------------------------------------------

# obtain all ensembl ids
all_pos_cor_genes <- map(poscor_genes_per_project, "Ensembl_Gene_ID") %>%
    unlist() %>%
    unique() %>%
    unname()
length(all_pos_cor_genes) 
all_pos_cor_genes <- data.frame(Ensembl_Gene_ID = all_pos_cor_genes)
str(all_pos_cor_genes)

all_neg_cor_genes <- map(negcor_genes_per_project, "Ensembl_Gene_ID") %>%
    unlist() %>%
    unique() %>%
    unname()
length(all_neg_cor_genes) 
all_neg_cor_genes <- data.frame(Ensembl_Gene_ID = all_neg_cor_genes)
str(all_neg_cor_genes)

# create a column for each tumor type indicating if these genes were correlated with GOI in this tumor type
for (i in seq_along(poscor_genes_per_project)) {
    genes <- poscor_genes_per_project[[i]] %>% pull(Ensembl_Gene_ID)
    new_col <- ifelse(all_pos_cor_genes$Ensembl_Gene_ID %in% genes, TRUE, FALSE)
    all_pos_cor_genes <- cbind(all_pos_cor_genes, project = new_col)
}
colnames(all_pos_cor_genes) <- c("Ensembl_Gene_ID", projects)
str(all_pos_cor_genes)

for (i in seq_along(negcor_genes_per_project)) {
    genes <- negcor_genes_per_project[[i]] %>% pull(Ensembl_Gene_ID)
    new_col <- ifelse(all_neg_cor_genes$Ensembl_Gene_ID %in% genes, TRUE, FALSE)
    all_neg_cor_genes <- cbind(all_neg_cor_genes, project = new_col)
}
names(all_neg_cor_genes) <- c("Ensembl_Gene_ID", projects)
str(all_neg_cor_genes)

# obtain the total: number of tumor types in which the gene was found correlated with MTHFD2

all_pos_cor_genes <- all_pos_cor_genes %>%
    dplyr::mutate(total = rowSums(.[2:32])) %>%
    arrange(desc(total))
nrow(all_pos_cor_genes) 

all_neg_cor_genes <- all_neg_cor_genes %>%
    dplyr::mutate(total = rowSums(.[2:32])) %>%
    arrange(desc(total))
nrow(all_neg_cor_genes) 

# obtain those genes which are present in at least 10 tumor types
all_pos_cor_genes_10 <- all_pos_cor_genes %>%
    filter(total >= 10)
nrow(all_pos_cor_genes_10) 

all_neg_cor_genes_10 <- all_neg_cor_genes %>%
    filter(total >= 10)
nrow(all_neg_cor_genes_10) 

# write output
write_tsv(all_pos_cor_genes, "1CM_cor_genes/core_MTHFD2_poscor_genes_all.tsv")
write_tsv(all_neg_cor_genes, "1CM_cor_genes/core_MTHFD2_negcor_genes_all.tsv")
write_tsv(all_pos_cor_genes_10, "1CM_cor_genes/core_MTHFD2_poscor_genes_10.tsv")
write_tsv(all_neg_cor_genes_10, "1CM_cor_genes/core_MTHFD2_negcor_genes_10.tsv")