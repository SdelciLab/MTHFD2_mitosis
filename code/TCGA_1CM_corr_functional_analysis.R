# libraries
library(tidyverse)
library(gprofiler2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(KEGGgraph)
library(ggpubr)

# 1. Get input data -------------------------------------------------------------------------------------------------

# get correlated genes
poscor_files <- list.files(path="../1CM_cor_genes/Cor_genes/", pattern = "poscor.*10", full.names = T)
read_files <- map(poscor_files, read_tsv)

# get dimensions
map(read_files, dim)

# obtain background
background <- read_tsv("../functional_analysis_cor_genes/background_allgenes.tsv")

# correlated genes and background
pos_core_cor_genes <- map(read_files, ~ .[["Ensembl_Gene_ID"]])
background_genes <- background$Ensembl_Gene_ID

# remove version
ensID_rm_version <- function(x) {
  for (i in seq_along(x)) {
    x[i] <- str_split(x[i], "\\.")[[1]][1]
  }
  x
}
pos_core_cor_genes_nov <- map(pos_core_cor_genes, ensID_rm_version)

# 2. Perform GO BP enrichment with Cluster Profiler ----------------------------------------------------------------------------

# perform BP for all of them
obtainBP <- function(pos_genes) {
    enrichGO(
        gene = pos_genes,
        universe = background_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.01,
        qvalueCutoff  = 0.05,
        readable      = TRUE)
}

BP_results <- map(pos_core_cor_genes_nov, obtainBP)
BP_results_simp <- map(BP_results, simplify)

# no results in 3, 12, 15, 17, 18

# filter significant results
filter_BP <- function(x) {
    df <- x@result %>%
        filter(qvalue < 0.05 & pvalue < 0.01)
    df
}
BP_results_filter <- map(BP_results, filter_BP)
BP_results_filter_simp <- map(BP_results_simp, filter_BP)

# export table of BP for each gene
gene_names <- list("AHCY", "ATIC", "CHDH", "DHFR", "DNMT1", "DNMT3A", "DNMT3B", 
                   "GART", "MAT2A", "MAT2B", "MTHFD1", "MTHFD1L", "MTHFD2", "MTHFD2L",
                   "MTHFR", "MTR", "PEMT", "SHMT1", "SHMT2", "TYMS")

export_BP <- function(df, name) {
    write.csv2(df, paste0("../1CM_cor_genes/BP_filtered/", name, "_BP_filtered.csv"))
}
walk2(BP_results_filter, gene_names, export_BP)

export_BP_simp <- function(df, name) {
    write.csv2(df, paste0("../1CM_cor_genes/BP_filtered_simp/", name, "_BP_filtered_simp.csv"))
}
walk2(BP_results_filter_simp, gene_names, export_BP_simp)

# get dimensions
map(BP_results_filter, nrow)
map(BP_results_filter_simp, nrow)

# plot results from MTHFD2

# get first 20 BP terms
MTHFD2_BP_results_filter_simp <- BP_results_filter_simp[[13]]
MTHFD2_BP_results_filter_simp <- MTHFD2_BP_results_filter_simp[1:20, ]

# remove redundant terms
# remove regulation of chromosome segregation: exact same genes as regulation of chromosome separation
# add function and gene ratio
MTHFD2_BP_results_filter_simp_circus <- MTHFD2_BP_results_filter_simp %>%
  filter(!(Description %in% c("regulation of chromosome segregation"))) %>%
  mutate(cluster = c("mitosis", "mitosis", "mitosis", "mitosis", "mitosis", "mitosis", "cell cycle", "mitosis", "mitosis", "other",
                     "cell cycle", "mitosis", "mitosis", "cell cycle", "other", "mitosis", "mitosis", "mitosis", "other")) %>%
  mutate(GeneRatio_real = Count/75) 
MTHFD2_BP_results_filter_simp_circus$cluster <- factor(MTHFD2_BP_results_filter_simp_circus$cluster, 
                                                       levels = c("mitosis", "cell cycle", "other"))

# filter by gene ratio and reorder
MTHFD2_BP_results_filter_simp_circus <- MTHFD2_BP_results_filter_simp_circus %>%
  filter(GeneRatio_real >= 0.10) %>%
  arrange(-GeneRatio_real) %>%
  arrange(cluster) %>%
  mutate(order = seq(1:14))

# circus plot
ggplot(MTHFD2_BP_results_filter_simp_circus,
       aes(x=fct_reorder(Description, order), y=GeneRatio_real, fill=cluster, size = `-log10(qvalue)`)) +
  geom_point(stroke = NA, shape = 21) +
  coord_polar() +
  scale_fill_manual(values = c("#559A73", "#9EFFC8", "darkgrey")) +
  scale_size(range = c(2, 6)) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(size  = 0.25, linetype = 1, color = "grey"),
        axis.text = element_text(size = 9, color = "black"), axis.title = element_text(size = 12))
ggsave("../functional_analysis_cor_genes/clusterProfiler/GOBP_terms_circus_point.pdf", device = "pdf", 
       width = 3.7, height = 3.7)

# plot results from all the other folate enzymes 

# adapt results from each enzyme
# add function and gene ratio
ATIC_results <- BP_results_filter_simp[[2]]
 
ATIC_results_add <- ATIC_results %>%
    mutate(Enzyme = rep("ATIC", nrow(ATIC_results)),
           cluster = ifelse(ID %in% c("GO:0033044"), "mitosis", 
                            ifelse(ID %in% c(""), "cell cycle", "other")),
           GeneRatio_real = Count/88,
           `-log10(qvalue)` = -log10(qvalue))

DHFR_results <- BP_results_filter_simp[[4]]

DHFR_results_add <- DHFR_results %>%
    mutate(Enzyme = rep("DHFR", nrow(DHFR_results)),
           cluster = ifelse(ID %in% c("GO:0007059", "GO:0000280", "GO:0098813", "GO:0140014",
                                      "GO:0000819", "GO:0033044", "GO:0007051", "GO:0034502",
                                      "GO:0051983", "GO:0007088", "GO:0051303", "GO:0050000",
                                      "GO:0051783", "GO:0000910", "GO:0051653", "GO:0040001"), "mitosis", 
                            ifelse(ID %in% c("GO:0006260", "GO:0000075", "GO:1901987", "GO:0051321",
                                             "GO:0006275", "GO:0044839", "GO:0044843", "GO:0051302"), 
                                   "cell cycle", "other")),
           GeneRatio_real = Count/1353,
           `-log10(qvalue)` = -log10(qvalue))

GART_results <- BP_results_filter_simp[[8]]

GART_results_add <- GART_results %>%
    mutate(Enzyme = rep("GART", nrow(GART_results)),
           cluster = ifelse(ID %in% c("GO:0007059", "GO:0033044", "GO:0098813", "GO:0000819",
                                      "GO:0007051", "GO:0044772", "GO:0051983", "GO:0050000",
                                      "GO:0045839", "GO:0000910", "GO:0051653", "GO:0032465",
                                      "GO:0040001"), "mitosis", 
                            ifelse(ID %in% c("GO:0006260", "GO:0000075", "GO:1901987", "GO:0006275",
                                             "GO:0044839", "GO:0051321", "GO:0044843"), 
                                   "cell cycle", "other")),
           GeneRatio_real = Count/2751,
           `-log10(qvalue)` = -log10(qvalue))

MTHFD1_results <- BP_results_filter_simp[[11]]

MTHFD1_results_add <- MTHFD1_results %>%
    mutate(Enzyme = rep("MTHFD1", nrow(MTHFD1_results)),
           cluster = ifelse(ID %in% c("GO:0007059", "GO:0098813", "GO:0000819", "GO:0044772",
                                      "GO:0033044", "GO:0007051", "GO:0051983", "GO:0034502",
                                      "GO:0007088", "GO:0051783", "GO:0051303", "GO:0050000",
                                      "GO:0000910", "GO:0034508", "GO:0032465", "GO:0090224",
                                      "GO:0051653"), "mitosis", 
                            ifelse(ID %in% c("GO:0006260", "GO:1901987", "GO:0000075", "GO:0051321",
                                             "GO:0044839", "GO:0006275", "GO:0044843", "GO:0051302"), 
                                   "cell cycle", "other")),
           GeneRatio_real = Count/843,
           `-log10(qvalue)` = -log10(qvalue))

MTR_results <- BP_results_filter_simp[[16]]

MTR_results_add <- MTR_results %>%
    mutate(Enzyme = rep("MTR", nrow(MTR_results)),
           cluster = ifelse(ID %in% c("GO:0033044", "GO:0007059", "GO:0007062", "GO:0044772",
                                      "GO:0098813", "GO:0000910", "GO:0007346", "GO:0032465"), "mitosis", 
                            ifelse(ID %in% c("GO:0006260", "GO:0044843", "GO:0045787", "GO:0006275"), 
                                   "cell cycle", "other")),
           GeneRatio_real = Count/2821,
           `-log10(qvalue)` = -log10(qvalue))

SHMT2_results <- BP_results_filter_simp[[19]]

SHMT2_results_add <- SHMT2_results %>%
    mutate(Enzyme = rep("SHMT2", nrow(SHMT2_results)),
           cluster = ifelse(ID %in% c("GO:0033044", "GO:0000819", "GO:0051988", "GO:0098813",
                                      "GO:0007052", "GO:0045931", "GO:0140014"), "mitosis", 
                            ifelse(ID %in% c("GO:0000086", "GO:0044839", "GO:1901987"), 
                                   "cell cycle", "other")),
           GeneRatio_real = Count/24,
           `-log10(qvalue)` = -log10(qvalue))


TYMS_results <- BP_results_filter_simp[[20]]

TYMS_results_add <- TYMS_results %>%
    mutate(Enzyme = rep("TYMS", nrow(TYMS_results)),
           cluster = ifelse(ID %in% c("GO:0007059", "GO:0000280", "GO:0140014", "GO:0098813",
                                      "GO:0007346", "GO:0007051", "GO:0051783", "GO:0007088", 
                                      "GO:0051983", "GO:0033044", "GO:0034508", "GO:0050000", 
                                      "GO:0034502", "GO:0000910", "GO:0051653", "GO:0040001", 
                                      "GO:0090224","GO:0040020", "GO:0031145"), "mitosis", 
                            ifelse(ID %in% c("GO:0006260", "GO:1901987", "GO:0051321", "GO:0044839",
                                             "GO:0090329", "GO:0044843", "GO:0051302", "GO:0051445",
                                             "GO:1901993"), 
                                   "cell cycle", "other")),
           GeneRatio_real = Count/297,
           `-log10(qvalue)` = -log10(qvalue))

# put data all enzymes together
all_enzymes <- ATIC_results_add[1:10,] %>%
    rbind(DHFR_results_add[1:10,]) %>%
    rbind(GART_results_add[1:10,]) %>%
    rbind(MTHFD1_results_add[1:10,]) %>%
    rbind(MTR_results_add[1:10,]) %>%
    rbind(SHMT2_results_add[1:10,]) %>%
    rbind(TYMS_results_add[1:10,])

all_enzymes$cluster <- factor(all_enzymes$cluster, levels = c("mitosis", "cell cycle", "other"))

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
    if (!is.list(within)) {
        within <- list(within)
    }
    
    new_x <- do.call(paste, c(list(x, sep = sep), within))
    stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., labels = reorder_func) {
    
    ggplot2::scale_x_discrete(labels = labels, ...)
}
reorder_func <- function(x, sep = "___") {
    reg <- paste0(sep, ".+$")
    gsub(reg, "", x)
}

# final plot
ggplot(all_enzymes, aes(reorder_within(Description, -`-log10(qvalue)`, Enzyme), `-log10(qvalue)`, 
                        col = cluster, size = GeneRatio_real)) +
    geom_point() +
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "lightgrey", linewidth = 0.3) +
    labs(x = "", size = "Gene Ratio") +
    facet_wrap(~Enzyme, scales = "free", ncol = 2) +
    scale_x_reordered() +
    scale_color_manual(values = c("#559A73", "#9EFFC8", "darkgrey")) +
    coord_flip() +
    theme_bw() +
    theme(legend.title = element_text(size = 11), legend.text = element_text(size = 8)) +
    theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 8)) +
    theme(title = element_text(size = 11), plot.title = element_text(hjust = 0.5),
          panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.03))
ggsave("../1CM_cor_genes/clusterProfiler/all_enzymes_topbp.pdf", device = "pdf", height = 6.5, width = 9.5) 

# 3. Get colors for mitosis and cell cycle enrichment  ----------------------------------------------------------------------------

# get colors for 1CM schematic
mitosis_cell_cycle_summary <- read.csv2("../1CM_cor_genes/1CM_corgenes_mitosis_cellcycle.csv")

# Scale your values to range between 0 and 1
mitosis_rr <- range(mitosis_cell_cycle_summary$Mitotic.terms....)
mitosis_sc <- (mitosis_cell_cycle_summary$Mitotic.terms.... -mitosis_rr[1])/diff(mitosis_rr)

mitosis_color_range <- colorRamp(c("#FFFFFF", "#1F653D"))
mitosis_color <- rgb(mitosis_color_range(mitosis_sc)/255)
mitosis_color

# Check that it works
image(seq_along(mitosis_sc), 1, as.matrix(seq_along(mitosis_sc)), col=mitosis_color,
      axes=FALSE, xlab="", ylab="")

# Scale your values to range between 0 and 1
cellcycle_rr <- range(mitosis_cell_cycle_summary$Cell.cycle....)
cellcycle_sc <- (mitosis_cell_cycle_summary$Cell.cycle.... -cellcycle_rr[1])/diff(cellcycle_rr)

cellcycle_color_range <- colorRamp(c("#FFFFFF", "#239E58"))
cellcycle_color <- rgb(cellcycle_color_range(cellcycle_sc)/255)
cellcycle_color

# Check that it works
image(seq_along(cellcycle_sc), 1, as.matrix(seq_along(cellcycle_sc)), col=cellcycle_color,
      axes=FALSE, xlab="", ylab="")
