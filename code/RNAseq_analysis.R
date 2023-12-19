# load libraries
library(tidyverse)
library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(biomaRt)
library(apeglm)
library(ggrepel)
library(VennDiagram)
library(pheatmap)
library(gprofiler2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(KEGGgraph)
library(pathview)
library(compEpiTools)
library(tidytext)
library(progeny)
library(dorothea)


# 1. Get input data -----------------------------------------------------------------------------

# obtain input files
MTHFD2_counts <- read_delim('counts/MTHFD2_KOvsWT_featureCount.txt', delim = "\t", comment = "#")
MTHFD2_sampleinfo <- read_delim("counts/MTHFD2_KOvsWT_ss.csv", delim = ";")

# modify counts file 
# get only counts
MTHFD2_counts <- MTHFD2_counts[, c(1, 7:18)]
colnames(MTHFD2_counts) <- c("GeneID", "K1T0R1", "K1T0R2", "K1T0R3", "K1T7R1", "K1T7R2", "K1T7R3", 
                             "WTT0R1", "WTT0R2", "WTT0R3", "WTT7R1", "WTT7R2", "WTT7R3")

# export background
background <- MTHFD2_counts %>%
    mutate(GeneID2 = str_remove(GeneID, "\\..*$")) %>%
    dplyr::select(GeneID2)
write_csv(background, "background_RNA.csv")

# modify both df so that rownames(sampleinfo) = colnames(counts)
MTHFD2_counts <- MTHFD2_counts %>% column_to_rownames(var = "GeneID") %>% as.matrix()
MTHFD2_sampleinfo <- MTHFD2_sampleinfo %>% column_to_rownames(var = "Name")

# convert chr to factor in sampleinfo df
MTHFD2_sampleinfo$Condition <- factor(MTHFD2_sampleinfo$Condition, levels = c("WT", "KO1"))
MTHFD2_sampleinfo$Time <- factor(MTHFD2_sampleinfo$Time, levels = c("initial", "final"))

# 2. Pre-filtering -----------------------------------------------------------------------------

# keep rows with cpm more than 1 in at least 2 samples
keep <- rowSums(cpm(MTHFD2_counts)>1) >= 2
MTHFD2_counts_filtered <- MTHFD2_counts[keep,]

# 3. Quality check -----------------------------------------------------------------------------

# obtain library sizes
library_sizes <- colSums(MTHFD2_counts_filtered)
library_sizes_df <- as.data.frame(library_sizes) %>%
    rownames_to_column(var = "Sample")

# plot library sizes
ggplot(library_sizes_df, aes(x = Sample, y = library_sizes)) +
    geom_bar(stat = "identity") +
    ylab("Library Sizes") +
    coord_flip()

ggsave("plots/MTHFD2_library_sizes_barplot.png")

# rlog normalization
rlog_MTHFD2_counts <- rlog(MTHFD2_counts_filtered)
rownames(rlog_MTHFD2_counts) <- rownames(MTHFD2_counts_filtered)

# obtain tidy df to plot
MTHFD2_sampleinfo_df <- MTHFD2_sampleinfo %>% rownames_to_column("Sample")

rlog_MTHFD2_counts_df_tidy <- as.data.frame(rlog_MTHFD2_counts) %>%
    rownames_to_column(var = "GeneID") %>%
    pivot_longer(names_to = "Sample", values_to = "rlog_counts", -GeneID) %>%
    left_join(MTHFD2_sampleinfo_df, by = "Sample")

# plot count distribution
ggplot(rlog_MTHFD2_counts_df_tidy, aes(x = Sample, y = rlog_counts, fill = Condition)) +
    geom_boxplot() +
    xlab("Sample") +
    ylab("rlog Counts") +
    scale_fill_discrete(labels=c("WT", "KO1")) +
    coord_flip()

ggsave("plots/MTHFD2_count_distribution_boxplot.png")

# PCA WT vs KO
# perform pca analysis
set.seed(1234)
PCA_compute_KO <- prcomp(t(rlog_MTHFD2_counts), center = T, scale. = T)

# add sample sheet data to PCA data
MTHFD2_PCA_df <- data.frame(PCA_compute_KO$x)
if (all(rownames(MTHFD2_PCA_df) == rownames(MTHFD2_sampleinfo))) {
    MTHFD2_PCA_df <- cbind(MTHFD2_PCA_df, MTHFD2_sampleinfo)
}

ggplot(MTHFD2_PCA_df, aes(x = PC1, y = PC2, colour = Condition, shape = Time)) +
    geom_point(size = 5, alpha = 0.7) +
    ggtitle("PCA RNAseq MTHFD2 KO vs WT") +
    labs(x = paste("PCA1 (",round(summary(PCA_compute_KO2)$importance[2,1]*100),"%)", sep=""), 
         y = paste("PCA2 (",round(summary(PCA_compute_KO2)$importance[2,2]*100),"%)", sep=""),
         colour = "Condition", shape = "Time") +
    theme_bw() +
    scale_color_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58")) +
    theme(legend.title = element_text(size = 11), legend.text = element_text(size = 11)) +
    theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10)) +
    theme(title = element_text(size = 11), plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste0("plots/PCA_MTHFD2_KO1vsWT.pdf"), device = "pdf", height = 3, width = 4) 

# 4. DESeq2 analysis -----------------------------------------------------------------------------

# check that our rows and columns match
all(rownames(MTHFD2_sampleinfo) == colnames(MTHFD2_counts_filtered)) # T

# create the design formula with reference WTT0
MTHFD2_sampleinfo$Comparison <- factor(MTHFD2_sampleinfo$Comparison, levels = c("WTT0", "WTT7", "KO1T0", "KO1T7"))
design_WTT0 <- as.formula(~ Comparison)

# create the DESeqDataSet object
dds_obj_MTHFD2_WTT0 <- DESeqDataSetFromMatrix(countData = MTHFD2_counts_filtered, colData = MTHFD2_sampleinfo, design = design_WTT0)

# rebuild DDS object
dds_obj_MTHFD2_WTT0 <- DESeq(dds_obj_MTHFD2_WTT0)

# see dif contrasts available
resultsNames(dds_obj_MTHFD2_WTT0)

# obtain results 
get_results_comparisons_WTT0 <- function(comparison) {
    results <- results(dds_obj_MTHFD2_WTT0, alpha = 0.05, name = comparison)
    results
}

comparisons_WTT0 <- c("Comparison_WTT7_vs_WTT0", "Comparison_KO1T0_vs_WTT0", "Comparison_KO1T7_vs_WTT0")
DESeq2_results_WTT0 <- map(comparisons_WTT0, get_results_comparisons_WTT0)

# check results
map(DESeq2_results_WTT0, summary)

# create the design formula with reference WTT7
MTHFD2_sampleinfo$Comparison <- factor(MTHFD2_sampleinfo$Comparison, levels = c("WTT7", "WTT0", "KO1T0", "KO1T7"))
design_WTT7 <- as.formula(~ Comparison)

# create the DESeqDataSet object
dds_obj_MTHFD2_WTT7 <- DESeqDataSetFromMatrix(countData = MTHFD2_counts_filtered, colData = MTHFD2_sampleinfo, design = design_WTT7)

# rebuit DDS object
dds_obj_MTHFD2_WTT7 <- DESeq(dds_obj_MTHFD2_WTT7)

# see dif contrasts available
resultsNames(dds_obj_MTHFD2_WTT7)

# obtain results 
get_results_comparisons_WTT7 <- function(comparison) {
    results <- results(dds_obj_MTHFD2_WTT7, alpha = 0.05, name = comparison)
    results
}

comparisons_WTT7 <- c("Comparison_WTT0_vs_WTT7", "Comparison_KO1T0_vs_WTT7", "Comparison_KO1T7_vs_WTT7")
DESeq2_results_WTT7 <- map(comparisons_WTT7, get_results_comparisons_WTT7)

# check results
map(DESeq2_results_WTT7, summary)

# perform shrink
perform_shrink_WTT0 <- function(coef) {
    lfcShrink(dds_obj_MTHFD2_WTT0, coef=coef, type = "apeglm")
}

comparisons_WTT0 <- c("Comparison_WTT7_vs_WTT0", "Comparison_KO1T0_vs_WTT0", "Comparison_KO1T7_vs_WTT0")
shrink_results_WTT0 <- map(comparisons_WTT0, perform_shrink_WTT0)

perform_shrink_WTT7 <- function(coef) {
    lfcShrink(dds_obj_MTHFD2_WTT7, coef=coef, type = "apeglm")
}

comparisons_WTT7 <- c("Comparison_WTT0_vs_WTT7", "Comparison_KO1T0_vs_WTT7", "Comparison_KO1T7_vs_WTT7")
shrink_results_WTT7 <- map(comparisons_WTT7, perform_shrink_WTT7)

# annotate shrinkage results
get_annotation <- function(results) {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    # set filter type
    filter_type <- "ensembl_gene_id_version"
    # set filter values
    filter_values <- rownames(results)
    # Set the list of attributes
    attribute_names <- c('ensembl_gene_id_version', 'hgnc_symbol', 'entrezgene_description', 'gene_biotype')
    # run the query
    getBM(attributes = attribute_names, filters = filter_type, values = filter_values, mart = ensembl)
}

annot_biomart_WTT0 <- map(shrink_results_WTT0, get_annotation)
annot_biomart_WTT7 <- map(shrink_results_WTT7, get_annotation)

# change column names
change_colnames <- function(df) {
    colnames(df) <- c("GeneID", "HGNC_symbol", "Description", "Biotype") 
    df
}

annot_biomart_WTT0_2 <- map(annot_biomart_WTT0, change_colnames)
annot_biomart_WTT7_2 <- map(annot_biomart_WTT7, change_colnames)

# Add annotation to the results table
add_annot <- function(results, annot) {
    annot_df <- as.data.frame(results) %>%
        rownames_to_column("GeneID") %>%
        left_join(annot, "GeneID")
    annot_df
}

DESeq2_shrink_WTT0_results_annot <- map2(shrink_results_WTT0, annot_biomart_WTT0_2, add_annot)
DESeq2_shrink_WTT7_results_annot <- map2(shrink_results_WTT7, annot_biomart_WTT7_2, add_annot)

map(DESeq2_shrink_WTT0_results_annot, nrow) # 14800
map(DESeq2_shrink_WTT7_results_annot, nrow) # 14800

# get interesting comparisons
DEG_interesting_shrink <- list(DESeq2_shrink_WTT0_results_annot[[2]], DESeq2_shrink_WTT7_results_annot[[3]])

# get significant genes with padj < 0.05 and log2FC > 0.58
get_significant_genes <- function(df) {
    df_filtered <- df %>%
        filter(padj < 0.05 & abs(log2FoldChange) >= 0.58)
    df_filtered
}

DEG_interesting_shrink_sign <- map(DEG_interesting_shrink, get_significant_genes)
map(DEG_interesting_shrink_sign, nrow) 

# write output
write_csv(DEG_interesting_shrink_sign[[1]], "DESeq2/DEG_KO1T0_vs_WTT0_shrink_annotated.csv")
write_csv(DEG_interesting_shrink_sign[[2]], "DESeq2/DEG_KO1T7_vs_WTT7_shrink_annotated.csv") 

# add top genes to dataframe
add_top_log10padj <- function(shrink_annot) {
    cutoff <- sort(shrink_annot$padj)[50]
    df <- shrink_annot %>%
        mutate(TopGeneLabel = ifelse(padj<=cutoff, HGNC_symbol, ""),
               `-log10(padj)` = -log10(padj)) 
    df
}

# interesting KO1T0 vs WTT0, KO1T7 vs WTT7, KO2T0 vs WTT0, KO2T7 vs WTT7
shrink_int_df_plot <- map(DEG_interesting_shrink, add_top_log10padj)


# plot volcano of WT vs KO - initial and final with common genes with surrounding
shrink_WT_KO1_t0 <- shrink_int_df_plot[[1]] %>%
    mutate(sign = ifelse(padj < 0.05 & abs(log2FoldChange) > 1.5, TRUE, FALSE),
           timepoint = rep("Early", nrow(shrink_int_df_plot[[1]])))
shrink_WT_KO1_t7 <- shrink_int_df_plot[[2]] %>%
    mutate(sign = ifelse(padj < 0.05 & abs(log2FoldChange) > 1.5, TRUE, FALSE),
           timepoint = rep("Late", nrow(shrink_int_df_plot[[2]])))

# add common genes
shrink_WT_KO1_t0_sign <- shrink_WT_KO1_t0 %>%
    filter(sign == TRUE) %>%
    pull(GeneID)
shrink_WT_KO1_t7_sign <- shrink_WT_KO1_t7 %>%
    filter(sign == TRUE) %>%
    pull(GeneID)
length(shrink_WT_KO1_t0_sign) # 213
length(shrink_WT_KO1_t7_sign) # 239
common_KO1 <- intersect(shrink_WT_KO1_t0_sign,shrink_WT_KO1_t7_sign)
length(intersect(shrink_WT_KO1_t0_sign,shrink_WT_KO1_t7_sign)) # 141
shrink_WT_KO1_t0_t7 <- shrink_WT_KO1_t0 %>%
    rbind(shrink_WT_KO1_t7) %>%
    mutate(common = ifelse(GeneID %in% common_KO2, "common", 
                           ifelse(sign == TRUE, "sign", "nosign")))

# put data together
ggplot(shrink_WT_KO1_t0_t7, aes(x = log2FoldChange, y=`-log10(padj)`)) + 
    geom_point(aes(colour=common, fill = sign), size=2, alpha = 0.7, shape = 21) +
    scale_fill_manual(values = c("TRUE" = "#559A73", "FALSE" = "lightgrey")) +
    scale_color_manual(values = c("common" = "black", "sign" = "#559A73", "nosign" = "lightgrey")) +
    facet_grid(~timepoint) +
    xlim(c(-10,10)) +
    ylim(c(0, 300)) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12),
          axis.text = element_text(size = 10), legend.text = element_text(size = 12),  
          strip.text = element_text(size = 12), strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/volcano_KO1vsWT_early_late.pdf", device = "pdf", width = 6.5, height =3)         

# check numbers
nrow(shrink_WT_KO1_t0[shrink_WT_KO1_t0$padj < 0.05 & shrink_WT_KO1_t0$log2FoldChange > 1.5,]) # 68
nrow(shrink_WT_KO1_t0[shrink_WT_KO1_t0$padj < 0.05 & shrink_WT_KO1_t0$log2FoldChange < -1.5,]) # 145

nrow(shrink_WT_KO1_t7[shrink_WT_KO1_t7$padj < 0.05 & shrink_WT_KO1_t7$log2FoldChange > 1.5,]) # 98
nrow(shrink_WT_KO1_t7[shrink_WT_KO1_t7$padj < 0.05 & shrink_WT_KO1_t7$log2FoldChange < -1.5,]) # 141

# get df from shrink results
get_shrink_df <- function(shrink) {
    df <- as.data.frame(shrink)
    df <- df %>%
        rownames_to_column("GeneID")
}

DESeq2_shrink_WTT0_df <- map(shrink_results_WTT0, get_shrink_df)
DESeq2_shrink_WTT7_df <- map(shrink_results_WTT7, get_shrink_df)

# check names columns to see if they are the same
map(DESeq2_shrink_WTT0_df, colnames)
map(DESeq2_shrink_WTT7_df, colnames)

# put together KO1 t0 and t7
table(DESeq2_shrink_WTT0_df[[2]]$GeneID == DESeq2_shrink_WTT7_df[[3]]$GeneID) # T

shrink_KO1_df <- DESeq2_shrink_WTT0_df[[2]] %>% 
    left_join(DESeq2_shrink_WTT7_df[[3]], by = "GeneID") %>%
    dplyr::select(GeneID, log2FoldChange.x, padj.x, log2FoldChange.y, padj.y) %>%
    mutate(FCdif = log2FoldChange.x - log2FoldChange.y,
           sign = ifelse(padj.x < 0.05 & padj.y < 0.05, "both", ifelse(padj.x < 0.05 & padj.y >= 0.05, "early", ifelse(padj.x >= 0.05 & padj.y < 0.05, "late", "NA"))))
colnames(shrink_KO1_df) <- c("GeneID", "log2FoldChange_t0", "padj_t0", "log2FoldChange_t7", "padj_t7", "FC_dif", "sign")

# select significant genes in at least one condition
shrink_KO1_df_filt <- shrink_KO1_df %>%
    filter((abs(log2FoldChange_t0) >= 0.58 & padj_t0 < 0.05) | (abs(log2FoldChange_t7) >= 0.58 & padj_t7 < 0.05))
nrow(shrink_KO1_df) # 14752
nrow(shrink_KO1_df_filt) # 3437
table(shrink_KO1_df_filt$sign) # both 2706, early 486, late 245

# get common genes of t0 and t7 and with same direction of FC
annot <- DESeq2_shrink_WTT0_results_annot[[1]] %>%
    dplyr::select(GeneID, HGNC_symbol, Description, Biotype)
nrow(annot) # 14800 -> all genes

shrink_KO1_df_filt_common <- shrink_KO1_df_filt %>%
    filter(sign == "both") %>% 
    filter((log2FoldChange_t0 > 0 & log2FoldChange_t7 > 0) | (log2FoldChange_t0 < 0 & log2FoldChange_t7 < 0)) %>%
    filter((log2FoldChange_t0 > 0.58 | log2FoldChange_t7 > 0.58) | (log2FoldChange_t0 < -0.58 | log2FoldChange_t7 < -0.58)) 
nrow(shrink_KO1_df_filt) # 3437
nrow(shrink_KO1_df_filt_common) # 3437 -> 2688

shrink_KO1_df_filt_common_annot <- shrink_KO1_df_filt_common %>%
    left_join(annot, by = "GeneID") 
length(unique(shrink_KO1_df_filt_common_annot$GeneID)) # 2688
nrow(shrink_KO1_df_filt_common_annot) # 2695

write_csv(shrink_KO1_df_filt_common_annot, "DESeq2/DEG_KO1_vs_WT_t0t7_common_annotated_shrink.csv")

# 5. GO BP analysis -----------------------------------------------------------------------------

# get DEG
DEG_KO1 <- read_delim("DESeq2/DEG_KO1_vs_WT_t0t7_common_annotated_shrink.csv", delim = ",")

# number of genes
nrow(DEG_KO1) # 2695

# keep only protein coding
keep_protein_coding <- function(df) {
    df2 <- df %>%
        filter(Biotype == "protein_coding")
    df2
}
DEG_protein_coding <- keep_protein_coding(DEG_KO1)

# number of genes
nrow(DEG_protein_coding)

# duplicates
table(duplicated(DEG_protein_coding$GeneID)) 

# check manually eliminate duplicates
DEG_protein_coding[duplicated(DEG_protein_coding$GeneID),]

remove_duplicates <- function(df) {
    df2 <- df %>%
        filter(!(Description %in% c("ADAMTSL4 antisense RNA 1", "uncharacterized LOC105374836", "uncharacterized LOC105377067", 
                                    "STT3A antisense RNA 1", "uncharacterized LOC105369632", "uncharacterized LOC100506403",
                                    "TBC1D7-LOC100130357 readthrough")))
    df2
}
DEG_protein_coding_unique <- remove_duplicates(DEG_protein_coding)

# check duplicates
table(duplicated(DEG_protein_coding_unique$GeneID)) # F

# get list of UP and DOWN
DEG_KO1_protein_up <- DEG_protein_coding_unique %>%
    filter(log2FoldChange_t0 > 0 & log2FoldChange_t7 > 0)
nrow(DEG_KO1_protein_up) # 1009

DEG_KO1_protein_down <- DEG_protein_coding_unique %>%
    filter(log2FoldChange_t0 < 0 & log2FoldChange_t7 < 0)
nrow(DEG_KO1_protein_down) # 1049

# background
background <- MTHFD2_counts$Geneid

# up and down genes - protein coding
DEG_KO1_protein_up_gene <- DEG_KO1_protein_up$GeneID
DEG_KO1_protein_down_gene <- DEG_KO1_protein_down$GeneID

# remover version
ensID_rm_version <- function(x) {
    for (i in seq_along(x)) {
        x[i] <- str_split(x[i], "\\.")[[1]][1]
    }
    x
}
KO1_up_nov <- ensID_rm_version(DEG_KO1_protein_up_gene)
KO1_down_nov <- ensID_rm_version(DEG_KO1_protein_down_gene)
background_nov <- ensID_rm_version(background)

# function for GO analysis
obtainGOBP <- function(genes) {
    enrichGO(
        gene = genes,
        universe = background_nov,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.10,
        readable      = TRUE)
}

# function to obtain df to plot
results_plot <- function(df) {
    df_r <- df %>%
        mutate(`-log10(qvalue)` = -log10(qvalue)) %>%
        arrange(desc(`-log10(qvalue)`)) %>%
        filter(qvalue < 0.05)
    df_r$qvalue <- formatC(df_r$qvalue, format = "e", digits = 2)
    df_r$qvalue <- as.double(df_r$qvalue)
    df_r
}

KO1_up_results <- obtainGOBP(KO1_up_nov)
KO1_up_results_simp <- simplify(KO1_up_results)
KO1_up_results_simp_result <- results_plot(KO1_up_results_simp@result)

KO1_down_results <- obtainGOBP(KO1_down_nov)
KO1_down_results_simp <- simplify(KO1_down_results)
KO1_down_results_simp_result <- results_plot(KO1_down_results_simp@result)

# put together up and down
KO1_up_results_filter <- KO1_up_results_simp_result[c(1:13),] 
KO1_up_results_add <- KO1_up_results_filter %>%
    mutate(DE = rep("Up", nrow(KO1_up_results_filter)),
           GeneRatio_real = Count/313)
KO1_down_results_filter <- KO1_down_results_simp_result[c(1:13),] 
KO1_down_results_add <- KO1_down_results_filter %>%
    mutate(DE = rep("Down", nrow(KO1_down_results_filter)),
           GeneRatio_real = Count/531)

table(colnames(KO1_up_results_add) == colnames(KO1_down_results_add))
KO1_up_down_simp <- rbind(KO1_up_results_add, KO1_down_results_add)

# plot
ggplot(KO1_up_down_simp, aes(reorder(Description, -`-log10(qvalue)`), `-log10(qvalue)`, color = DE, size = GeneRatio_real)) +
    geom_point(alpha = 1, stroke = 0) +
    facet_wrap(~DE, scales = "free", nrow = 1) +
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "lightgrey", size = 0.3) +
    labs(x = "", col = "") +
    scale_color_manual(values = c("Down" = "#aaafb0", "Up" = "#e56e58")) +
    coord_flip() +
    theme_bw() +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8)) +
    theme(axis.title = element_text(size = 10), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 8)) +
    theme(title = element_text(size = 11), plot.title = element_text(hjust = 0.5),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          strip.text = element_text(size = 10), 
          strip.background = element_rect(colour="black", fill="#EDEDED"))
ggsave("plots/GOBP_KOcommon_plots.pdf", device = "pdf", height = 2.5, width = 9)

# 5. Heatmap centromere-kinetochore genes  -----------------------------------------------------------------------------

# obtain input files
MTHFD2_counts <- read_delim('counts/MTHFD2_KOvsWT_featureCount.txt', delim = "\t", comment = "#")
MTHFD2_sampleinfo <- read_delim("counts/MTHFD2_KOvsWT_ss.csv", delim = ";")

# modify counts file - get only counts
# get only WT and KO2 data
MTHFD2_counts <- MTHFD2_counts[, c(1, 7:18)]
colnames(MTHFD2_counts) <- c("GeneID", "K1T0R1", "K1T0R2", "K1T0R3", "K1T7R1", "K1T7R2", "K1T7R3", 
                             "WTT0R1", "WTT0R2", "WTT0R3", "WTT7R1", "WTT7R2", "WTT7R3")

# convert to matrix
MTHFD2_counts_m <- MTHFD2_counts %>%
    column_to_rownames("GeneID") %>%
    as.matrix()

# normalise counts
MTHFD2_counts_norm <- rlog(MTHFD2_counts_m)

# heatmap of genes involved in centromeres
# get ens ID from genes
symbol_centromeres <- c("CENPA", "CENPC", "CENPN", "CENPT", "CENPU", "CENPH", "CENPM", "CENPS",
                        "CENPR", "CENPO", "CENPI", "CENPL", "CENPP", "CENPK", "CENPF", "CENPE",
                        "INCENP", "AURKB", "KIF2C", "NDC80", "MIS12", "KNL1", "PLK1", "ZWINT", "TTK",
                        "KNTC1", "BUB1", "BUB3", "BUB1B", "MAD1L1", "MAD2L1", "CDC20")

length(symbol_centromeres) == length(unique(symbol_centromeres)) # T
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensID_centromeres <- getBM(attributes = c("ensembl_gene_id_version", "hgnc_symbol"), 
                           filters = "hgnc_symbol", values = symbol_centromeres, mart = ensembl)

# get counts of interesting genes
centromeres_counts_norm_df <- MTHFD2_counts_norm %>%
    as.data.frame() %>%
    rownames_to_column("GeneID") %>%
    filter(GeneID %in% ensID_centromeres$ensembl_gene_id_version) %>%
    left_join(ensID_centromeres, by = c("GeneID" = "ensembl_gene_id_version")) %>%
    dplyr::select(-GeneID) %>%
    arrange(hgnc_symbol)
nrow(centromeres_counts_norm_df) # 28

centromeres_counts_norm_m <- centromeres_counts_norm_df %>%
    column_to_rownames("hgnc_symbol") %>%
    as.matrix()

# Get some nice colours
mypalette <- brewer.pal(11, "RdBu")
morecols <- colorRampPalette(c("navy", "white", "red"))(50)

# annotation
MTHFD2_sampleinfo_annot2 <- MTHFD2_sampleinfo %>%
    column_to_rownames(var = "Name") %>%
    filter(Condition %in% c("WT", "KO2")) %>%
    dplyr::select(-c(Comparison,Replicate))
ann_colors <- list(Condition = c("WT" = "#aaafb0", "KO1" = "#e56e58"), Time = c("initial" = "#86dff2", "final" = "#5998f0"))

# heatmap with clustering to obtain the order of the rows (genes)
set.seed(1234)
pdf(file= paste0("plots/heatmap_KO1vsWT_centromeres.pdf"), width = 5, height = 4)
pheatmap(centromeres_counts_norm_m,
         scale = "row",
         show_rownames =  T, 
         show_colnames =  T, 
         border_color = NA,
         cluster_rows = T, cluster_cols = F, color = morecols, 
         annotation_col = MTHFD2_sampleinfo_annot2, annotation_colors = ann_colors,
         na_col = "whitesmoke", fontsize = 9, fontsize_row = 8)
dev.off()

# check if differentially expressed
DEG_KO1_centromeric_genes <- DEG_KO1 %>%
    filter(HGNC_symbol %in% symbol_centromeres )
DEG_KO1_centromeric_genes
