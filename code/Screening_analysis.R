# libraries
library(tidyverse)
library(MAGeCKFlute)
library(pheatmap)
library(RColorBrewer)

# Objective: Analyse whole-genome CRISPR genetic screening MTHFD2 WT vs KO cells

# 1. Get input data -------------------------------------------------------

# open count files
brunello_MTHFD2_count <- read.delim("mageck_count/brunello_MTHFD2.count.txt", sep = "\t")
head(brunello_MTHFD2_count)

# open summary
brunello_MTHFD2_summary <- read.delim("mageck_count/brunello_MTHFD2.countsummary.txt", sep = "\t")
head(brunello_MTHFD2_summary)

# 2. Quality control of screening --------------------------------------------------------

# mapping eficiency
mapping_info <- brunello_MTHFD2_summary[, c("Label", "Reads", "Mapped")]
mapping_info <- mapping_info %>%
    mutate(Sample = case_when(
        Label == "WT_day0" ~ "WT0",
        Label == "WT_dayf" ~ "WTF",
        Label == "K1_day0" ~ "K10",
        Label == "K1_dayf" ~ "K1F",
        Label == "K2_day0" ~ "K20",
        Label == "K2_dayf" ~ "K2F"),
        Unmapped = Reads-Mapped,
        Mapped_perc = Mapped/Reads * 100,
        Unmapped_perc = Unmapped/Reads * 100)

mapping_numbers <- mapping_info[, c("Sample", "Mapped", "Unmapped")] %>%
    pivot_longer(-Sample, names_to = "Mapping", values_to = "Reads")

mapping_percentage <- mapping_info[, c("Sample", "Mapped_perc", "Unmapped_perc")] %>%
    pivot_longer(-Sample, names_to = "Mapping", values_to = "Percentage")

mapping_numbers$Percentage <- mapping_percentage$Percentage

mapping_numbers$Sample <- factor(mapping_info_plot$Sample, levels = c("WT0",  "K10", "K20", "WTF", "K1F","K2F"))
mapping_numbers$Mapping <- factor(mapping_info_plot$Mapping, levels = c("Unmapped", "Mapped"))

ggplot(mapping_numbers, aes(x = Sample, y = Reads, fill = Mapping, label = paste0(sprintf("%0.2f", round(Percentage, digits = 2)), "%"))) +
    geom_col(width = 0.75) +
    ggtitle("Mapping Ratio") +
    scale_fill_manual(values = c("Unmapped" = "#B3E1EB", "Mapped" = "#51B0C4")) + 
    geom_text(size = 3, position = position_stack(vjust = 0.5))  +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5))
ggsave("plots/mapping_efficiency.pdf")

# check Gini index
gini_index_info <- brunello_MTHFD2_summary[, c("Label", "GiniIndex")]

gini_index_info <- gini_index_info %>%
    mutate(Sample = case_when(
        Label == "WT_day0" ~ "WT0",
        Label == "WT_dayf" ~ "WTF",
        Label == "K1_day0" ~ "K10",
        Label == "K1_dayf" ~ "K1F",
        Label == "K2_day0" ~ "K20",
        Label == "K2_dayf" ~ "K2F"),
        Condition = case_when(
            str_detect(Label, "WT") ~ "WT",
            str_detect(Label, "K1") ~ "K1",
            str_detect(Label, "K2") ~ "K2"))

gini_index_info$Sample <- factor(gini_index_info$Sample, levels = c("WT0",  "K20", "K10", "WTF", "K2F","K1F"))
gini_index_info$Condition <- factor(gini_index_info$Condition, levels = c("WT", "K1", "K2"))

ggplot(gini_index_info, aes(x = Sample, y = GiniIndex, fill = Condition)) +
    geom_col(width = 0.5) +
    ggtitle("Evenness of sgRNA reads") +
    ylab("Gini Index") +
    scale_fill_manual(values = c("WT" = "#aaafb0", "K1" = "#e56e58", "K2" = "#f9beb3")) +
    theme_classic() + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12))
ggsave("plots/gini_index.pdf", device = "pdf", width = 4.5, height = 3)

# check sgrna with 0 counts
WT0_guides0 <- sum(brunello_MTHFD2_count$WT_day0 == 0)
WTF_guides0 <- sum(brunello_MTHFD2_count$WT_dayf == 0)
K10_guides0 <- sum(brunello_MTHFD2_count$K1_day0 == 0)
K1F_guides0 <- sum(brunello_MTHFD2_count$K1_dayf == 0)
K20_guides0 <- sum(brunello_MTHFD2_count$K2_day0 == 0)
K2F_guides0 <- sum(brunello_MTHFD2_count$K2_dayf == 0)

# create data frame
info_sg_0 <- data.frame(Samples = c("WT0","K10", "K20", "WTF", "K1F", "K2F"), 
                        sg_0 = c(WT0_guides0, K10_guides0, K20_guides0, 
                                 WTF_guides0, K1F_guides0, K2F_guides0))
info_sg_0$Samples <- factor(info_sg_0$Samples, levels = c("WT0", "K10", "K20", "WTF", "K1F", "K2F"))
info_sg_0 <- info_sg_0 %>%
    mutate(log10_missed_guides = log10(sg_0)) %>%
    mutate(Perc_missed_guides = (sg_0/77442)*100) %>%
    mutate(Condition = case_when(
        str_detect(Samples, "WT") ~ "WT",
        str_detect(Samples, "K1") ~ "K1",
        str_detect(Samples, "K2") ~ "K2"
    ))

# plot
ggplot(info_sg_0, aes(x = Samples, y = log10_missed_guides, fill = Condition)) +
    geom_col(width = 0.75) +
    ylab("log10 missing sgRNAs") +
    ggtitle("Missed sgRNAs") +
    scale_fill_manual(values = c("WT" = "#aaafb0", "K1" = "#f9beb3", "K2" = "#e56e58")) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5))
ggsave("plots/missed_sgRNAs.pdf")

ggplot(info_sg_0, aes(x = Samples, y = Perc_missed_guides, fill = Condition)) +
    geom_col(width = 0.75) +
    ylab("% missing sgRNAs") +
    ggtitle("Missed sgRNAs") +
    scale_fill_manual(values = c("WT" = "#aaafb0", "K1" = "#f9beb3", "K2" = "#e56e58")) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5))
ggsave("plots/missed_sgRNAs_perc.pdf")

# check genes with 0 sgRNAs
brunello_MTHFD2_genes <- brunello_MTHFD2_count %>%
    group_by(Gene) %>%
    summarise(WT_day0_sum = sum(WT_day0),
              K1_day0_sum = sum(K1_day0),
              K2_day0_sum = sum(K2_day0),
              WT_dayf_sum = sum(WT_dayf),
              K1_dayf_sum = sum(K1_dayf),
              K2_dayf_sum = sum(K2_dayf))

WT0_genes0 <- sum(brunello_MTHFD2_genes$WT_day0_sum == 0)
WTF_genes0 <- sum(brunello_MTHFD2_genes$WT_dayf_sum == 0)
K10_genes0 <- sum(brunello_MTHFD2_genes$K1_day0_sum == 0)
K1F_genes0 <- sum(brunello_MTHFD2_genes$K1_dayf_sum == 0)
K20_genes0 <- sum(brunello_MTHFD2_genes$K2_day0_sum == 0)
K2F_genes0 <- sum(brunello_MTHFD2_genes$K2_dayf_sum == 0)

# create data frame
info_genes_0 <- data.frame(Samples = c("WT0", "K10", "K20", "WTF", "K1F", "K2F"), 
                           sg_0 = c(WT0_genes0, K10_genes0, K20_genes0, 
                                    WTF_genes0, K1F_genes0, K2F_genes0))
info_genes_0$Samples <- factor(info_genes_0$Samples, levels = c("WT0", "K10", "K20", "WTF", "K1F", "K2F"))

# 0 -> so not plot

# plot distribution of counts
WT0_guides_per_gene <- brunello_MTHFD2_count %>%
    dplyr::select(sgRNA, Gene, WT_day0)
WT0_mean <- mean(WT0_guides_per_gene$WT_day0)
K10_guides_per_gene <- brunello_MTHFD2_count %>%
    dplyr::select(sgRNA, Gene, K1_day0)
K10_mean <- mean(K10_guides_per_gene$K1_day0)
K20_guides_per_gene <- brunello_MTHFD2_count %>%
    dplyr::select(sgRNA, Gene, K2_day0)
K20_mean <- mean(K20_guides_per_gene$K2_day0)

ggplot(WT0_guides_per_gene, aes(x = WT_day0)) +
    geom_histogram(binwidth = 100, fill = "#aaafb0", col="black") +
    geom_vline(xintercept = WT0_mean) + 
    ggtitle("Count distribution - WT0") +
    ylab("n") +
    xlab("hits") +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.25))
ggsave("plots/WT_count_hist_distribution.pdf", device = "pdf", width = 5, height = 4)

ggplot(K10_guides_per_gene, aes(x = K1_day0)) +
    geom_histogram(binwidth = 100, fill = "#e56e58", col="black") +
    geom_vline(xintercept = K10_mean) + 
    ggtitle("Count distribution - K10") +
    ylab("n") +
    xlab("hits") +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.25))
ggsave("plots/K1_count_hist_distribution.pdf", device = "pdf", width = 5, height = 4)

ggplot(K20_guides_per_gene, aes(x = K2_day0)) +
    geom_histogram(binwidth = 100, fill = "#f9beb3", col="black") +
    geom_vline(xintercept = K20_mean) + 
    ggtitle("Count distribution - K20") +
    ylab("n") +
    xlab("hits") +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.25))
ggsave("plots/K2_count_hist_distribution.pdf", device = "pdf", width = 5, height = 4)

# 3. Analyse MLE results with MAGeCKFlute -------------------------------------------------------

# Analysis of MLE results putting together WT and KO
FluteMLE("mageck_mle/MLE_WT_K1.gene_summary.txt", treatname = "K1", ctrlname = "WT", incorporateDepmap = TRUE, proj = "FluteMLE_WT_K1")
FluteMLE("mageck_mle/MLE_WT_K2.gene_summary.txt", treatname = "K2", ctrlname = "WT", incorporateDepmap = TRUE, proj = "FluteMLE_WT_K2")

# 4. Check cell cycle normalization  -------------------------------------------------------

# read files before normalization
mle_WT <- read_delim("mageck_mle/MLE_WT.gene_summary.txt", delim = "\t")
mle_K1 <- read_delim("mageck_mle/MLE_K1.gene_summary.txt", delim = "\t")
mle_K2 <- read_delim("mageck_mle/MLE_K2.gene_summary.txt", delim = "\t")

# read files after normalization
WT_K1_flute <- read_delim("MAGeCKFlute_FluteMLE_WT_K1/QC/FluteMLE_WT_K1_processed_data.txt", delim = "\t")
WT_K2_flute <- read_delim("MAGeCKFlute_FluteMLE_WT_K2/QC/FluteMLE_WT_K2_processed_data.txt", delim = "\t")

# get essential genes from MAGeCKFlute
essent_genes_data <- read_delim("brunello/MAGeCKFlute_essential_genes.txt", delim = "\t")

# scores before normalization
mle_WT_allessential_score <- mle_WT %>%
    filter(Gene %in% essent_genes_data$hgnc_symbol) %>%
    dplyr::select(Gene, `HCT116|beta`) 
colnames(mle_WT_allessential_score)[2] <- c("WT")

mle_K1_allessential_score <- mle_K1 %>%
    filter(Gene %in% essent_genes_data$hgnc_symbol) %>%
    dplyr::select(Gene, `HCT116|beta`) 
colnames(mle_K1_allessential_score)[2] <- c("K1")

mle_K2_allessential_score <- mle_K2 %>%
    filter(Gene %in% essent_genes_data$hgnc_symbol) %>%
    dplyr::select(Gene, `HCT116|beta`) 
colnames(mle_K2_allessential_score)[2] <- c("K2")

# put together WT and KO1, and WT and KO2
mle_WT_KO1_allessential_score <- mle_K1_allessential_score %>%
    left_join(mle_WT_allessential_score) %>%
    mutate(Normalization = rep("before", nrow(mle_K1_allessential_score)))

mle_WT_KO2_allessential_score <- mle_K2_allessential_score %>%
    left_join(mle_WT_allessential_score) %>%
    mutate(Normalization = rep("before", nrow(mle_K2_allessential_score)))

# scores after normalization
WT_K1_flute_allessential <- WT_K1_flute %>%
    filter(Gene %in% essent_genes_data$hgnc_symbol) %>%
    dplyr::select(Gene, K1, WT)

WT_K1_flute_allessential <- WT_K1_flute_allessential %>%
    mutate(Normalization = rep("after", nrow(WT_K1_flute_allessential)))

WT_K2_flute_allessential <- WT_K2_flute %>%
    filter(Gene %in% essent_genes_data$hgnc_symbol) %>%
    dplyr::select(Gene, K2, WT)

WT_K2_flute_allessential <- WT_K2_flute_allessential %>%
    mutate(Normalization = rep("after", nrow(WT_K2_flute_allessential)))

# put data together
WT_K1_allessential_scores <- mle_WT_KO1_allessential_score %>%
    rbind(WT_K1_flute_allessential) 

WT_K2_allessential_scores <- mle_WT_KO2_allessential_score %>%
    rbind(WT_K2_flute_allessential) 

# plot
ggplot(WT_K1_allessential_scores, aes(x=WT, y=K1, color = Normalization)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_color_manual(values = c("#990099", "#FF66FF")) +
    theme_classic() +
    theme(
        axis.title = element_text(size = 13), axis.text = element_text(size = 12))
ggsave("plots/cellcycle_normal_WTK1_allessential.pdf", device = "pdf", width = 5, height = 3.8)

ggplot(WT_K2_allessential_scores, aes(x=WT, y=K2, color = Normalization)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_color_manual(values = c("#990099", "#FF66FF")) +
    theme_classic() +
    theme(
        axis.title = element_text(size = 13), axis.text = element_text(size = 12))
ggsave("plots/cellcycle_normal_WTK2_allessential.pdf", device = "pdf", width = 5, height = 3.8)

# 5. Get top hits  -----------------------------------------------------------------------------------------

# read files
WT_K1_flute <- read_delim("MAGeCKFlute_FluteMLE_WT_K1/QC/FluteMLE_WT_K1_processed_data.txt", delim = "\t")
WT_K2_flute <- read_delim("MAGeCKFlute_FluteMLE_WT_K2/QC/FluteMLE_WT_K2_processed_data.txt", delim = "\t")

# get intersting info
WT_K1_flute_int <- WT_K1_flute[, c("Gene", "WT", "K1", "Diff")]
WT_K1_flute_int_2 <- WT_K1_flute_int %>%
    arrange(Diff) %>%
    mutate(Condition = rep("KO1", nrow(WT_K1_flute_int)),
           Rank = seq(1: nrow(WT_K1_flute_int)))
colnames(WT_K1_flute_int_2)[3] <- "KO"
WT_K2_flute_int <- WT_K2_flute[, c("Gene", "WT", "K2", "Diff")]
WT_K2_flute_int_2 <- WT_K2_flute_int %>%
    arrange(Diff) %>%
    mutate(Condition = rep("KO2", nrow(WT_K2_flute_int)),
           Rank = seq(1: nrow(WT_K2_flute_int)))  
colnames(WT_K2_flute_int_2)[3] <- "KO"

# common significant genes
K1_WT_top <- WT_K1_flute_int$Gene[WT_K1_flute_int$Diff > 1]
K2_WT_top <- WT_K2_flute_int$Gene[WT_K2_flute_int$Diff > 1]
K1_K2_common_top <- intersect(K1_WT_top, K2_WT_top) 

K1_WT_bottom <- WT_K1_flute_int$Gene[WT_K1_flute_int$Diff < -1]
K2_WT_bottom <- WT_K2_flute_int$Gene[WT_K2_flute_int$Diff < -1]
K1_K2_common_bottom <- intersect(K1_WT_bottom, K2_WT_bottom) 

# put together both
WT_K1_K2_flute_int <- rbind(WT_K1_flute_int_2, WT_K2_flute_int_2)

# label common top up and down
# up TOP2A, GP1BA, NRM, HPSE, SPATA5L1,
# down XRN2, QTRT1, HSPD1, LIAS, WAPAL
genes_label <- c("TOP2A", "GP1BA", "NRM", "HPSE", "SPATA5L1",
                 "XRN2", "QTRT1", "HSPD1", "LIAS", "WAPAL")

# add significance
WT_K1_K2_flute_int_sign <- WT_K1_K2_flute_int %>%
    mutate(sign = ifelse(Diff > 1, "up", ifelse(Diff < -1, "down", "no")),
           common = ifelse(Gene %in% K1_K2_common_top | Gene %in% K1_K2_common_bottom, "both", ""),
           label = ifelse(Gene %in% genes_label, Gene, ""))

# plot beta difference vs ranl
ggplot(WT_K1_K2_flute_int_sign, aes(label = label, x = Rank, y = Diff)) +
    geom_point(data = WT_K1_K2_flute_int_sign %>% subset(common == ""), aes(x = Rank, y = Diff, fill = sign), 
               shape = 21, stroke = NA, size = 3) +
    geom_point(data = WT_K1_K2_flute_int_sign %>% subset(common == "both"), aes(x = Rank, y = Diff, fill = sign), 
               shape = 21, color = "black", size = 3, stroke = 0.5) +
    facet_wrap(~Condition) +
    ylab("beta(KO) - beta(WT)") +
    geom_text_repel(size = 3, max.overlaps = 400) +
    scale_fill_manual(values = c("up" = "#4343FF", "down" = "red", "no" = "#C0C0C0")) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_text(size = 14), axis.text = element_text(size = 12),
          strip.text = element_text(size = 14), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          title = element_text(size = 16), plot.title = element_text(hjust = 0.5))
ggsave("plots/beta_dif_rank_KO1_KO2.pdf", device = "pdf", width = 7, height = 4)

