# libraries
library(tidyverse)
library(ggpubr)
library(gtools)
library(rstatix)

# 1. Read input ------------------------------------------------------------------------------------------------------------------

# get RNA expression
expression_MTHFD2 <- read.csv("Expression_Public_22Q4_subsetted.csv")
head(expression_MTHFD2)
str(expression_MTHFD2)
colnames(expression_MTHFD2) <- c("X", "MTHFD2_RNA")

# get CRISPR data
CRISPR_MTHFD2 <- read.csv("CRISPR_(DepMap_Public_22Q4+Score,_Chronos)_subsetted.csv")
head(CRISPR_MTHFD2)
str(CRISPR_MTHFD2)
colnames(CRISPR_MTHFD2) <- c("X", "MTHFD2_CRISPR")

# get all proteomics
proteomics_data <- read.csv("Proteomics.csv")
head(proteomics_data)
str(proteomics_data)

# filter proteomics MTHFD2
proteomics_MTHFD2 <- proteomics_data[, c("X", "MTHFD2..P13995.")]
head(proteomics_MTHFD2)
str(proteomics_MTHFD2)
colnames(proteomics_MTHFD2) <- c("X", "MTHFD2_proteomics")

# get aneuplody info
data_aneuploidy <- read.csv("NIHMS1706666-supplement-Supplementary_Table_1.csv")
head(data_aneuploidy)
str(data_aneuploidy)

# put together all data
data_all <- data_aneuploidy %>% 
    left_join(expression_MTHFD2, by = c("DepMap_ID" = "X")) %>%
    left_join(CRISPR_MTHFD2, by = c("DepMap_ID" = "X")) %>%
    left_join(proteomics_MTHFD2, by = c("DepMap_ID" = "X")) 
head(data_all)
str(data_all)

# check NA
nrow(data_all) # 997
length(unique(data_all$DepMap_ID)) # 997
table(is.na(data_all$Aneuploidy_score))
table(is.na(data_all$Aneuploidy_group_classification))
table(is.na(data_all$MTHFD2_RNA)) # data of 55 missing
table(is.na(data_all$MTHFD2_proteomics)) # data of 626 missing
table(is.na(data_all$MTHFD2_CRISPR)) # data of 355 missing

#  aneuploidy factor
data_all$Aneuploidy_group_classification <- factor(data_all$Aneuploidy_group_classification,
                                                   c("Near-Euploid", "Intermediate", "Highly-aneuploid"))

# get short name
short_name <- str_split(data_all$CCLE_name, "_")
short_name <- unlist(map(short_name, ~ .x[1]))
length(short_name)

data_all$short_name <- short_name

# 2. Analyse MTHFD2 RNA, protein and CRISPR in specific cell lines ----------------------------------------------------------------

# get cell lines of interest
cells <- c("MCF7_BREAST", "T47D_BREAST", "SKBR3_BREAST", "MDAMB231_BREAST", "BT549_BREAST",
           "NCIH358_LUNG", "A549_LUNG", "NCIH1437_LUNG", "EBC1_LUNG", "NCIH226_LUNG",
           "HCT116_LARGE_INTESTINE", "RKO_LARGE_INTESTINE", "HT29_LARGE_INTESTINE",
           "SW620_LARGE_INTESTINE", "SW480_LARGE_INTESTINE")

# filter
data_cells <- data_all %>%
    filter(CCLE_name %in% cells) %>%
    mutate(type = case_when(
        str_detect(CCLE_name, "BREAST") ~ "Breast",
        str_detect(CCLE_name, "LUNG") ~ "Lung",
        str_detect(CCLE_name, "INTESTINE") ~ "Colon",
    ))
data_cells$type <- factor(data_cells$type, levels = c("Breast", "Lung", "Colon"))

# plot RNA levels
ggplot(data_cells, aes(x=reorder(short_name, -MTHFD2_RNA), y=MTHFD2_RNA, fill = type)) +
    geom_bar(width = 0.60, stat="identity") +
    facet_wrap(~type, scales = "free_y", ncol = 1) +
    scale_fill_manual(values=c("Breast" = "#E27AD8", "Lung" = "#669DDB", "Colon" = "#DFA247")) +
    labs(x = "", y = "MTHFD2 RNA Levels", fill = "Cancer Type") +
    theme_bw() +
    coord_flip() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none")
ggsave("plots/cancer_cells_WB_RNA_MTHFD2.pdf", device = "pdf", width = 3, height = 6)

# plot protein levels
ggplot(data_cells, aes(x=reorder(short_name, -MTHFD2_proteomics), y=MTHFD2_proteomics, fill = type)) +
    geom_bar(width = 0.60, stat="identity") +
    facet_wrap(~type, scales = "free_y", ncol = 1) +
    scale_fill_manual(values=c("Breast" = "#E27AD8", "Lung" = "#669DDB", "Colon" = "#DFA247")) +
    labs(x = "", y = "MTHFD2 Protein Levels", fill = "Cancer Type") +
    theme_bw() +
    coord_flip() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none")
ggsave("plots/cancer_cells_WB_protein_MTHFD2.pdf", device = "pdf", width = 3, height = 6)

# plot CRISPR values
ggplot(data_cells, aes(x=reorder(short_name, -MTHFD2_CRISPR), y=MTHFD2_CRISPR, fill = type)) +
    geom_bar(width = 0.60, stat="identity") +
    facet_wrap(~type, scales = "free_y", ncol = 1) +
    labs(x = "", y = "MTHFD2 Essentiality", fill = "Cancer Type") +
    scale_fill_manual(values=c("Breast" = "#E27AD8", "Lung" = "#669DDB", "Colon" = "#DFA247")) +
    theme_bw() +
    coord_flip() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none")
ggsave("plots/cancer_cells_WB_CRISPR_MTHFD2.pdf", device = "pdf", width = 3, height = 6)

# 3. Analyse relationship between MTHFD2 and ploidy ----------------------------------------------------------------

# RNA expression vs aneuploidy class

# filter cell lines with RNA info
data_all_RNA <- data_all %>%
    filter(!is.na(MTHFD2_RNA))
nrow(data_all_RNA) # 942

# check duplicates names
table(duplicated(data_all_RNA$short_name))
table(duplicated(data_all_RNA$CCLE_name))
data_all_RNA$short_name[duplicated(data_all_RNA$short_name)]

# statistics -  not normal - wilcox
aneuploidy_low <- data_all_RNA$MTHFD2_RNA[data_all_RNA$Aneuploidy_group_classification == "Near-Euploid"]
aneuploidy_med <- data_all_RNA$MTHFD2_RNA[data_all_RNA$Aneuploidy_group_classification == "Intermediate"]
aneuploidy_high <- data_all_RNA$MTHFD2_RNA[data_all_RNA$Aneuploidy_group_classification == "Highly-aneuploid"]

shapiro.test(aneuploidy_low)
shapiro.test(aneuploidy_med)
shapiro.test(aneuploidy_high)

# check RNA expression level in different categories
my_comp <- list(c("Near-Euploid", "Intermediate"), c("Near-Euploid", "Highly-aneuploid"), 
                c("Intermediate", "Highly-aneuploid"))

ggplot(data_all_RNA, aes(x= Aneuploidy_group_classification, y = MTHFD2_RNA,
                         fill = Aneuploidy_group_classification)) +
    geom_boxplot(width = 0.5, fatten = 3) +
    labs(x = "", y = "MTHFD2 RNA expression", fill="Aneuploidy group") +
    stat_compare_means(comparisons = my_comp, method = "wilcox") +
    scale_fill_manual(values = c("#FBB13C", "#FE6847", "#B66D0D")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(size = 12), axis.text = element_text(size = 11))
ggsave("plots/RNA_MTHFD2_aneuploidy_group_nice.pdf", device = "pdf", width = 4.2, height = 4)

# compare RNA expression top/bottom with aneuplody score

# RNA top bottom
MTH_RNA_high <- data_all_RNA %>%
    filter(MTHFD2_RNA >= quantile(MTHFD2_RNA, 0.75))
MTH_RNA_high <- MTH_RNA_high %>%
    mutate(Level = rep("High", nrow(MTH_RNA_high)))
MTH_RNA_low <- data_all_RNA %>%
    filter(MTHFD2_RNA <= quantile(MTHFD2_RNA, 0.25))
MTH_RNA_low <- MTH_RNA_low %>%
    mutate(Level = rep("Low", nrow(MTH_RNA_low)))
high_low_RNA <- rbind(MTH_RNA_high, MTH_RNA_low)
nrow(high_low_RNA) # 1194
high_low_RNA$Level <- factor(high_low_RNA$Level, levels = c("Low", "High"))

# plot aneuploidy score 
ggplot(high_low_RNA, aes(x= Level, y = Aneuploidy_score)) +
    geom_violin(aes(fill = Level), col = NA, show.legend = T) +
    geom_boxplot(width = 0.3, alpha = 0, fatten = 3) +
    ylab("Aneuploidy Score") +
    scale_fill_manual(values = c("High" = "#0B7E8A", "Low" = "#7FC2C9")) +
    stat_compare_means(method = "wilcox", comparisons = list(c("Low","High"))) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(size = 12), axis.text = element_text(size = 11))
ggsave("plots/RNA_MTHFD2_aneuploidy_score.pdf", device = "pdf", width = 2, height = 3)

# compare MTHFD2 essentiality with aneuploidy class

# get high and low individuals in each aneuploidy type
data_aneuploidy <- split(data_all_RNA, data_all_RNA$Aneuploidy_group_classification)

get_high_MTH <- function(df) {
    df_MTHFD2_high <- df %>% 
        filter(MTHFD2_RNA >= quantile(MTHFD2_RNA, 0.75))
}
data_aneuploidy_high <- map(data_aneuploidy, get_high_MTH)

get_low_MTH <- function(df) {
    df_MTHFD2_low <- df %>% 
        filter(MTHFD2_RNA <= quantile(MTHFD2_RNA, 0.25))
}
data_aneuploidy_low <- map(data_aneuploidy, get_low_MTH)

# put together datasets
high_each_aneu <- do.call(rbind, data_aneuploidy_high)
high_each_aneu <- high_each_aneu %>%
    mutate(Level = rep("High", nrow(high_each_aneu)))

low_each_aneu <- do.call(rbind, data_aneuploidy_low)
low_each_aneu <- low_each_aneu %>%
    mutate(Level = rep("Low", nrow(low_each_aneu)))

nrow(high_each_aneu) # 96
nrow(low_each_aneu) # 96
high_each_aneu_RNA_all <- rbind(high_each_aneu, low_each_aneu)
nrow(high_each_aneu_RNA_all) # 192
high_each_aneu_RNA_all$Level <- factor(high_each_aneu_RNA_all$Level, levels = c("Low", "High"))

# plot only high expression
high_each_aneu_RNA_all_high <- high_each_aneu_RNA_all[high_each_aneu_RNA_all$Level == "High",]

# statistics -  not normal - wilcox
aneuploidy_low_CR <- high_each_aneu_RNA_all_high$MTHFD2_CRISPR[high_each_aneu_RNA_all_high$Aneuploidy_group_classification == "Near-Euploid"]
aneuploidy_med_CR <- high_each_aneu_RNA_all_high$MTHFD2_CRISPR[high_each_aneu_RNA_all_high$Aneuploidy_group_classification == "Intermediate"]
aneuploidy_high_CR <- high_each_aneu_RNA_all_high$MTHFD2_CRISPR[high_each_aneu_RNA_all_high$Aneuploidy_group_classification == "Highly-aneuploid"]

shapiro.test(aneuploidy_low_CR)
shapiro.test(aneuploidy_med_CR)
shapiro.test(aneuploidy_high_CR) # only this normal

my_comp <- list(c("Near-Euploid", "Intermediate"), c("Near-Euploid", "Highly-aneuploid"), 
                c("Intermediate", "Highly-aneuploid"))

ggplot(high_each_aneu_RNA_all_high, aes(x= Aneuploidy_group_classification, y = MTHFD2_CRISPR,
                                        fill = Aneuploidy_group_classification)) +
    geom_boxplot(width = 0.5, fatten = 3) +
    labs(x = "", y = "MTHFD2 CRISPR Essentiality", fill="Aneuploidy group") +
    stat_compare_means(comparisons = my_comp, method = "wilcox") +
    scale_fill_manual(values = c("#FBB13C", "#FE6847", "#B66D0D")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(size = 12), axis.text = element_text(size = 11))
ggsave("plots/CRISPR_MTHFD2high_aneuploidy_group_nice.pdf", device = "pdf", width = 4.2, height = 4)
