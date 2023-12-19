# libraries
library(tidyverse)
library(gtools)
library(ggpubr)

# 1. Read input ------------------------------------------------------------------------------------------------------------------

# obtain MTHFD2 sample sheet
sample_sheet <- read_tsv("../files/sample_sheet_normal_vs_tumor.tsv")
nrow(sample_sheet) # 7,048
MTHFD2_sample_sheet <- read_tsv("../files/MTHFD2_expression_normal_vs_tumor.tsv")
str(MTHFD2_sample_sheet) # 7,048 × 11
colnames(MTHFD2_sample_sheet)

MTHFD2_sample_sheet$SampleType <- factor(MTHFD2_sample_sheet$SampleType, 
                                         levels = c("Solid Tissue Normal", "Primary Tumor"))
table(MTHFD2_sample_sheet$SampleType)

# 2. Statistical analysis of PAIRED DATA  ------------------------------------------------------------------------------------------

# get paired data
MTHFD2_sample_sheet_paired <- MTHFD2_sample_sheet %>%
    filter(paired == "paired") %>%
    mutate(ProjectID_SampleType = paste0(ProjectID, "-", SampleType))
nrow(MTHFD2_sample_sheet_paired) # 1266

# store results from different tests in parameters for paired data

projects <- sort(unique(MTHFD2_sample_sheet_paired$ProjectID))
parameters_paired <- tibble(group = projects)

# split data in projects and sample types to check normality in each group

MTHFD2_paired_per_project_sample <- split(MTHFD2_sample_sheet_paired, 
                                          MTHFD2_sample_sheet_paired$ProjectID_SampleType)
length(MTHFD2_paired_per_project_sample) # 30, 15 project ID: normal and tumor

# check normal distribution with shapiro-wilk normality test

shapiro_results_paired <- map(MTHFD2_paired_per_project_sample, ~ shapiro.test(.x[["MTHFD2_expression"]]))
shapiro_paired_pvalues <- map(shapiro_results_paired, "p.value")
shapiro_paired_pvalues_normal <- unlist(shapiro_paired_pvalues[str_detect(names(shapiro_paired_pvalues), "Solid Tissue Normal")])
shapiro_paired_pvalues_tumor <- unlist(shapiro_paired_pvalues[str_detect(names(shapiro_paired_pvalues), "Primary Tumor")])

parameters_paired$shapiro.pval_normal <- shapiro_paired_pvalues_normal
parameters_paired$shapiro.pval_tumor <- shapiro_paired_pvalues_tumor
parameters_paired <- parameters_paired %>%
    mutate(normality = ifelse(shapiro.pval_normal > 0.05 & shapiro.pval_tumor > 0.05, "normal", "not normal"))

# TCGA-BLCA both normal
# TCGA-BRCA healthy not normal -> non-parametric
# TCGA-COAD both normal
# TCGA-ESCA both normal
# TCGA-HNSC tumor not normal -> non parametric
# TCGA-KICH both normal
# TCGA-KIRC healthy not normal -> non-parametric
# TCGA-KIRP tumor not normal -> non parametric
# TCGA-LIHC both normal
# TCGA-LUAD both normal
# TCGA-LUSC both normal
# TCGA-PRAD tumor not normal -> non parametric
# TCGA-STAD both normal
# TCGA-THCA both normal
# TCGA-UCEC both normal

# split normal data in projects to check variance in each group

paired_normal_groups_name <- parameters_paired %>%
    filter(normality == "normal") %>%
    select(group) %>%
    unlist()

MTHFD2_paired_per_project <- MTHFD2_sample_sheet_paired %>%
    filter(ProjectID %in% paired_normal_groups_name) %>%
    split(.$ProjectID)
length(MTHFD2_paired_per_project) # 10, 15 - 5 not normal

# check if groups that are normal have the same variances or not with var.test

var_test_results_paired <- map(MTHFD2_paired_per_project,
                               ~ var.test(MTHFD2_expression ~ SampleType, data = .x))
var_test_results_paired
# TCGA-BLCA equal.var = T
# TCGA-COAD equal.var = T
# TCGA-ESCA equal.var = T
# TCGA-KICH equal.var = F
# TCGA-LIHC equal.var = T
# TCGA-LUAD bequal.var = F
# TCGA-LUSC equal.var = F
# TCGA-STAD equal.var = T
# TCGA-THCA equal.var = T
# TCGA-UCEC equal.var = T

var_test_paired_pvalues <- map(var_test_results_paired, "p.value")
var_test_paired_pvalues <- as_tibble(var_test_paired_pvalues)
var_test_paired_pvalues <- var_test_paired_pvalues %>%
    pivot_longer(cols = 1:ncol(var_test_paired_pvalues), names_to = "group", values_to = "var.test.pval")

parameters_paired <- parameters_paired %>%
    left_join(var_test_paired_pvalues, by = "group")

# 3. Wilcox.test for all ---------------------------------------------------------------------------------------------------------

# split MTHFD2 sample sheet paired in projects 
wilcox_test_all_per_project <- MTHFD2_sample_sheet_paired %>%
    split(.$ProjectID)
length(wilcox_test_all_per_project) # 15

# perform wilcox test
wilcox_test_paired_all_results <- map(wilcox_test_all_per_project, 
                                      ~ wilcox.test(MTHFD2_expression ~ SampleType, data = .x,
                                                                                  paired = TRUE))
wilcox_test_paired_all_results
# TCGA-BLCA pval < 0.05
# TCGA-BRCA pval < 0.05
# TCGA-COAD pval < 0.05
# TCGA-ESCA pval < 0.05
# TCGA-HNSC pval < 0.05
# TCGA-KICH pval > 0.05 -> Not significative
# TCGA-KIRC pval < 0.05
# TCGA-KIRP pval < 0.05
# TCGA-LIHC pval < 0.05
# TCGA-LUAD pval < 0.05
# TCGA-LUSC pval < 0.05
# TCGA-PRAD pval < 0.05
# TCGA-STAD pval < 0.05
# TCGA-THCA pval > 0.05 -> Not significative
# TCGA-UCEC pval < 0.05

# add results to parameters
all_wilcox_test_pvalues <- map(wilcox_test_paired_all_results, "p.value")
all_wilcox_test_pvalues <- as.tibble(all_wilcox_test_pvalues)
all_wilcox_test_pvalues <- all_wilcox_test_pvalues %>%
    pivot_longer(cols = 1:ncol(all_wilcox_test_pvalues), names_to = "group", values_to = "all.wilcox.test.pval")

parameters_paired <- parameters_paired %>%
    left_join(all_wilcox_test_pvalues, by = "group")

# use compare means function to check
wilcox_paired_all_compare_means <- map(wilcox_test_all_per_project, 
                                       ~ compare_means(MTHFD2_expression ~ SampleType, data = .x,
                                                       method = "wilcox.test", paired = TRUE))
wilcox_paired_all_compare_means # same values

# create df with all pvalues
all_wilcox_pvalues_df <- tibble(ProjectID = projects)
all_wilcox_pvalues_df$group1 <- c(rep("Solid Tissue Normal", nrow(all_wilcox_pvalues_df)))
all_wilcox_pvalues_df$group2 <- c(rep("Primary Tumor", nrow(all_wilcox_pvalues_df)))
all_wilcox_pvalues_df$test <- c(rep("wilcox.test", nrow(all_wilcox_pvalues_df)))
all_wilcox_pvalues_df$pval <- parameters_paired$all.wilcox.test.pval
all_wilcox_pvalues_df$pval <- formatC(all_wilcox_pvalues_df$pval, format = "e", digits = 2)
all_wilcox_pvalues_df$sig <- stars.pval(as.double(all_wilcox_pvalues_df$pval))
all_wilcox_pvalues_df$pval_sig <- paste0(all_wilcox_pvalues_df$pval, all_wilcox_pvalues_df$sig)

# plot for publication
MTHFD2_sample_sheet_paired %>%
    ggplot(aes(SampleType, MTHFD2_expression)) +
    geom_violin(aes(fill = SampleType), col = NA, show.legend = T) +
    geom_boxplot(width = 0.3, alpha = 0, fatten = 3) +
    stat_pvalue_manual(all_wilcox_pvalues_df, label = "pval_sig", y.position = 9.1, xmax = "group2") +
    facet_grid(~ ProjectID) +
    labs(subtitle = "Paired Data\n", x = "", y = "\nMTHFD2 expression\nlog2(TPM+0.01)\n", fill = "") +
    ggtitle("\nTCGA Gene Expression Data - RNA-seq\n") +
    scale_fill_manual(labels = c("Normal", "Tumor"), values = c("#C0C0C0", "#9588B6")) +
    ylim(-3, 10) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 15), axis.text = element_text(size = 13), axis.text.x = element_blank(), 
          legend.text = element_text(size = 16), axis.ticks = element_blank(), 
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          title = element_text(size = 16), plot.title = element_text(hjust = 0.5),  
          plot.subtitle = element_text(size = 16, hjust = 0.5),
          panel.grid.major.x = element_blank())
ggsave("../plots/Violin_MTHFD2exp_paired_pvals_wilcox_nicer.pdf", device = "pdf", 
       width = 24, height = 5)

# save parameters paired results to output
write_tsv(parameters_paired, "MTHFD2_normal_vs_tumor_parameters_paired_results.tsv")

# plot with all tissues together 
MTHFD2_sample_sheet_paired$all_data <- rep("All_cancer", nrow(MTHFD2_sample_sheet_paired))

# pval all together
wilcox.test(MTHFD2_expression ~ SampleType, data = MTHFD2_sample_sheet_paired, paired = TRUE)
# p-value < 2.2e-16

# plot all tissues together
MTHFD2_sample_sheet_paired %>%
    ggplot(aes(SampleType, MTHFD2_expression)) +
    geom_violin(aes(fill = SampleType), col = NA, show.legend = T) +
    geom_boxplot(width = 0.3, alpha = 0, fatten = 3) +
    stat_compare_means(paired = TRUE, method = "wilcox") +
    facet_grid(~ all_data) +
    labs(subtitle = "Paired Data\n", x = "", y = "\nMTHFD2 expression\nlog2(TPM+0.01)\n", fill = "") +
    ggtitle("\nTCGA Gene Expression Data - RNA-seq\n") +
    scale_fill_manual(labels = c("Normal", "Tumor"), values = c("#C0C0C0", "#9588B6")) +
    ylim(-3, 10) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 15), axis.text = element_text(size = 13), axis.text.x = element_blank(), 
          legend.text = element_text(size = 16), axis.ticks = element_blank(), 
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          title = element_text(size = 16), plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(size = 16, hjust = 0.5),
          panel.grid.major.x = element_blank())
ggsave("../plots/Violin_MTHFD2exp_paired_pvals_wilcox_all_nicer.pdf", device = "pdf", 
       width = 4, height = 5)

