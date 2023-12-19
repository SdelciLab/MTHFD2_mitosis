# libraries
library(tidyverse)
library(ggpubr)

# Objective: Calculate correlation between MTHFD2 and spindle assembly checkpoint genes expression

# 1. Get input files -----------------------------------------------------------------------------

# get expression data of MTHFD2 and other mitotic genes
expression_MTHFD2_mitotic <- read.delim("../files/MTHFD2_othergenes_expression_normal_vs_tumor.tsv",
                                        sep = "\t")
head(expression_MTHFD2_mitotic)
str(expression_MTHFD2_mitotic)
colnames(expression_MTHFD2_mitotic)[12:23] <- c("MAD1L1_expression", "MAD2L1_expression",
                                                "BUB1_expression", "BUB1B_expression", 
                                                "BUB3_expression", "CDC20_expression",
                                                "TKK_expression", "CCNB1_expression",
                                                "MD2BP_expression", "AURB_expression",
                                                "ZW10_expression", "ZWILC_expression")

# filter patient data
table(expression_MTHFD2_mitotic$SampleType)

expression_MTHFD2_mitotic_cancer <- expression_MTHFD2_mitotic %>%
  filter(SampleType == "Primary Tumor")
nrow(expression_MTHFD2_mitotic)
nrow(expression_MTHFD2_mitotic_cancer)

# 2. Plot correlation -----------------------------------------------------------------------------

# check coexpression
ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "MAD1L1_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "MAD1L1 expression",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE) + # Add confidence interval  
  stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_MAD1L1.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "MAD2L1_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "MAD2L1 expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
  stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_MAD2L1.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "BUB1_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "BUB1 expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
  stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_BUB1.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, stroke = NA, x = "MTHFD2_expression", y = "BUB1B_expression",
          alpha = 0.2, xlab = "MTHFD2 expression", ylab = "BUB1B expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
  stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_BUB1B.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "BUB3_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "BUB3 expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
  stat_cor(method = "pearson", label.x = 0, label.y = 10)
ggsave("mitotic_genes_coex/coex_MTHFD2_BUB3.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "CDC20_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "CDC20 expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
  stat_cor(method = "pearson", label.x = 0, label.y = 10)
ggsave("mitotic_genes_coex/coex_MTHFD2_CDC20.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "TKK_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "TKK expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
  stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_TKK.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "CCNB1_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "CCNB1 expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
    stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_CCNB1.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "MD2BP_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "MD2BP expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
    stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_MD2BP.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "AURB_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "AURB expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
    stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_AURB.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "ZW10_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "ZW10 expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
    stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_ZW10.pdf", device = "pdf", width = 4, height = 4)

ggscatter(expression_MTHFD2_mitotic, x = "MTHFD2_expression", y = "ZWILC_expression",
          alpha = 0.2, stroke = NA, xlab = "MTHFD2 expression", ylab = "ZWILC expression",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + # Add confidence interval  
    stat_cor(method = "pearson", label.x = 0, label.y = 8)
ggsave("mitotic_genes_coex/coex_MTHFD2_ZWILC.pdf", device = "pdf", width = 4, height = 4)
