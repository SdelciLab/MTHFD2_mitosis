# libraries
library(tidyverse)
library(ggpubr)

# get file
MTH_quant <- read.delim("MTHFD2rb hct mcf7 h358 40X 181022__2022-10-18T17_44_44-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected Single.txt", skip = 9)
head(MTH_quant)
colnames(MTH_quant)

# add conditions
MTH_quant_cond <- MTH_quant %>%
    mutate(Condition = case_when(
        Row == 2 ~ "HCT116",
        Row == 3 ~ "KO1",
        Row == 4 ~ "KO2"))

# select only cell signalling antibody
MTH_quant_cond_2nd <- MTH_quant_cond %>%
    filter(Column == 4 | Column == 5)

# add logs
MTH_quant_cell_log <- MTH_quant_cond_2nd %>%
    mutate(logIntensityCyt = log2(All.Nuclei.Selected.Single...Intensity.Cytoplasm.Alexa.488.Mean),
           logIntensityNuc = log2(All.Nuclei.Selected.Single...Intensity.Nucleus.Alexa.488.Mean),
           ratioMI = logIntensityNuc/logIntensityCyt,
           logIntDenCyt =  log2(All.Nuclei.Selected.Single...Intensity.Cytoplasm.Alexa.488.Mean * All.Nuclei.Selected.Single...Area.Cytoplasm.Area..µm..),
           logIntDenNuc =  log2(All.Nuclei.Selected.Single...Intensity.Nucleus.Alexa.488.Mean * All.Nuclei.Selected.Single...Area.Nuclei.Area..µm..),
           ratioID = logIntDenNuc/logIntDenCyt)

# factor
MTH_quant_cell_log$Condition <- factor(MTH_quant_cell_log$Condition, levels = c("HCT116", "KO1", "KO2"))

# comparisons
my_comp <- list(c("HCT116", "KO2"), c("HCT116", "KO1"))

# remember change KO1 and KO2 name

# plot cytosolic signal
ggplot(MTH_quant_cell_log, aes(x=Condition, y=logIntensityCyt, fill = Condition)) +
    geom_boxplot() +
    labs(x = "", y = "log2(MTHFD2 cytosolic mean intensity)", fill = "") +
    scale_fill_manual(values=c("HCT116" = "#aaafb0", "KO1" = "#e56e58", "KO2" = "#f9beb3")) +
    stat_compare_means(test = "wilcox", comparisons = my_comp, size = 3) +
    theme_classic() +
    theme(legend.position = "none")
ggsave("plots/logInt_cyt_HCT_KO.pdf", device = "pdf", width = 3, height = 3)

# plot nuclear signal
ggplot(MTH_quant_cell_log, aes(x=Condition, y=logIntensityNuc, fill = Condition)) +
    geom_boxplot() +
    labs(x = "", y = "log2(MTHFD2 nuclear mean intensity)", fill = "") +
    stat_compare_means(test = "wilcox", comparisons = my_comp, size = 3) +
    scale_fill_manual(values=c("HCT116" = "#aaafb0", "KO1" = "#e56e58", "KO2" = "#f9beb3")) +
    theme_classic() +
    theme(legend.position = "none")
ggsave("plots/logInt_nuc_HCT_KO.pdf", device = "pdf", width = 3, height = 3)

# numbers
table(MTH_quant_cell_log_WT_KO$Condition)