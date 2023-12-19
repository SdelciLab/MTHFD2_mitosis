# libraries
library(tidyverse)
library(car)
library(ggpubr)

# 1. Get input data ----------------------------------------------------------------------------------

# get files
centromere_quant <- read_delim("summary_counts.csv", delim = ",")

# each file corresponds to one condition, one time point, one replicate
centromere_quant <- centromere_quant %>%
    mutate(condition = case_when(
        str_detect(sample, "WT") ~ "WT",
        str_detect(sample, "K1") ~ "KO1",
        time = case_when(
            str_detect(sample, "T0") ~ "T0",
            str_detect(sample, "T7") ~ "T7"),
        replicate = case_when(
            str_detect(sample, "R1") ~ "1",
            str_detect(sample, "R2") ~ "2",
            str_detect(sample, "R3") ~ "3"))) %>%
            mutate(condition_time = paste(condition, time, sep = "_"))
        centromere_quant$condition_time <- factor(centromere_quant$condition_time, levels = c("WT_T0", "WT_T7", "KO1_T0", "KO1_T7"))
        centromere_quant$condition <- factor(centromere_quant$condition, levels = c("WT", "KO1"))
        

# 2. Analyse and plot centromeric intersect --------------------------------------------------------------

# analyse percentage of intersects in all chromosomes
centromere_quant_all <- centromere_quant %>%
    group_by(condition, time, replicate) %>%
    summarise(per_sum = sum(per), count_sum = sum(counts)) %>%
    mutate(condition_time = paste(condition, time, sep = "_"))

centromere_quant_all$condition_time <- factor(centromere_quant_all$condition_time, levels = c("WT_T0", "WT_T7", "KO1_T0", "KO1_T7"))
centromere_quant_all$condition <- factor(centromere_quant_all$condition, levels = c("WT", "KO1"))

# get summaries of condition and time
centromere_quant_all_cond_time_global_stats <- centromere_quant_all %>%
    group_by(condition_time) %>%
    summarise(n = n(), count_mean = mean(count_sum), perc_mean = mean(per_sum), 
              perc_sd = sd(per_sum), perc_se = sd(per_sum)/sqrt(3))
centromere_quant_all_cond_time_global_stats$condition_time <- factor(centromere_quant_all_cond_time_global_stats$condition_time, 
                                                                     levels = c("WT_T0", "WT_T7", "KO1_T0", "KO1_T7"))
write_csv(centromere_quant_all_cond_time_global_stats, "all_intersect_cond_time_global_stats.csv")

# get summaries of condition 
centromere_quant_all_cond_global_stats <- centromere_quant_all %>%
    group_by(condition) %>%
    summarise(n = n(), count_mean = mean(count_sum), perc_mean = mean(per_sum), 
              perc_sd = sd(per_sum), perc_se = sd(per_sum)/sqrt(6))
centromere_quant_all_cond_global_stats$condition <- factor(centromere_quant_all_cond_global_stats$condition, 
                                                           levels = c("WT", "KO1"))
write_csv(centromere_quant_all_cond_global_stats, "all_intersect_cond_global_stats.csv")

# plot
wilcox_pvals_KO <- compare_means(per_sum ~ condition, data = centromere_quant_all, 
                                 ref.group = "WT", method = "wilcox.test") %>%
    mutate(y.position = c(12.5e-05))

set.seed(123)
ggplot() +
    geom_boxplot(data = centromere_quant_all, aes(x = condition, y = per_sum, col = condition), width = 0.4, lwd=0.5) +
    geom_jitter(data = centromere_quant_all, aes(x = condition, y = per_sum, col = condition), size = 1.25) +
    scale_color_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58")) +
    xlab("") +
    ylab("\n% of intersects / total mapped reads") + 
    ylim(c(0.000075,0.00013)) +
    stat_pvalue_manual(wilcox_pvals_KO, label = "p.signif", size = 5, 
                       coord.flip = FALSE, tip.length = 0.01,  bracket.size = 0.2, vjust = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5, size = 11),
          axis.text.y = element_text(vjust = 1, hjust=1, size = 8), 
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), size = 11),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position="none") 
ggsave("plots/all_centrom_quant_cond_point_pval_KO.pdf", width = 3.4, height = 3)