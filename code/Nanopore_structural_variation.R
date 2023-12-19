# libraries
library(tidyverse)
library(nortest)
library(gtools)
library(ggpubr)

# 1. Get input data of all structural variation  ---------------------------------------------------------------------------------

# get data
del_WT <- read.delim("del_plot.ok.bg.plus", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))
del_KO <- read.delim("del_plot.ok.bg.minus", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))

ins_WT <- read.delim("ins_plot.ok.bg.plus", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))
ins_KO <- read.delim("ins_plot.ok.bg.minus", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))

dup_WT <- read.delim("dup_plot.ok.bg.plus", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))
dup_KO <- read.delim("dup_plot.ok.bg.minus", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))

inv_WT <- read.delim("inv_plot.ok.bg.plus", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))
inv_KO <- read.delim("inv_plot.ok.bg.minus", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))

# add info to combine data tables
del_WT_add <- del_WT %>%
    mutate(Condition = rep("WT", nrow(del_WT)),
           Alteration = rep("Deletion", nrow(del_WT)))
del_KO_add <- del_KO %>%
    mutate(Condition = rep("KO", nrow(del_KO)),
           Alteration = rep("Deletion", nrow(del_KO)))

ins_WT_add <- ins_WT %>%
    mutate(Condition = rep("WT", nrow(ins_WT)),
           Alteration = rep("Insertion", nrow(ins_WT)))
ins_KO_add <- ins_KO %>%
    mutate(Condition = rep("KO", nrow(ins_KO)),
           Alteration = rep("Insertion", nrow(ins_KO)))

dup_WT_add <- dup_WT %>%
    mutate(Condition = rep("WT", nrow(dup_WT)),
           Alteration = rep("Duplication", nrow(dup_WT)))
dup_KO_add <- dup_KO %>%
    mutate(Condition = rep("KO", nrow(dup_KO)),
           Alteration = rep("Duplication", nrow(dup_KO)))

inv_WT_add <- inv_WT %>%
    mutate(Condition = rep("WT", nrow(inv_WT)),
           Alteration = rep("Inversion", nrow(inv_WT)))
inv_KO_add <- inv_KO %>%
    mutate(Condition = rep("KO", nrow(inv_KO)),
           Alteration = rep("Inversion", nrow(inv_KO)))

# combine data tables
table(colnames(del_WT_add) == colnames(del_KO_add))
table(colnames(del_WT_add) == colnames(ins_WT_add))
table(colnames(del_WT_add) == colnames(ins_KO_add))
table(colnames(del_WT_add) == colnames(dup_WT_add))
table(colnames(del_WT_add) == colnames(dup_KO_add))
table(colnames(del_WT_add) == colnames(inv_WT_add))
table(colnames(del_WT_add) == colnames(inv_KO_add))

SV_all <- rbind(del_WT_add, del_KO_add, ins_WT_add, ins_KO_add,
                dup_WT_add, dup_KO_add, inv_WT_add, inv_KO_add)
table(SV_all$Condition, SV_all$Alteration)

# 2. Analyse and plot ---------------------------------------------------------------------------------------------

# plot n
SV_all_summary <- SV_all %>%
    group_by(Alteration, Condition) %>%
    summarise(number = n())

SV_all_summary$Alteration <- factor(SV_all_summary$Alteration,
                                    levels = c("Insertion", "Deletion",
                                               "Duplication", "Inversion"))

SV_all_summary$Condition <- factor(SV_all_summary$Condition,
                                   levels = c("WT", "KO"))

ggplot(SV_all_summary, aes(x = Condition, y = number, fill = Condition)) +
    geom_col(width = 0.5) +
    facet_wrap(~Alteration, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO" = "#e56e58")) +
    labs(y = "Number of alterations") +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 11), legend.position = "none",
          strip.text = element_text(size = 11), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/SV_specific_number.pdf", device = "pdf", width = 5.5, height = 2)

# plot sizes
SV_all_sizes <- SV_all %>%
    mutate(size_abs = abs(size))

SV_all_sizes$Alteration <- factor(SV_all_sizes$Alteration,
                                  levels = c("Insertion", "Deletion",
                                             "Duplication", "Inversion"))

SV_all_sizes$Condition <- factor(SV_all_sizes$Condition,
                                 levels = c("WT", "KO"))

# summarise
SV_all_sizes_summary <- SV_all_sizes %>%
    group_by(Alteration, Condition) %>%
    summarise(mean_size = mean(size_abs),
              median_size = median(size_abs),
              sd_size = sd(size_abs),
              se_size = sd_size/sqrt(n()))

SV_all_sizes_summary$Alteration <- factor(SV_all_sizes_summary$Alteration,
                                          levels = c("Insertion", "Deletion",
                                                     "Duplication", "Inversion"))

SV_all_sizes_summary$Condition <- factor(SV_all_sizes_summary$Condition,
                                         levels = c("WT", "KO"))

# # check normality
ins_WT_size <- SV_all_sizes[SV_all_sizes$Condition == "WT" 
                            & SV_all_sizes$Alteration == "Insertion", "size_abs"]
ins_KO_size <- SV_all_sizes[SV_all_sizes$Condition == "KO" 
                            & SV_all_sizes$Alteration == "Insertion", "size_abs"]

lillie.test(ins_WT_size)
lillie.test(ins_KO_size)

del_WT_size <- SV_all_sizes[SV_all_sizes$Condition == "WT" 
                            & SV_all_sizes$Alteration == "Deletion", "size_abs"]
del_KO_size <- SV_all_sizes[SV_all_sizes$Condition == "KO" 
                            & SV_all_sizes$Alteration == "Deletion", "size_abs"]

lillie.test(del_WT_size)
lillie.test(del_KO_size)

dup_WT_size <- SV_all_sizes[SV_all_sizes$Condition == "WT" 
                            & SV_all_sizes$Alteration == "Duplication", "size_abs"]
dup_KO_size <- SV_all_sizes[SV_all_sizes$Condition == "KO" 
                            & SV_all_sizes$Alteration == "Duplication", "size_abs"]

lillie.test(dup_WT_size)
lillie.test(dup_KO_size)

inv_WT_size <- SV_all_sizes[SV_all_sizes$Condition == "WT" 
                            & SV_all_sizes$Alteration == "Inversion", "size_abs"]
inv_KO_size <- SV_all_sizes[SV_all_sizes$Condition == "KO" 
                            & SV_all_sizes$Alteration == "Inversion", "size_abs"]

shapiro.test(inv_WT_size) # only this normal
shapiro.test(inv_KO_size)

# use wilcox test
wilcox_ins <- wilcox.test(ins_WT_size, ins_KO_size) # pval 0.0003198
wilcox_del <- wilcox.test(del_WT_size, del_KO_size) # pval 0.000526
wilcox_dup <- wilcox.test(dup_WT_size, dup_KO_size) # pval 0.9244 no
wilcox_inv <- wilcox.test(inv_WT_size, inv_KO_size) # pval 0.6646 no

all_wilcox_pvalues_df <- tibble(Alteration = c("Insertion", "Deletion", "Duplication", "Inversion"))
all_wilcox_pvalues_df$Alteration <- factor(all_wilcox_pvalues_df$Alteration,
                                           levels = c("Insertion", "Deletion",
                                                      "Duplication", "Inversion"))
all_wilcox_pvalues_df$group1 <- c(rep("WT", nrow(all_wilcox_pvalues_df)))
all_wilcox_pvalues_df$group2 <- c(rep("KO", nrow(all_wilcox_pvalues_df)))
all_wilcox_pvalues_df$test <- c(rep("wilcox.test", nrow(all_wilcox_pvalues_df)))
all_wilcox_pvalues_df$pval <- c(wilcox_ins$p.value, wilcox_del$p.value, 
                                wilcox_dup$p.value, wilcox_inv$p.value)
all_wilcox_pvalues_df$pval <- round(all_wilcox_pvalues_df$pval, 4)
all_wilcox_pvalues_df$sig <- stars.pval(as.double(all_wilcox_pvalues_df$pval))
all_wilcox_pvalues_df$pval_sig <- paste0(all_wilcox_pvalues_df$pval, all_wilcox_pvalues_df$sig)
all_wilcox_pvalues_df$y.position <- c(11500, 85000, 105000, 95000)

# plot all
ggplot() +
    geom_jitter(data = SV_all_sizes, aes(x = Condition, y = size_abs, color = Condition), 
                alpha = 0.7, stroke = NA) +
    geom_point(data = SV_all_sizes_summary, aes(x = Condition, y = mean_size),
               color = "black") +
    geom_errorbar(data = SV_all_sizes_summary, aes(x = Condition, ymin=mean_size-sd_size,
                                                   ymax=mean_size+sd_size),
                  color = "black", width=.25) +
    facet_wrap(~Alteration, scales = "free_y", nrow = 1) +
    stat_pvalue_manual(all_wilcox_pvalues_df, label = "pval_sig", xmax = "group2") +
    labs(y = "Size of alteration (bp)") +
    scale_color_manual(values = c("WT" = "#aaafb0", "KO" = "#e56e58")) +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 11), legend.position = "none",
          strip.text = element_text(size = 11), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/alterations_sizes_all.pdf", device = "pdf",  width = 7, height = 2.8)

# 3. Get input data of centromeric structural variation  ---------------------------------------------------------------------------------

# get data
del_WT_centr <- read.delim("del_plot.ok.bg.plus.centro", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))
del_KO_centr <- read.delim("del_plot.ok.bg.minus.centro", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))

ins_WT_centr <- read.delim("ins_plot.ok.bg.plus.centro", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))
ins_KO_centr <- read.delim("ins_plot.ok.bg.minus.centro", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))

dup_WT_centr <- read.delim("dup_plot.ok.bg.plus.centro", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))
dup_KO_centr <- read.delim("dup_plot.ok.bg.minus.centro", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))

inv_WT_centr <- read.delim("inv_plot.ok.bg.plus.centro", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))
inv_KO_centr <- read.delim("inv_plot.ok.bg.minus.centro", sep = "\t", col.names = c("chr", "init", "end", "size", "other"))

# add info to combine data tables
del_WT_add <- del_WT_centr %>%
    mutate(Condition = rep("WT", nrow(del_WT_centr)),
           Alteration = rep("Deletion", nrow(del_WT_centr)))
del_KO_add <- del_KO_centr %>%
    mutate(Condition = rep("KO", nrow(del_KO_centr)),
           Alteration = rep("Deletion", nrow(del_KO_centr)))

ins_WT_add <- ins_WT_centr %>%
    mutate(Condition = rep("WT", nrow(ins_WT_centr)),
           Alteration = rep("Insertion", nrow(ins_WT_centr)))
ins_KO_add <- ins_KO_centr %>%
    mutate(Condition = rep("KO", nrow(ins_KO_centr)),
           Alteration = rep("Insertion", nrow(ins_KO_centr)))

dup_WT_add <- dup_WT_centr %>%
    mutate(Condition = rep("WT", nrow(dup_WT_centr)),
           Alteration = rep("Duplication", nrow(dup_WT_centr)))
dup_KO_add <- dup_KO_centr %>%
    mutate(Condition = rep("KO", nrow(dup_KO_centr)),
           Alteration = rep("Duplication", nrow(dup_KO_centr)))

inv_WT_add <- inv_WT_centr %>%
    mutate(Condition = rep("WT", nrow(inv_WT_centr)),
           Alteration = rep("Inversion", nrow(inv_WT_centr)))
inv_KO_add <- inv_KO_centr %>%
    mutate(Condition = rep("KO", nrow(inv_KO_centr)),
           Alteration = rep("Inversion", nrow(inv_KO_centr)))

# combine data tables
table(colnames(del_WT_add) == colnames(del_KO_add))
table(colnames(del_WT_add) == colnames(ins_WT_add))
table(colnames(del_WT_add) == colnames(ins_KO_add))
table(colnames(del_WT_add) == colnames(dup_WT_add))
table(colnames(del_WT_add) == colnames(dup_KO_add))
table(colnames(del_WT_add) == colnames(inv_WT_add))
table(colnames(del_WT_add) == colnames(inv_KO_add))

SV_all_centr <- rbind(del_WT_add, del_KO_add, ins_WT_add, ins_KO_add,
                      dup_WT_add, dup_KO_add, inv_WT_add, inv_KO_add)
table(SV_all_centr$Condition, SV_all_centr$Alteration)

# 4. Analyse and plot ---------------------------------------------------------------------------------------------

# plot n
SV_all_centr_summary <- SV_all_centr %>%
    group_by(Alteration, Condition) %>%
    summarise(number = n())

SV_all_centr_summary$Alteration <- factor(SV_all_centr_summary$Alteration,
                                          levels = c("Insertion", "Deletion",
                                                     "Duplication", "Inversion"))

SV_all_centr_summary$Condition <- factor(SV_all_centr_summary$Condition,
                                         levels = c("WT", "KO"))

ggplot(SV_all_centr_summary, aes(x = Condition, y = number, fill = Condition)) +
    geom_col(width = 0.5) +
    facet_wrap(~Alteration, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO" = "#e56e58")) +
    labs(y = "Number of alterations") +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 11), legend.position = "none",
          strip.text = element_text(size = 11), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/SV_specific_number_centro.pdf", device = "pdf", width = 5.5, height = 2)
