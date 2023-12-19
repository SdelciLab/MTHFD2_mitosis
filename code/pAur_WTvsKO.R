# libraries
library(tidyverse)
library(car)
library(ggpubr)
library(rstatix)

# get data
pAur_rep1 <- read.delim("WTvsKO-crest-ABC_2_rep1__2023-12-13T10_38_26-Measurement 1/Evaluation1/Objects_Population - Nuclei Selected.txt", skip = 9)
head(pAur_rep1)

pAur_rep1_2 <- read.delim("WTvsKO-crest-ABC_3_rep1__2023-12-13T12_37_09-Measurement 1/Evaluation1/Objects_Population - Nuclei Selected.txt", skip = 9)
head(pAur_rep1_2)

pAur_rep2 <- read.delim("WTvsKO-crest-ABC_2_rep2__2023-12-13T10_25_18-Measurement 1/Evaluation1/Objects_Population - Nuclei Selected.txt", skip = 9)
head(pAur_rep2)

pAur_rep2_2 <- read.delim("WTvsKO-crest-ABC_3_rep2__2023-12-13T11_46_31-Measurement 1/Evaluation1/Objects_Population - Nuclei Selected.txt", skip = 9)
head(pAur_rep2_2)

pAur_rep3 <- read.delim("WTvsKO-crest-ABC_2_rep3__2023-12-13T10_51_19-Measurement 1/Evaluation1/Objects_Population - Nuclei Selected.txt", skip = 9)
head(pAur_rep3)

pAur_rep3_2 <- read.delim("WTvsKO-crest-ABC_3_rep3__2023-12-13T11_06_41-Measurement 1/Evaluation1/Objects_Population - Nuclei Selected.txt", skip = 9)
head(pAur_rep3_2)

# check colnames are the same
table(colnames(pAur_rep1) == colnames(pAur_rep1_2))
table(colnames(pAur_rep1) == colnames(pAur_rep2))
table(colnames(pAur_rep1) == colnames(pAur_rep2_2))
table(colnames(pAur_rep1) == colnames(pAur_rep3))
table(colnames(pAur_rep1) == colnames(pAur_rep3_2))

# check
colnames(pAur_rep1)

# distribution
table(pAur_rep1$Row,  pAur_rep1$Column) 
table(pAur_rep1_2$Row,  pAur_rep1_2$Column) 
table(pAur_rep2$Row,  pAur_rep2$Column) 
table(pAur_rep2_2$Row,  pAur_rep2_2$Column) 
table(pAur_rep3$Row,  pAur_rep3$Column) 
table(pAur_rep3_2$Row,  pAur_rep3_2$Column) 

# add replicate and plate
pAur_rep1 <- pAur_rep1 %>%
    mutate(Replicate = rep(1, nrow(pAur_rep1)),
           Plate = rep(1, nrow(pAur_rep1)))
pAur_rep1_2 <- pAur_rep1_2 %>%
    mutate(Replicate = rep(1, nrow(pAur_rep1_2)),
           Plate = rep(2, nrow(pAur_rep1_2)))
pAur_rep2 <- pAur_rep2 %>%
    mutate(Replicate = rep(2, nrow(pAur_rep2)),
           Plate = rep(1, nrow(pAur_rep2)))
pAur_rep2_2 <- pAur_rep2_2 %>%
    mutate(Replicate = rep(2, nrow(pAur_rep2_2)),
           Plate = rep(2, nrow(pAur_rep2_2)))
pAur_rep3 <- pAur_rep3 %>%
    mutate(Replicate = rep(3, nrow(pAur_rep3)),
           Plate = rep(1, nrow(pAur_rep3)))
pAur_rep3_2 <- pAur_rep3_2 %>%
    mutate(Replicate = rep(3, nrow(pAur_rep3_2)),
           Plate = rep(2, nrow(pAur_rep3_2)))

# put all data together -> since colnames are the same, do rbind
pAur_alldata <- rbind(pAur_rep1, pAur_rep1_2, pAur_rep2, pAur_rep2_2, 
                      pAur_rep3, pAur_rep3_2)

# check
table(pAur_alldata$Replicate, pAur_alldata$Row) 
table(pAur_alldata$Replicate, pAur_alldata$Column) 
table(pAur_alldata$Plate, pAur_alldata$Row) 
table(pAur_alldata$Plate, pAur_alldata$Column) 

# add condition
pAur_alldata_cond <- pAur_alldata %>%
    mutate(Condition = case_when(
        Row == 2 ~ "WT",
        Row == 3 ~ "KO"
    ))

# filter single nuclei
pAur_alldata_cond_singl <- pAur_alldata_cond %>%
    filter(`Nuclei.Selected...Nucleus.Morphological.Properties.Area..µm..` > 60 &
               `Nuclei.Selected...Nucleus.Morphological.Properties.Area..µm..` < 300 &
               Nuclei.Selected...Nucleus.Morphological.Properties.Roundness > 0.7)

nrow(pAur_alldata_cond) # 141280
nrow(pAur_alldata_cond_singl) # 133882
table(pAur_alldata_cond_singl$Condition)

# get logs
pAur_alldata_cond_singl_log <- pAur_alldata_cond_singl %>%
    mutate(log_488_nuclear = log2(Nuclei.Selected...Intensity.Nucleus.Alexa.488.Mean),
           log_488_cyt = log2(Nuclei.Selected...Intensity.Cytoplasm.Alexa.488.Mean))

# separate WT and KO
pAur_alldata_cond_singl_WT <- pAur_alldata_cond_singl_log %>%
    filter(Condition == "WT")
pAur_alldata_cond_singl_KO <- pAur_alldata_cond_singl_log %>%
    filter(Condition == "KO")

# check aur ABC nucleus in histogram
ggplot(pAur_alldata_cond_singl_WT, aes(x=log_488_nuclear)) +
    geom_histogram(bins=100) +
    facet_grid(Replicate ~ Plate) +
    xlab("log2(pAurA/B/C mean intensity)") +
    theme_bw()
ggsave("plots/histogram_pAurABC_nucleus_WT.pdf")

ggplot(pAur_alldata_cond_singl_KO, aes(x=log_488_nuclear)) +
    geom_histogram(bins=100) +
    facet_grid(Replicate ~ Plate) +
    xlab("log2(pAurA/B/C mean intensity)") +
    theme_bw()
ggsave("plots/histogram_pAurABC_nucleus_KO.pdf")

# check aur ABC cytosol in histogram
ggplot(pAur_alldata_cond_singl_WT, aes(x=log_488_cyt)) +
    geom_histogram(bins=100) +
    facet_grid(Replicate ~ Plate) +
    xlab("log2(pAurA/B/C mean intensity)") +
    theme_bw()
ggsave("plots/histogram_pAurABC_cyt_WT.pdf")

ggplot(pAur_alldata_cond_singl_KO, aes(x=log_488_cyt)) +
    geom_histogram(bins=100) +
    facet_grid(Replicate ~ Plate) +
    xlab("log2(pAurA/B/C mean intensity)") +
    theme_bw()
ggsave("plots/histogram_pAurABC_cyt_KO.pdf")

# get higher values
pAur_alldata_cond_singl_log_high <- pAur_alldata_cond_singl_log %>%
    filter(log_488_nuclear > 13.75)

# histogram higher values
pAur_alldata_cond_singl_high_WT <- pAur_alldata_cond_singl_log_high %>%
    filter(Condition == "WT")
pAur_alldata_cond_singl_high_KO <- pAur_alldata_cond_singl_log_high %>%
    filter(Condition == "KO")

# check aur ABC nucleus in histogram
ggplot(pAur_alldata_cond_singl_high_WT, aes(x=log_488_nuclear)) +
    geom_histogram(bins=100) +
    facet_grid(Replicate ~ Plate) +
    xlab("log2(pAurA/B/C mean intensity)") +
    theme_bw()
ggsave("plots/histogram_pAurABC_nucleus_WT_high.pdf")

ggplot(pAur_alldata_cond_singl_high_KO, aes(x=log_488_nuclear)) +
    geom_histogram(bins=100) +
    facet_grid(Replicate ~ Plate) +
    xlab("log2(pAurA/B/C mean intensity)") +
    theme_bw()
ggsave("plots/histogram_pAurABC_nucleus_KO_high.pdf")

# check aur ABC cytosol in histogram
ggplot(pAur_alldata_cond_singl_high_WT, aes(x=log_488_cyt)) +
    geom_histogram(bins=100) +
    facet_grid(Replicate ~ Plate) +
    xlab("log2(pAurA/B/C mean intensity)") +
    theme_bw()
ggsave("plots/histogram_pAurABC_cyt_WT_high.pdf")

ggplot(pAur_alldata_cond_singl_high_KO, aes(x=log_488_cyt)) +
    geom_histogram(bins=100) +
    facet_grid(Replicate ~ Plate) +
    xlab("log2(pAurA/B/C mean intensity)") +
    theme_bw()
ggsave("plots/histogram_pAurABC_cyt_KO_high.pdf")

# get background
pAur_alldata_cond_singl_background <- pAur_alldata_cond_singl_log %>%
    group_by(Condition) %>%
    summarise(n= n(),
              mean_aur_nucleus = mean(log_488_nuclear),
              median_aur_nucleus = median(log_488_nuclear),
              mean_aur_cyt = mean(log_488_cyt, na.rm = TRUE),
              median_aur_cyt = median(log_488_cyt, na.rm = TRUE))
pAur_alldata_cond_singl_background

# background
background_WT <- pAur_alldata_cond_singl_background %>%
    filter(Condition == "WT") %>%
    pull(mean_aur_nucleus)
background_KO <- pAur_alldata_cond_singl_background %>%
    filter(Condition == "KO") %>%
    pull(mean_aur_nucleus)

# normalise by mean all cells (background)
pAur_alldata_cond_singl_log_norm <- pAur_alldata_cond_singl_log %>%
    mutate(nuclear_488_norm = case_when(
        Condition == "WT" ~ log_488_nuclear/background_WT,
        Condition == "KO" ~ log_488_nuclear/background_KO
    ))

# select WT and KO
pAur_alldata_cond_singl_log_norm_WT <- pAur_alldata_cond_singl_log_norm %>%
    filter(Condition == "WT")
pAur_alldata_cond_singl_log_norm_KO <- pAur_alldata_cond_singl_log_norm %>%
    filter(Condition == "KO")

# get values different top %
quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.95) # 13.76394
quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.96) # 13.77129
quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.97) # 13.78116
quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.98) # 13.79742
quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.99) # 13.83711

quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.95) # 13.75665
quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.96) # 13.76265
quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.97) # 13.77039
quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.98) # 13.78306
quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.99) # 13.81261

# top 5% - WT and KO
pAur_alldata_cond_WT_top5 <- pAur_alldata_cond_singl_log_norm_WT %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.95))
pAur_alldata_cond_KO_top5 <- pAur_alldata_cond_singl_log_norm_KO %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.95))

# check n
nrow(pAur_alldata_cond_WT_top5) # 3088
nrow(pAur_alldata_cond_KO_top5) # 3608

# put together
pAur_alldata_cond_top5 <- rbind(pAur_alldata_cond_WT_top5, pAur_alldata_cond_KO_top5)
nrow(pAur_alldata_cond_top5) # 6696

# factor
pAur_alldata_cond_top5$Condition <- factor(pAur_alldata_cond_top5$Condition, levels = c("WT", "KO"))

# comparisons
my_comp = list(c("WT", "KO"))

# plot
ggplot(pAur_alldata_cond_top5, aes(x=Condition, y = nuclear_488_norm, fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO" = "#e56e58")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(pAurA/B/C mean intensity) \n normalized to background", fill = "") +
    theme_classic() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          legend.position = "none"
    )
ggsave("plots/log2_aur_nuclear_normalised_background_top5.pdf", device = "pdf", width = 4, height = 4)

# top 4% - WT and KO
pAur_alldata_cond_WT_top4 <- pAur_alldata_cond_singl_log_norm_WT %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.96))
pAur_alldata_cond_KO_top4 <- pAur_alldata_cond_singl_log_norm_KO %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.96))

# check n
nrow(pAur_alldata_cond_WT_top4) # 2471
nrow(pAur_alldata_cond_KO_top4) # 2885

# put together
pAur_alldata_cond_top4 <- rbind(pAur_alldata_cond_WT_top4, pAur_alldata_cond_KO_top4)
nrow(pAur_alldata_cond_top4) # 5356

# factor
pAur_alldata_cond_top4$Condition <- factor(pAur_alldata_cond_top4$Condition, levels = c("WT", "KO"))

# comparisons
my_comp = list(c("WT", "KO"))

# plot
ggplot(pAur_alldata_cond_top4, aes(x=Condition, y = nuclear_488_norm, fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO" = "#e56e58")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(pAurA/B/C mean intensity) \n normalized to background", fill = "") +
    theme_classic() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          legend.position = "none"
    )
ggsave("plots/log2_aur_nuclear_normalised_background_top4.pdf", device = "pdf", width = 4, height = 4)

# top 3% - WT and KO
pAur_alldata_cond_WT_top3 <- pAur_alldata_cond_singl_log_norm_WT %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.97))
pAur_alldata_cond_KO_top3 <- pAur_alldata_cond_singl_log_norm_KO %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.97))

# check n
nrow(pAur_alldata_cond_WT_top3) # 1853
nrow(pAur_alldata_cond_KO_top3) # 2164

# put together
pAur_alldata_cond_top3 <- rbind(pAur_alldata_cond_WT_top3, pAur_alldata_cond_KO_top3)
nrow(pAur_alldata_cond_top3) # 4017

# factor
pAur_alldata_cond_top3$Condition <- factor(pAur_alldata_cond_top3$Condition, levels = c("WT", "KO"))

# comparisons
my_comp = list(c("WT", "KO"))

# plot
ggplot(pAur_alldata_cond_top3, aes(x=Condition, y = nuclear_488_norm, fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO" = "#e56e58")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(pAurA/B/C mean intensity) \n normalized to background", fill = "") +
    theme_classic() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          legend.position = "none"
    )
ggsave("plots/log2_aur_nuclear_normalised_background_top3.pdf", device = "pdf", width = 4, height = 4)

# top 2% - WT and KO
pAur_alldata_cond_WT_top2 <- pAur_alldata_cond_singl_log_norm_WT %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.98))
pAur_alldata_cond_KO_top2 <- pAur_alldata_cond_singl_log_norm_KO %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.98))

# check n
nrow(pAur_alldata_cond_WT_top2) # 1236
nrow(pAur_alldata_cond_KO_top2) # 1443

# put together
pAur_alldata_cond_top2 <- rbind(pAur_alldata_cond_WT_top2, pAur_alldata_cond_KO_top2)
nrow(pAur_alldata_cond_top2) # 2679

# factor
pAur_alldata_cond_top2$Condition <- factor(pAur_alldata_cond_top2$Condition, levels = c("WT", "KO"))

# comparisons
my_comp = list(c("WT", "KO"))

# plot
ggplot(pAur_alldata_cond_top2, aes(x=Condition, y = nuclear_488_norm, fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO" = "#e56e58")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(pAurA/B/C mean intensity) \n normalized to background", fill = "") +
    theme_classic() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          legend.position = "none"
    )
ggsave("plots/log2_aur_nuclear_normalised_background_top2.pdf", device = "pdf", width = 4, height = 4)

# top 1% - WT and KO
pAur_alldata_cond_WT_top1 <- pAur_alldata_cond_singl_log_norm_WT %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_WT$log_488_nuclear, 0.99))
pAur_alldata_cond_KO_top1 <- pAur_alldata_cond_singl_log_norm_KO %>%
    filter(log_488_nuclear >= quantile(pAur_alldata_cond_singl_log_norm_KO$log_488_nuclear, 0.99))

# check n
nrow(pAur_alldata_cond_WT_top1) # 618
nrow(pAur_alldata_cond_KO_top1) # 722

# put together
pAur_alldata_cond_top1 <- rbind(pAur_alldata_cond_WT_top1, pAur_alldata_cond_KO_top1)
nrow(pAur_alldata_cond_top1) # 1340

# factor
pAur_alldata_cond_top1$Condition <- factor(pAur_alldata_cond_top1$Condition, levels = c("WT", "KO"))

# comparisons
my_comp = list(c("WT", "KO"))

# plot
ggplot(pAur_alldata_cond_top1, aes(x=Condition, y = nuclear_488_norm, fill = Condition)) +
    geom_boxplot(width = 0.75, fatten = 2) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO" = "#e56e58")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(pAurA/B/C mean intensity) \n normalized to background", fill = "") +
    theme_classic() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          legend.position = "none"
    )
ggsave("plots/log2_aur_nuclear_normalised_background_top1.pdf", device = "pdf", width = 4, height = 4)