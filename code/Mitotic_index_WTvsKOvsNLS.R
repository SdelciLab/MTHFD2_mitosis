# libraries
library(tidyverse)
library(car)
library(ggpubr)
library(rstatix)

# Objective: Analyse mitotic index between MTHFD2-WT and MTHFD2-KO and MTHFD2-NLS

# 1. Get input data -------------------------------------------------------------

# get data
h3pser_rep1a <- read.delim("rep3/h3Pser atub WTvsNLS40X 010323__2023-03-01T17_45_22-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep1a)

h3pser_rep1b <- read.delim("rep3/h3Pser crest WTvsNLS 40X 010323__2023-03-01T18_03_50-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep1b)

h3pser_rep2a <- read.delim("rep4 old/h3Pser atub WTvsNLS40X old 010323__2023-03-01T18_41_44-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep2a)

h3pser_rep2b <- read.delim("rep4 old/h3Pser crest WTvsNLS 40X old 010323__2023-03-01T18_33_16-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep2b)

h3pser_rep3a <- read.delim("rep5/h3Pser atub WTvsNLS40X 080323__2023-03-08T15_30_07-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep3a)

h3pser_rep3b <- read.delim("rep5/h3Pser crest WTvsNLS 40X 080321__2023-03-08T15_43_00-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep3b)

# check colnames are the same
table(colnames(h3pser_rep1a) == colnames(h3pser_rep1b))
table(colnames(h3pser_rep1a) == colnames(h3pser_rep2a))
table(colnames(h3pser_rep1a) == colnames(h3pser_rep2b))
table(colnames(h3pser_rep1a) == colnames(h3pser_rep3a))
table(colnames(h3pser_rep1a) == colnames(h3pser_rep3b))

# add replicate
h3pser_rep1a <- h3pser_rep1a %>%
    mutate(Replicate = rep(1, nrow(h3pser_rep1a)))
h3pser_rep1b <- h3pser_rep1b %>%
    mutate(Replicate = rep(1, nrow(h3pser_rep1b)))
h3pser_rep2a <- h3pser_rep2a %>%
    mutate(Replicate = rep(2, nrow(h3pser_rep2a)))
h3pser_rep2b <- h3pser_rep2b %>%
    mutate(Replicate = rep(2, nrow(h3pser_rep2b)))
h3pser_rep3a <- h3pser_rep3a %>%
    mutate(Replicate = rep(3, nrow(h3pser_rep3a)))
h3pser_rep3b <- h3pser_rep3b %>%
    mutate(Replicate = rep(3, nrow(h3pser_rep3b)))

# put all data together -> since colnames are the same, do rbind
h3pser_alldata <- rbind(h3pser_rep1a, h3pser_rep1b, h3pser_rep2a, h3pser_rep2b, 
                        h3pser_rep3a, h3pser_rep3b)

# check
table(h3pser_alldata$Replicate, h3pser_alldata$Row) 
table(h3pser_alldata$Replicate, h3pser_alldata$Column) 

# 2. Process data --------------------------------------------------------------

# add condition
h3pser_alldata_cond <- h3pser_alldata %>%
    mutate(Condition = case_when(
        Row == 2 ~ "WT",
        Row == 3 ~ "KO1",
        Row == 4 ~ "KO2",
        Row == 5 ~ "NLS",
    ))

# filter single nuclei
h3pser_alldata_singl <- h3pser_alldata_cond %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200)

nrow(h3pser_alldata_cond) 
nrow(h3pser_alldata_singl)
table(h3pser_alldata_singl$Condition)

# get interesting data and log2 intensity
logh3pser_alldata <- h3pser_alldata_singl %>%
    select(Condition, Replicate, All.Nuclei.Selected...H3phosphoSer10.Mean) %>%
    mutate(log2H3pSer = log2(All.Nuclei.Selected...H3phosphoSer10.Mean),
           log10H3pSer = log10(All.Nuclei.Selected...H3phosphoSer10.Mean))

# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_H3pSer_data <- function(n) {
    repdata <- logh3pser_alldata %>%
        filter(Replicate == n) 
    repdata_scale <- repdata %>%
        mutate(H3pSer_sc = range01(All.Nuclei.Selected...H3phosphoSer10.Mean),
               log2H3pSer_sc = range01(log2H3pSer))
}

reps <- 1:3
h3pser_alldata_scaled_reps <- map(reps, scale_H3pSer_data)

# build data frame
h3pser_alldata_scaled <- bind_rows(h3pser_alldata_scaled_reps)
table(h3pser_alldata_scaled$Condition)

# plot again histograms
ggplot(h3pser_alldata_scaled, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("reps_plots/histogram_all_values_scaled.pdf")

mitosis_quant_low <- h3pser_alldata_scaled %>%
    filter(H3pSer_sc <= 0.25)

ggplot(mitosis_quant_low, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("reps_plots/histogram_low_values_scaled.pdf")

mitosis_quant_med <- h3pser_alldata_scaled %>%
    filter(H3pSer_sc > 0.25 & H3pSer_sc < 0.75)

ggplot(mitosis_quant_med, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("reps_plots/histogram_med_values_scaled.pdf")

mitosis_quant_high <- h3pser_alldata_scaled %>%
    filter(H3pSer_sc >= 0.75)

ggplot(mitosis_quant_high, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity scaled") +
    theme_classic()
ggsave("reps_plots/histogram_high_values_scaled.pdf")

# 3. Get % of mitotic cells ----------------------------------------------

# get mitotic cells with threshold
mitosis_quant_thresmit <- h3pser_alldata_scaled %>%
    mutate(Mitotic = ifelse(H3pSer_sc > 0.22, T, F))
table(mitosis_quant_thresmit$Condition, mitosis_quant_thresmit$Mitotic)

# convert condition to factor
mitosis_quant_thresmit$Condition <- factor(mitosis_quant_thresmit$Condition, levels = c("WT", "KO1", "KO2", "NLS"))

# average and sd for each treatment
mitosis_quant_avg_sc <- mitosis_quant_thresmit %>%
    group_by(Condition, Replicate) %>%
    summarise(
        count = n(),
        mean = mean(H3pSer_sc),
        median = median(H3pSer_sc),
        sd = sd(H3pSer_sc),
        se = sd/sqrt(count),
        perc_mitotic = sum(Mitotic)/count*100)

mitosis_quant_avg_sc$Condition <- factor(mitosis_quant_avg_sc$Condition, levels = c("WT", "KO1", "KO2", "NLS"))

# decide threshold
mean(unlist(mitosis_quant_avg_sc[mitosis_quant_avg_sc$Condition == "WT", "median"])) # 0.002947436
mean(unlist(mitosis_quant_avg_sc[mitosis_quant_avg_sc$Condition == "WT", "sd"])) # 0.07246142

# plot
ggplot(mitosis_quant_thresmit, aes(x = Condition, y = H3pSer_sc, col = Mitotic)) +
    geom_jitter(alpha = 0.4, stroke = NA) +
    geom_hline(yintercept=0.22) +
    facet_wrap(~Replicate, nrow = 1) +
    scale_color_manual(values = c("TRUE" = "#559A73", "FALSE" = "lightgrey")) +
    ylab("Scaled H3pSer10 mean intensity") +
    theme_bw() +
    theme(axis.text = element_text(size = 10), 
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    )
ggsave("reps_plots/H3pS3_quantif_scaled.pdf", width = 7, height = 3)

# 4. Statistical analysis and plot ------------------------------------------------

# data of perc_mitotic
mitotic_perc_info_sc <- as.data.frame(mitosis_quant_avg_sc[, c("Condition", "Replicate", "perc_mitotic")])

# statistic summary
mitotic_perc_stats_sc <- mitotic_perc_info_sc %>%
    group_by(Condition) %>%
    summarise(
        count = n(),
        mean_mit = mean(perc_mitotic),
        median_mit= median(perc_mitotic),
        sd_mit = sd(perc_mitotic),
        se_mit = sd_mit/sqrt(count))

# check normality of groups -> all normal 
group_WT <- mitotic_perc_info_sc %>%
    filter(Condition == "WT")
group_KO1 <- mitotic_perc_info_sc %>%
    filter(Condition == "KO1")
group_KO2 <- mitotic_perc_info_sc %>%
    filter(Condition == "KO2")
group_NLS <- mitotic_perc_info_sc %>%
    filter(Condition == "NLS")
shapiro.test(group_WT$perc_mitotic)
shapiro.test(group_KO1$perc_mitotic)
shapiro.test(group_KO2$perc_mitotic)
shapiro.test(group_NLS$perc_mitotic)

# check homogeneity of variance -> yes
leveneTest(perc_mitotic ~ Condition, data=mitotic_perc_info_sc)

# t-test between WT, KO1 and WT, KO2
group_WT_KO1 <- mitotic_perc_info_sc %>%
    filter(Condition %in% c("WT", "KO1"))
group_WT_KO2 <- mitotic_perc_info_sc %>%
    filter(Condition %in% c("WT", "KO2"))
group_WT_NLS <- mitotic_perc_info_sc %>%
    filter(Condition %in% c("WT", "NLS"))

t.test(perc_mitotic ~ Condition, data=group_WT_KO1, var.equal = TRUE) 
t.test(perc_mitotic ~ Condition, data=group_WT_KO2, var.equal = TRUE) 
t.test(perc_mitotic ~ Condition, data=group_WT_NLS, var.equal = TRUE) 

# test
ttest_pvals_all <- compare_means(perc_mitotic ~ Condition, data = mitotic_perc_info_sc, ref.group = "WT", 
                                 method = "t.test",  var.equal = TRUE) %>%
    mutate(y.position = c(3.2, 3.4, 3.6))

# ggplot data
mitotic_perc_info_sc$Condition <- factor(mitotic_perc_info_sc$Condition, levels = c("WT", "NLS", "KO1", "KO2"))
mitotic_perc_stats_sc$Condition <- factor(mitotic_perc_stats_sc$Condition, levels = c("WT", "NLS", "KO1", "KO2"))

ggplot() +
    geom_errorbar(data = mitotic_perc_stats_sc, aes(x=Condition, ymin=mean_mit-sd_mit, ymax=mean_mit+sd_mit), 
                  width=.2,
                  position=position_dodge(0.05)) +
    geom_col(data = mitotic_perc_stats_sc, aes(x=Condition, y=mean_mit, fill = Condition), width = 0.5)+
    geom_jitter(data = mitotic_perc_info_sc, aes(x=Condition, y=perc_mitotic), height = 0, width = 0.1, alpha = 1, size = 1.5)+
    stat_pvalue_manual(ttest_pvals_all, label = "p.format", size = 3) +
    xlab("") +
    ylab("Percentage of mitotic cells (%)") +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58", "KO2" = "#f9beb3", "NLS" = "#890808")) +
    theme_classic() +
    theme(axis.text = element_text(size = 11), 
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 12))
ggsave("reps_plots/perc_mitotic_cells_bar_pval_scaled.pdf", width = 5, height = 3.5)
