# libraries
library(tidyverse)
library(car)
library(ggpubr)
library(rstatix)

# Objective: analyse mitotic index MTHFD2-WT vs MTHFD2-KO cells

# 1. Get input data ---------------------------------------------------------------

# get data
h3pser_rep1 <- read.delim("H3pS3 tub MTHFD2 WT KO KD 40X__2022-02-10T16_01_26-Measurement 1/Evaluation9/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep1)

h3pser_rep2 <- read.delim("H3pS10 atub CREST WT KO KD 40X 030322__2022-03-03T10_30_34-Measurement 1/Evaluation4/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep2)

h3pser_rep3a <- read.delim("h3pser3 atub WT KO KD 070422__2022-04-07T14_40_23-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep3a)

h3pser_rep3b <- read.delim("h3pser3 crest WT KO KD 070422__2022-04-07T15_00_20-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep3b)

h3pser_rep4a <- read.delim("h3pser3 atub WT KO KD 210422__2022-04-21T11_21_08-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep4a)

h3pser_rep4b <- read.delim("h3pser3 crest WT KO KD 210422__2022-04-21T11_35_02-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep4b)

h3pser_rep5a <- read.delim("H3pS3 tub MTHFD2 WT KO KD 40X 110522__2022-05-11T14_41_28-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep5a)

h3pser_rep5b <- read.delim("h3pser3 crest WT KO KD 110522__2022-05-11T14_34_12-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep5b)

# check colnames are the same
table(colnames(h3pser_rep1) == colnames(h3pser_rep2))
table(colnames(h3pser_rep2) == colnames(h3pser_rep3a))
table(colnames(h3pser_rep2) == colnames(h3pser_rep3b))
table(colnames(h3pser_rep2) == colnames(h3pser_rep4a))
table(colnames(h3pser_rep2) == colnames(h3pser_rep4b))
table(colnames(h3pser_rep2) == colnames(h3pser_rep5a))
table(colnames(h3pser_rep2) == colnames(h3pser_rep5b))

# change colnames of replicate 1 - H3phosphoSer3 -> H3phosphoSer10
colnames(h3pser_rep1) <- colnames(h3pser_rep2)
colnames(h3pser_rep1) 

# add replicate
h3pser_rep1 <- h3pser_rep1 %>%
    mutate(Replicate = rep("rep1", nrow(h3pser_rep1)))
h3pser_rep2 <- h3pser_rep2 %>%
    mutate(Replicate = rep("rep2", nrow(h3pser_rep2)))
h3pser_rep3a <- h3pser_rep3a %>%
    mutate(Replicate = rep("rep3", nrow(h3pser_rep3a)))
h3pser_rep3b <- h3pser_rep3b %>%
    mutate(Replicate = rep("rep3", nrow(h3pser_rep3b)))
h3pser_rep4a <- h3pser_rep4a %>%
    mutate(Replicate = rep("rep4", nrow(h3pser_rep4a)))
h3pser_rep4b <- h3pser_rep4b %>%
    mutate(Replicate = rep("rep4", nrow(h3pser_rep4b)))
h3pser_rep5a <- h3pser_rep5a %>%
    mutate(Replicate = rep("rep5", nrow(h3pser_rep5a)))
h3pser_rep5b <- h3pser_rep5b %>%
    mutate(Replicate = rep("rep5", nrow(h3pser_rep5b)))

# put all data together -> since colnames are the same, do rbind
h3pser_alldata <- rbind(h3pser_rep1, h3pser_rep2, h3pser_rep3a, h3pser_rep3b, h3pser_rep4a, h3pser_rep4b, h3pser_rep5a, h3pser_rep5b)

# check
table(h3pser_alldata$Replicate, h3pser_alldata$Row) 
table(h3pser_alldata$Replicate, h3pser_alldata$Column) 

# 2. Process data -------------------------------------------------------------------------------------------

# add condition
h3pser_alldata_cond <- h3pser_alldata %>%
    mutate(Condition = case_when(
        Row == 2 ~ "WT",
        Row == 3 ~ "KO1",
        Row == 4 ~ "KO2"
    ))

# select WT, KO1 and KO2
h3pser_alldata_int <- h3pser_alldata_cond %>%
    filter(Condition %in% c("WT", "KO1", "KO2"))

# check
table(h3pser_alldata_int$Replicate, h3pser_alldata_int$Row) # rows 1-WT 2-KO1 3-KO2 
table(h3pser_alldata_int$Replicate, h3pser_alldata_int$Column)

# filter single nuclei
h3pser_alldata_singl <- h3pser_alldata_int %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200)

nrow(h3pser_alldata_int) # 97576
nrow(h3pser_alldata_singl) # 88223
table(h3pser_alldata_singl$Condition)

# get interesting data and log2 intensity
logh3pser_alldata <- h3pser_alldata_singl %>%
    select(Condition, Replicate, All.Nuclei.Selected...H3phosphoSer10.Mean) %>%
    mutate(log2H3pSer = log2(All.Nuclei.Selected...H3phosphoSer10.Mean),
           log10H3pSer = log10(All.Nuclei.Selected...H3phosphoSer10.Mean))

# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

rep1data <- logh3pser_alldata %>%
    filter(Replicate == "rep1") 
rep1data_scale <- rep1data %>%
    mutate(H3pSer_sc = range01(All.Nuclei.Selected...H3phosphoSer10.Mean),
           log2H3pSer_sc = range01(log2H3pSer))

rep2data <- logh3pser_alldata %>%
    filter(Replicate == "rep2") 
rep2data_scale <- rep2data %>%
    mutate(H3pSer_sc = range01(All.Nuclei.Selected...H3phosphoSer10.Mean),
           log2H3pSer_sc = range01(log2H3pSer))

rep3data <- logh3pser_alldata %>%
    filter(Replicate == "rep3") 
rep3data_scale <- rep3data %>%
    mutate(H3pSer_sc = range01(All.Nuclei.Selected...H3phosphoSer10.Mean),
           log2H3pSer_sc = range01(log2H3pSer))

rep4data <- logh3pser_alldata %>%
    filter(Replicate == "rep4") 
rep4data_scale <- rep4data %>%
    mutate(H3pSer_sc = range01(All.Nuclei.Selected...H3phosphoSer10.Mean),
           log2H3pSer_sc = range01(log2H3pSer))

rep5data <- logh3pser_alldata %>%
    filter(Replicate == "rep5") 
rep5data_scale <- rep5data %>%
    mutate(H3pSer_sc = range01(All.Nuclei.Selected...H3phosphoSer10.Mean),
           log2H3pSer_sc = range01(log2H3pSer))

# build data frame
h3pser_alldata_scaled <- rbind(rep1data_scale, rep2data_scale, rep3data_scale, rep4data_scale, rep5data_scale)
nrow(h3pser_alldata_scaled)
nrow(rep1data_scale)
nrow(rep2data_scale)
nrow(rep3data_scale)
nrow(rep4data_scale)
nrow(rep5data_scale)
table(h3pser_alldata_scaled$Condition)

# plot histograms
ggplot(h3pser_alldata_scaled, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("plots/histogram_all_values_scaled.pdf")

mitosis_quant_low <- h3pser_alldata_scaled %>%
    filter(H3pSer_sc <= 0.25)

ggplot(mitosis_quant_low, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("plots/histogram_low_values_scaled.pdf")

mitosis_quant_med <- h3pser_alldata_scaled %>%
    filter(H3pSer_sc > 0.25 & H3pSer_sc < 0.75)

ggplot(mitosis_quant_med, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("plots/histogram_med_values_scaled.pdf")

mitosis_quant_high <- h3pser_alldata_scaled %>%
    filter(H3pSer_sc >= 0.75)

ggplot(mitosis_quant_high, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity scaled") +
    theme_classic()
ggsave("plots/histogram_high_values_scaled.pdf")

# 3. Get % of mitotic cells -----------------------------------------------------------------

# get mitotic cells with threshold
mitosis_quant_thresmit <- h3pser_alldata_scaled %>%
    mutate(Mitotic = ifelse(H3pSer_sc > 0.23, T, F))
table(mitosis_quant_thresmit$Condition, mitosis_quant_thresmit$Mitotic)

# convert condition to factor
mitosis_quant_thresmit$Condition <- factor(mitosis_quant_thresmit$Condition, levels = c("WT", "KO1", "KO2"))

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

mitosis_quant_avg_sc$Condition <- factor(mitosis_quant_avg_sc$Condition, levels = c("WT", "KO1", "KO2"))

# decide threshold
mean(unlist(mitosis_quant_avg_sc[mitosis_quant_avg_sc$Condition == "WT", "median"])) # 0.001935692
mean(unlist(mitosis_quant_avg_sc[mitosis_quant_avg_sc$Condition == "WT", "sd"])) # 0.07618692

# plot
ggplot(mitosis_quant_thresmit, aes(x = Condition, y = H3pSer_sc, col = Mitotic)) +
    geom_jitter(alpha = 0.4, stroke = NA) +
    geom_hline(yintercept=0.23) +
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
ggsave("plots/H3pS3_quantif_scaled.pdf", width = 7, height = 3)

# 4. Statistical analysis and plot  -----------------------------------------------------------------

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
shapiro.test(group_WT$perc_mitotic)
shapiro.test(group_KO1$perc_mitotic)
shapiro.test(group_KO2$perc_mitotic)

# check homogeneity of variance -> yes
leveneTest(perc_mitotic ~ Condition, data=mitotic_perc_info_sc)

# t-test between WT, KO1 and WT, KO2
group_WT_KO1 <- mitotic_perc_info_sc %>%
    filter(Condition %in% c("WT", "KO1"))
group_WT_KO2 <- mitotic_perc_info_sc %>%
    filter(Condition %in% c("WT", "KO2"))

t.test(perc_mitotic ~ Condition, data=group_WT_KO1, var.equal = TRUE) # yes
t.test(perc_mitotic ~ Condition, data=group_WT_KO2, var.equal = TRUE) # yes

# test
ttest_pvals_all <- compare_means(perc_mitotic ~ Condition, data = mitotic_perc_info_sc,
                                 ref.group = "WT", method = "t.test",  var.equal = TRUE) %>%
    mutate(y.position = c(3, 3.2))

# ggplot data
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
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58", "KO2" = "#f9beb3")) +
    theme_classic() +
    theme(axis.text = element_text(size = 11), 
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 12))
ggsave("plots/perc_mitotic_cells_bar_pval_scaled.pdf", width = 5, height = 3.5)

