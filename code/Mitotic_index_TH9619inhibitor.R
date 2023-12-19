# libraries
library(tidyverse)
library(car)
library(ggpubr)
library(rstatix)

# Objective: Analyse mitotic index of HCT116 treated with TH9619 inhibitor

# 1. Get input data ----------------------------------------------------------

# get data
h3pser_rep1a <- read.delim("rep2 96h/h3Pser atub HAP1 MTHinh 230223__2023-02-23T13_09_29-Measurement 1/Evaluation3/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep1a)

h3pser_rep1b <- read.delim("rep2 96h/h3Pser crest HAP1 MTHinh__2023-02-23T13_37_56-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep1b)

h3pser_rep2a <- read.delim("rep3 96h/h3Pser atub HCT116 MTHinh 010323__2023-03-01T16_58_15-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep2c)

h3pser_rep2b <- read.delim("rep3 96h/h3Pser crest HCT116 MTHinh 010323__2023-03-01T17_14_36-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep2d)

h3pser_rep3a <- read.delim("rep4 96h/h3Pser atub HCT116 MTHinh 070323__2023-03-07T19_18_24-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep3c)

h3pser_rep3b <- read.delim("rep4 96h/h3Pser crest HCT116 MTHinh 070323__2023-03-07T19_29_25-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(h3pser_rep3d)

# check colnames are the same
table(colnames(h3pser_rep1a) == colnames(h3pser_rep1b))
table(colnames(h3pser_rep1a) == colnames(h3pser_rep2a))
table(colnames(h3pser_rep1a) == colnames(h3pser_rep2b))
table(colnames(h3pser_rep1a) == colnames(h3pser_rep3a))
table(colnames(h3pser_rep1a) == colnames(h3pser_rep3b))
# add replicate
h3pser_rep1a <- h3pser_rep1a %>%
    mutate(Replicate = rep(1, nrow(h3pser_rep1a)),
           Plate = rep("HAP1_HCT116", nrow(h3pser_rep1a)))
h3pser_rep1b <- h3pser_rep1b %>%
    mutate(Replicate = rep(1, nrow(h3pser_rep1b)),
           Plate = rep("HAP1_HCT116", nrow(h3pser_rep1b)))

h3pser_rep2a <- h3pser_rep2a %>%
    mutate(Replicate = rep(2, nrow(h3pser_rep2a)),
           Plate = rep("HCT116", nrow(h3pser_rep2a)))
h3pser_rep2b <- h3pser_rep2d %>%
    mutate(Replicate = rep(2, nrow(h3pser_rep2b)),
           Plate = rep("HCT116", nrow(h3pser_rep2b)))

h3pser_rep3a <- h3pser_rep3a %>%
    mutate(Replicate = rep(3, nrow(h3pser_rep3a)),
           Plate = rep("HCT116", nrow(h3pser_rep3a)))
h3pser_rep3b <- h3pser_rep3b %>%
    mutate(Replicate = rep(3, nrow(h3pser_rep3b)),
           Plate = rep("HCT116", nrow(h3pser_rep3b)))

# put all data together -> since colnames are the same, do rbind
h3pser_alldata <- rbind(h3pser_rep1a, h3pser_rep1b,
                        h3pser_rep2a, h3pser_rep2b, 
                        h3pser_rep3a, h3pser_rep3b)

# check
table(h3pser_alldata$Replicate, h3pser_alldata$Row) 
table(h3pser_alldata$Replicate, h3pser_alldata$Column) 

# 2. Process data ----------------------------------------------------------

# add condition
mitosis_quant_cond <- h3pser_alldata %>%
    mutate(Cell_Line = case_when(
        (Row == 2 | Row == 3) & Plate == "HCT116" ~ "HCT116",
        (Row == 2 | Row == 3) & Plate == "HAP1_HCT116" ~ "HAP1-WT",
        (Row == 4 | Row == 5) & Plate == "HAP1_HCT116" ~ "HAP1-mut",
        (Row == 6 | Row == 7) & Plate == "HAP1_HCT116" ~ "HCT116"),
        Drug = case_when(
            (Column == 2 | Column == 3) ~ 0,
            ((Column == 4 | Column == 5) & Replicate != 1) ~ 10,
            ((Column == 4 | Column == 5) & Replicate == 1) ~ 25,
            ((Column == 6 | Column == 7) & Replicate != 1) ~ 25,
            ((Column == 6  | Column == 7) & Replicate == 1) ~ 63,
            ((Column == 8  | Column == 9) & Replicate != 1) ~ 63,
            ((Column == 8 | Column == 9) & Replicate == 1) ~ 156,
            ((Column == 10 | Column == 11) & Replicate != 1) ~ 156,
            ((Column == 10 | Column == 11) & Replicate == 1) ~ 391))

# check
table(mitosis_quant_cond$Replicate, mitosis_quant_cond$Cell_Line) 
table(mitosis_quant_cond$Row, mitosis_quant_cond$Cell_Line) 
table(mitosis_quant_cond$Replicate, mitosis_quant_cond$Drug)  
table(mitosis_quant_cond$Column, mitosis_quant_cond$Drug)  

# concentrations 25 and 63
mitosis_quant_cond2 <- mitosis_quant_cond %>%
    filter(Drug %in% c(0,25,63))
table(mitosis_quant_cond2$Column, mitosis_quant_cond2$Drug)  

# select HCT116 data
mitosis_quant_cond_HCT <- mitosis_quant_cond2 %>%
    filter(Cell_Line == "HCT116")

# filter single nuclei
h3pser_alldata_singl_HCT <- mitosis_quant_cond_HCT %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200)

nrow(mitosis_quant_cond_HCT) # 177671
nrow(h3pser_alldata_singl_HCT) # 169923
table(h3pser_alldata_singl_HCT$Cell_Line, h3pser_alldata_singl_HCT$Drug)
table(h3pser_alldata_singl_HCT$Cell_Line, h3pser_alldata_singl_HCT$Replicate)

# get interesting data and log2 intensity
logh3pser_alldata_HCT <- h3pser_alldata_singl_HCT %>%
    select(Cell_Line, Drug, Replicate, All.Nuclei.Selected...H3phosphoSer10.Mean) %>%
    mutate(log2H3pSer = log2(All.Nuclei.Selected...H3phosphoSer10.Mean),
           log10H3pSer = log10(All.Nuclei.Selected...H3phosphoSer10.Mean))

# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_H3pSer_data <- function(n) {
    repdata <- logh3pser_alldata_HCT %>%
        filter(Replicate == n) 
    repdata_scale <- repdata %>%
        mutate(H3pSer_sc = range01(All.Nuclei.Selected...H3phosphoSer10.Mean),
               log2H3pSer_sc = range01(log2H3pSer))
}

reps <- 1:3
h3pser_alldata_scaled_reps_HCT <- map(reps, scale_H3pSer_data)

# build data frame
h3pser_alldata_scaled_HCT <- bind_rows(h3pser_alldata_scaled_reps_HCT)
table(h3pser_alldata_scaled_HCT$Cell_Line)
table(h3pser_alldata_scaled_HCT$Drug)

# plot histograms
ggplot(h3pser_alldata_scaled_HCT, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("global_plots_96h/histogram_all_values_scaled_HCT.pdf")

mitosis_quant_low <- h3pser_alldata_scaled_HCT %>%
    filter(H3pSer_sc <= 0.25)

ggplot(mitosis_quant_low, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("global_plots_96h/histogram_low_values_scaled_HCT.pdf")

mitosis_quant_med <- h3pser_alldata_scaled_HCT %>%
    filter(H3pSer_sc > 0.25 & H3pSer_sc < 0.75)

ggplot(mitosis_quant_med, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("global_plots_96h/histogram_med_values_scaled_HCT.pdf")

mitosis_quant_high <- h3pser_alldata_scaled_HCT %>%
    filter(H3pSer_sc >= 0.75)

ggplot(mitosis_quant_high, aes(x=H3pSer_sc)) +
    geom_histogram(bins=100) +
    xlab("scaled H3pSer10 mean intensity scaled") +
    theme_classic()
ggsave("global_plots_96h/histogram_high_values_scaled_HCT.pdf")

# 3. Get % of mitotic cells ---------------------------------------------------------------------------------------

# get mitotic cells with threshold
mitosis_quant_thresmit_HCT <- h3pser_alldata_scaled_HCT %>%
    mutate(Mitotic = ifelse(H3pSer_sc > 0.17, T, F))
table(mitosis_quant_thresmit_HCT$Cell_Line, mitosis_quant_thresmit_HCT$Mitotic)

# convert condition to factor
mitosis_quant_thresmit_HCT$Drug <- factor(mitosis_quant_thresmit_HCT$Drug)

# average and sd for each treatment
mitosis_quant_avg_sc_HCT <- mitosis_quant_thresmit_HCT %>%
    group_by(Cell_Line, Drug, Replicate) %>%
    summarise(
        count = n(),
        mean = mean(H3pSer_sc),
        median = median(H3pSer_sc),
        sd = sd(H3pSer_sc),
        se = sd/sqrt(count),
        perc_mitotic = sum(Mitotic)/count*100,
        perc_mitotic2 = sum(Mitotic2)/count*100,
        perc_mitotic3 = sum(Mitotic3)/count*100,
        perc_mitotic4 = sum(Mitotic4)/count*100)

mitosis_quant_avg_sc_HCT$Drug <- factor(mitosis_quant_avg_sc_HCT$Drug)

# decide threshold
mean(unlist(mitosis_quant_avg_sc_HCT[mitosis_quant_avg_sc_HCT$Cell_Line == "HCT116" & 
                                         mitosis_quant_avg_sc_HCT$Drug == 0 , "median"])) # 0.003723032

mean(unlist(mitosis_quant_avg_sc_HCT[mitosis_quant_avg_sc_HCT$Cell_Line == "HCT116" & 
                                         mitosis_quant_avg_sc_HCT$Drug == 0 , "sd"])) # 0.05534305

# plot HCT116 threshold
ggplot(mitosis_quant_thresmit_HCT, aes(x = Drug, y = H3pSer_sc, col = Mitotic)) +
    geom_jitter(alpha = 0.4) +
    geom_hline(yintercept=0.17) +
    facet_wrap(~Replicate) +
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
ggsave("global_plots_96h/H3pS3_quantif_scaled_HCT.pdf", width = 7, height = 3)

# 4. Statistics and plot ---------------------------------------------------------------------------------------

# data of perc_mitotic
mitotic_perc_info_sc_HCT <- as.data.frame(mitosis_quant_avg_sc_HCT[, c("Cell_Line", "Drug",
                                                                       "Replicate", "perc_mitotic")])

# normalise with DMSO condition of same replicate
rep1_DMSO <- mitotic_perc_info_sc_HCT$perc_mitotic[mitotic_perc_info_sc_HCT$Drug == 0 &
                                                       mitotic_perc_info_sc_HCT$Replicate == 1]
rep2_DMSO <- mitotic_perc_info_sc_HCT$perc_mitotic[mitotic_perc_info_sc_HCT$Drug == 0 &
                                                       mitotic_perc_info_sc_HCT$Replicate == 2]
rep3_DMSO <- mitotic_perc_info_sc_HCT$perc_mitotic[mitotic_perc_info_sc_HCT$Drug == 0 &
                                                       mitotic_perc_info_sc_HCT$Replicate == 3]
mitotic_perc_info_sc_HCT_norm <- mitotic_perc_info_sc_HCT %>%
    mutate(perc_mitosis_norm = case_when(
        Replicate == 1 ~ perc_mitotic/rep1_DMSO,
        Replicate == 2 ~ perc_mitotic/rep2_DMSO,
        Replicate == 3 ~ perc_mitotic/rep3_DMSO))

# statistic summary
mitotic_perc_stats_sc_HCT <- mitotic_perc_info_sc_HCT %>%
    group_by(Drug) %>%
    summarise(
        count = n(),
        mean_mit = mean(perc_mitotic),
        median_mit= median(perc_mitotic),
        sd_mit = sd(perc_mitotic),
        se_mit = sd_mit/sqrt(count))

mitotic_perc_stats_sc_HCT_norm <- mitotic_perc_info_sc_HCT_norm %>%
    group_by(Drug) %>%
    summarise(
        count = n(),
        mean_mit = mean(perc_mitosis_norm),
        median_mit= median(perc_mitosis_norm),
        sd_mit = sd(perc_mitosis_norm),
        se_mit = sd_mit/sqrt(count))

# check normality of groups 
group_25 <- mitotic_perc_info_sc_HCT_norm %>%
    filter(Drug == 25)
group_63 <- mitotic_perc_info_sc_HCT_norm %>%
    filter(Drug == 63)
shapiro.test(group_25$perc_mitosis_norm)
shapiro.test(group_63$perc_mitosis_norm)

# check homogeneity of variance -> yes
leveneTest(perc_mitosis_norm ~ Drug, data=mitotic_perc_info_sc_HCT_norm)

# one sample t.test
t.test(mitotic_perc_info_sc_HCT_norm$perc_mitosis_norm[mitotic_perc_info_sc_HCT_norm$Drug == 25], 
       mu = 1, alternative = "two.sided") # pval 0.05397
t.test(mitotic_perc_info_sc_HCT_norm$perc_mitosis_norm[mitotic_perc_info_sc_HCT_norm$Drug == 63], 
       mu = 1, alternative = "two.sided") # pval 0.006727

# ggplot data
ggplot() +
    geom_errorbar(data = mitotic_perc_stats_sc_HCT_norm, aes(x=Drug,
                                                             ymin=mean_mit-sd_mit,
                                                             ymax=mean_mit+sd_mit), 
                  width=.2,
                  position=position_dodge(0.05)) +
    geom_col(data = mitotic_perc_stats_sc_HCT_norm, aes(x=Drug, y=mean_mit), 
             fill = "#aaafb0", width = 0.5)+
    geom_jitter(data = mitotic_perc_info_sc_HCT_norm, aes(x=Drug, y=perc_mitosis_norm),
                height = 0, width = 0.1, alpha = 1, size = 1.5) +
    xlab("TH9619 Inhibitor (nM)") +
    ylab("Proportion of mitotic cells \n normalized to DMSO condition") +
    scale_y_continuous(expand=c(0,0), limits = c(0,1.5)) +
    theme_classic() +
    theme(axis.text = element_text(size = 11), 
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 12))
ggsave("global_plots_96h/perc_mitotic_cells_HCT_normDMSO.pdf", width = 4, height = 3)

