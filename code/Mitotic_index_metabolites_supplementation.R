# libraries
library(tidyverse)
library(car)
library(ggpubr)
library(rstatix)

# Objective: Calculate mitotic index in MTHFD2 WT and KO cells after metabolites supplementation

# 1. Get input data -------------------------------------------------------------

# get data low
metab_rep1aL <- read.delim("Low_levels_181122/h3Pser atub metab 40X low levels 171122__2022-11-17T18_50_08-Measurement 1/Evaluation4/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep1aL)

metab_rep1bL <- read.delim("Low_levels_181122/h3Pser crest metab 40X low levels 171122__2022-11-17T19_30_34-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep1bL)

metab_rep2aL <- read.delim("Low_levels_2ndrep/h3Pser atub metab 40X low 2ndrep__2022-12-14T16_01_11-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep2aL)

metab_rep2bL <- read.delim("Low_levels_2ndrep/h3Pser crest metab 40X low levels 2ndrep__2022-12-14T16_38_43-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep2bL)

metab_rep3aL <- read.delim("Low_levels_3rdrep/h3Pser atub metab 40X low 3rd rep__2022-12-15T14_52_03-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep3aL)

metab_rep3bL <- read.delim("Low_levels_3rdrep/h3Pser crest metab 40X low 3rd rep__2022-12-15T17_16_38-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep3bL)

# get data high
metab_rep1aH <- read.delim("High_levels_181122/h3Pser atub metab 40X high levels 171122__2022-11-17T18_09_38-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep1aH)

metab_rep1bH <- read.delim("High_levels_181122/h3Pser crest metab 40X high levels 171122__2022-11-17T20_18_07-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep1bH)

metab_rep2aH <- read.delim("High_levels_2ndrep/h3Pser atub metab 40X interm real 2ndrep__2022-12-14T13_56_12-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep2aH)

metab_rep2bH <- read.delim("High_levels_2ndrep/h3Pser crest metab 40X interm levels 2ndrep__2022-12-14T17_11_45-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep2bH)

metab_rep3aH <- read.delim("High_levels_3rdrep/h3Pser atub metab 40X medium 3rd rep__2022-12-15T16_05_45-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep3aH)

metab_rep3bH <- read.delim("High_levels_3rdrep/h3Pser crest metab 40X medium 3rd rep__2022-12-15T16_42_49-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(metab_rep3bH)

# check colnames are the same
table(colnames(metab_rep1aL) == colnames(metab_rep1bL))
table(colnames(metab_rep1aL) == colnames(metab_rep2aL))
table(colnames(metab_rep1aL) == colnames(metab_rep2bL))
table(colnames(metab_rep1aL) == colnames(metab_rep3aL))
table(colnames(metab_rep1aL) == colnames(metab_rep3bL))
table(colnames(metab_rep1aL) == colnames(metab_rep1aH))
table(colnames(metab_rep1aL) == colnames(metab_rep1bH))
table(colnames(metab_rep1aL) == colnames(metab_rep2aH))
table(colnames(metab_rep1aL) == colnames(metab_rep2bH))
table(colnames(metab_rep1aL) == colnames(metab_rep3aH))
table(colnames(metab_rep1aL) == colnames(metab_rep3bH))

# add replicate and level
metab_rep1aL <- metab_rep1aL %>%
    mutate(Replicate = rep(1, nrow(metab_rep1aL)),
           Level = rep("Low", nrow(metab_rep1aL)))
metab_rep1bL <- metab_rep1bL %>%
    mutate(Replicate = rep(1, nrow(metab_rep1bL)),
           Level = rep("Low", nrow(metab_rep1bL)))

metab_rep2aL <- metab_rep2aL %>%
    mutate(Replicate = rep(2, nrow(metab_rep2aL)),
           Level = rep("Low", nrow(metab_rep2aL)))
metab_rep2bL <- metab_rep2bL %>%
    mutate(Replicate = rep(2, nrow(metab_rep2bL)),
           Level = rep("Low", nrow(metab_rep2bL)))

metab_rep3aL <- metab_rep3aL %>%
    mutate(Replicate = rep(3, nrow(metab_rep3aL)),
           Level = rep("Low", nrow(metab_rep3aL)))
metab_rep3bL <- metab_rep3bL %>%
    mutate(Replicate = rep(3, nrow(metab_rep3bL)),
           Level = rep("Low", nrow(metab_rep3bL)))

metab_rep1aH <- metab_rep1aH %>%
    mutate(Replicate = rep(1, nrow(metab_rep1aH)),
           Level = rep("High", nrow(metab_rep1aH)))
metab_rep1bH <- metab_rep1bH %>%
    mutate(Replicate = rep(1, nrow(metab_rep1bH)),
           Level = rep("High", nrow(metab_rep1bH)))

metab_rep2aH <- metab_rep2aH %>%
    mutate(Replicate = rep(2, nrow(metab_rep2aH)),
           Level = rep("High", nrow(metab_rep2aH)))
metab_rep2bH <- metab_rep2bH %>%
    mutate(Replicate = rep(2, nrow(metab_rep2bH)),
           Level = rep("High", nrow(metab_rep2bH)))

metab_rep3aH <- metab_rep3aH %>%
    mutate(Replicate = rep(3, nrow(metab_rep3aH)),
           Level = rep("High", nrow(metab_rep3aH)))
metab_rep3bH <- metab_rep3bH %>%
    mutate(Replicate = rep(3, nrow(metab_rep3bH)),
           Level = rep("High", nrow(metab_rep3bH)))

# put all data together -> since colnames are the same, do rbind
metab_alldata <- rbind(metab_rep1aL, metab_rep1bL, metab_rep2aL, metab_rep2bL, 
                       metab_rep3aL, metab_rep3bL, metab_rep1aH, metab_rep1bH, 
                       metab_rep2aH, metab_rep2bH, metab_rep3aH, metab_rep3bH)

# check
table(metab_alldata$Replicate, metab_alldata$Row) 
table(metab_alldata$Replicate, metab_alldata$Column) 

# 2. Process data -------------------------------------------------------------

# add condition
metab_alldata_cond <- metab_alldata %>%
    mutate(Condition = case_when(
        (Column == 1 | Column == 2 | Column == 3 | Column == 4) ~ "WT",
        (Column == 5 | Column == 6 | Column == 7 | Column == 8) ~ "KO1",
        (Column == 9 | Column == 10 | Column == 11 | Column == 12) ~ "KO2"),
        Metabolite = case_when(
            Row == 2 ~ "DMSO",
            Row == 3 ~ "Formate",
            Row == 4 ~ "Folate",
            Row == 5 ~ "meTHF",
            Row == 6 ~ "SAM")
        ) %>% 
    filter(Metabolite %in% c("DMSO", "Formate", "Folate", "meTHF", "SAM"))

# check
table(metab_alldata_cond$Replicate, metab_alldata_cond$Row)
table(metab_alldata_cond$Replicate, metab_alldata_cond$Column)
table(metab_alldata_cond$Row, metab_alldata_cond$Metabolite)
table(metab_alldata_cond$Column, metab_alldata_cond$Condition)
table(metab_alldata_cond$Field, metab_alldata_cond$Condition)
table(metab_alldata_cond$Field, metab_alldata_cond$Metabolite)

# filter single nuclei
metab_alldata_singl <- metab_alldata_cond %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200)

nrow(metab_alldata_cond) 
nrow(metab_alldata_singl) 
table(metab_alldata_singl$Condition)

# get interesting data and log intensity
h3pser_alldata <- metab_alldata_singl %>%
    select(Condition, Metabolite, Level, Replicate, All.Nuclei.Selected...H3phosphoSer10.Mean) %>%
    mutate(log2H3pSer = log2(All.Nuclei.Selected...H3phosphoSer10.Mean),
           log10H3pSer = log10(All.Nuclei.Selected...H3phosphoSer10.Mean))

# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_H3pSer_data <- function(n, c) {
    repdata <- h3pser_alldata %>%
        filter(Replicate == n & Level == c) 
    repdata_scale <- repdata %>%
        mutate(H3pSer_sc = range01(All.Nuclei.Selected...H3phosphoSer10.Mean),
               log2H3pSer_sc = range01(log2H3pSer))
}

reps <- c(1,2,3,1,2,3)
levels <- c("Low", "Low", "Low", "High", "High", "High")
h3pser_alldata_scaled_reps <- map2(reps, levels, scale_H3pSer_data)

# build data frame
h3pser_alldata_scaled <- bind_rows(h3pser_alldata_scaled_reps)
table(h3pser_alldata_scaled$Condition)

# plot again histograms
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
    xlab("scaled H3pSer10 mean intensity") +
    theme_classic()
ggsave("plots/histogram_high_values_scaled.pdf")

# 3. Get % of mitotic cels ------------------------------------------------------------

# get mitotic cells with threshold
mitosis_quant_thresmit <- h3pser_alldata_scaled %>%
    mutate(Mitotic = ifelse(H3pSer_sc > 0.16, T, F))

table(mitosis_quant_thresmit$Condition, mitosis_quant_thresmit$Mitotic)

# convert condition, level and metabolite to factor
mitosis_quant_thresmit$Condition <- factor(mitosis_quant_thresmit$Condition, levels = c("WT", "KO1", "KO2"))
mitosis_quant_thresmit$Metabolite <- factor(mitosis_quant_thresmit$Metabolite, levels = c("DMSO", "Formate", "Folate", "meTHF", "SAM"))
mitosis_quant_thresmit$Level <- factor(mitosis_quant_thresmit$Level, levels = c("Low", "High"))

# average and sd for each treatment
mitosis_quant_avg_sc <- mitosis_quant_thresmit %>%
    group_by(Condition, Metabolite, Level, Replicate) %>%
    summarise(
        count = n(),
        mean = mean(H3pSer_sc),
        median = median(H3pSer_sc),
        sd = sd(H3pSer_sc),
        se = sd/sqrt(count),
        perc_mitotic = sum(Mitotic)/count*100)

mitosis_quant_avg_sc$Condition <- factor(mitosis_quant_avg_sc$Condition, levels = c("WT", "KO1", "KO2"))
mitosis_quant_avg_sc$Metabolite <- factor(mitosis_quant_avg_sc$Metabolite, levels = c("DMSO", "Formate", "Folate", "meTHF", "SAM"))
mitosis_quant_avg_sc$Level <- factor(mitosis_quant_avg_sc$Level, levels = c("Low", "High"))

# decide threshold
mean(unlist(mitosis_quant_avg_sc[mitosis_quant_avg_sc$Condition == "WT" & mitosis_quant_avg_sc$Metabolite == "DMSO" 
                                 & mitosis_quant_avg_sc$Level == "Low", "median"])) # 0.002528056
mean(unlist(mitosis_quant_avg_sc[mitosis_quant_avg_sc$Condition == "WT" & mitosis_quant_avg_sc$Metabolite == "DMSO" 
                                 & mitosis_quant_avg_sc$Level == "Low", "sd"])) # 0.04288337

mean(unlist(mitosis_quant_avg_sc[mitosis_quant_avg_sc$Condition == "WT" & mitosis_quant_avg_sc$Metabolite == "DMSO" & 
                                     mitosis_quant_avg_sc$Level == "High", "median"])) # 0.002434639
mean(unlist(mitosis_quant_avg_sc[mitosis_quant_avg_sc$Condition == "WT" & mitosis_quant_avg_sc$Metabolite == "DMSO" & 
                                     mitosis_quant_avg_sc$Level == "High", "sd"])) # 0.0522789

# plot
ggplot(mitosis_quant_thresmit[mitosis_quant_thresmit$Level == "Low", ], aes(x = Metabolite, y = H3pSer_sc, color = Mitotic)) +
    geom_jitter(alpha = 0.1) +
    geom_hline(yintercept=0.16) +
    facet_grid(rows = vars(Condition), cols = vars(Replicate)) +
    scale_color_manual(values = c("TRUE" = "#559A73", "FALSE" = "lightgrey")) +
    ylab("Scaled H3pSer10 mean intensity") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle=90, hjust=1), axis.title = element_text(size = 14),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 12))
ggsave("plots/H3pS3_quantif_all_reps_scaled_low.pdf", width = 7, height = 4)

ggplot(mitosis_quant_thresmit[mitosis_quant_thresmit$Level == "High", ], aes(x = Metabolite, y = H3pSer_sc, color = Mitotic)) +
    geom_jitter(alpha = 0.1) +
    geom_hline(yintercept=0.16) +
    facet_grid(rows = vars(Condition), cols = vars(Replicate)) +
    scale_color_manual(values = c("TRUE" = "#559A73", "FALSE" = "lightgrey")) +
    ylab("Scaled H3pSer10 mean intensity") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle=90, hjust=1), axis.title = element_text(size = 14),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 12))
ggsave("plots/H3pS3_quantif_all_reps_scaled_high.pdf", width = 7, height = 4)

# 4. Statistical analysis -----------------------------------------------------------------------------

# data of perc_mitotic
mitotic_perc_info_sc <- as.data.frame(mitosis_quant_avg_sc[, c("Condition", "Metabolite", "Replicate", "Level", "perc_mitotic")])

# statistic summary
mitotic_perc_stats_sc <- mitotic_perc_info_sc %>%
    group_by(Condition, Metabolite, Level) %>%
    summarise(
        count = n(),
        mean_mit = mean(perc_mitotic),
        median_mit= median(perc_mitotic),
        sd_mit = sd(perc_mitotic),
        se_mit = sd_mit/sqrt(count))

# check normality of groups -> WT
group_WT_DMSO_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "WT" & Metabolite == "DMSO")
group_WT_form_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "WT" & Metabolite == "Formate")
group_WT_fol_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "WT" & Metabolite == "Folate")
group_WT_THF_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "WT" & Metabolite == "meTHF")
group_WT_SAM_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "WT" & Metabolite == "SAM")
shapiro.test(group_WT_DMSO_L$perc_mitotic)
shapiro.test(group_WT_form_L$perc_mitotic)
shapiro.test(group_WT_fol_L$perc_mitotic)
shapiro.test(group_WT_THF_L$perc_mitotic)
shapiro.test(group_WT_SAM_L$perc_mitotic)

# check homogeneity of variance -> yes
leveneTest(perc_mitotic ~ Metabolite, data=mitotic_perc_info_sc[mitotic_perc_info_sc$Condition == "WT"
                                                                & mitotic_perc_info_sc$Level == "Low",])

# t-test between metabolites and DMSO
group_WT_DMSO_FormL <- mitotic_perc_info_sc %>%
    filter(Condition == "WT" & Level == "Low" & (Metabolite == "DMSO" | Metabolite == "Formate"))
group_WT_DMSO_FolL <- mitotic_perc_info_sc %>%
    filter(Condition == "WT" & Level == "Low" & (Metabolite == "DMSO" | Metabolite == "Folate"))
group_WT_DMSO_THFL <- mitotic_perc_info_sc %>%
    filter(Condition == "WT" & Level == "Low" & (Metabolite == "DMSO" | Metabolite == "meTHF"))
group_WT_DMSO_SAML <- mitotic_perc_info_sc %>%
    filter(Condition == "WT" & Level == "Low" & (Metabolite == "DMSO" | Metabolite == "SAM"))

t.test(perc_mitotic ~ Metabolite, data=group_WT_DMSO_FormL, var.equal = TRUE) 
t.test(perc_mitotic ~ Metabolite, data=group_WT_DMSO_FolL, var.equal = TRUE) 
t.test(perc_mitotic ~ Metabolite, data=group_WT_DMSO_THFL, var.equal = TRUE)
t.test(perc_mitotic ~ Metabolite, data=group_WT_DMSO_SAML, var.equal = TRUE) 

# check normality of groups -> KO1
group_K1_DMSO_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO1" & Metabolite == "DMSO")
group_K1_form_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO1" & Metabolite == "Formate")
group_K1_fol_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO1" & Metabolite == "Folate")
group_K1_THF_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO1" & Metabolite == "meTHF")
group_K1_SAM_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO1" & Metabolite == "SAM")
shapiro.test(group_K1_DMSO_L$perc_mitotic)
shapiro.test(group_K1_form_L$perc_mitotic)
shapiro.test(group_K1_fol_L$perc_mitotic)
shapiro.test(group_K1_THF_L$perc_mitotic)
shapiro.test(group_K1_SAM_L$perc_mitotic)

# check homogeneity of variance -> yes
leveneTest(perc_mitotic ~ Metabolite, data=mitotic_perc_info_sc[mitotic_perc_info_sc$Condition == "KO1" & mitotic_perc_info_sc$Level == "Low",])

# t-test between metabolites and DMSO
group_K1_DMSO_FormL <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO1" & (Metabolite == "DMSO" | Metabolite == "Formate"))
group_K1_DMSO_FolL <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO1" & (Metabolite == "DMSO" | Metabolite == "Folate"))
group_K1_DMSO_THFL <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO1" & (Metabolite == "DMSO" | Metabolite == "meTHF"))
group_K1_DMSO_SAML <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO1" & (Metabolite == "DMSO" | Metabolite == "SAM"))

t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_FormL, var.equal = TRUE) # no
t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_FolL, var.equal = TRUE) # no
t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_THFL, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_SAML, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_thyL, var.equal = TRUE) # no 

# check normality of groups -> KO2
group_K2_DMSO_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO2" & Metabolite == "DMSO")
group_K2_form_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO2" & Metabolite == "Formate")
group_K2_fol_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO2" & Metabolite == "Folate")
group_K2_THF_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO2" & Metabolite == "meTHF")
group_K2_SAM_L <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO2" & Metabolite == "SAM")
shapiro.test(group_K2_DMSO_L$perc_mitotic)
shapiro.test(group_K2_form_L$perc_mitotic) # not normal but close
shapiro.test(group_K2_fol_L$perc_mitotic)
shapiro.test(group_K2_THF_L$perc_mitotic)
shapiro.test(group_K2_SAM_L$perc_mitotic)

# check homogeneity of variance -> yes
leveneTest(perc_mitotic ~ Metabolite, data=mitotic_perc_info_sc[mitotic_perc_info_sc$Condition == "KO2" & mitotic_perc_info_sc$Level == "Low",])

# t-test between metabolites and DMSO
group_K2_DMSO_FormL <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO2" & (Metabolite == "DMSO" | Metabolite == "Formate"))
group_K2_DMSO_FolL <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO2" & (Metabolite == "DMSO" | Metabolite == "Folate"))
group_K2_DMSO_THFL <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO2" & (Metabolite == "DMSO" | Metabolite == "meTHF"))
group_K2_DMSO_SAML <- mitotic_perc_info_sc %>%
    filter(Level == "Low" & Condition == "KO2" & (Metabolite == "DMSO" | Metabolite == "SAM"))
t.test(perc_mitotic ~ Metabolite, data=group_K2_DMSO_FormL, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K2_DMSO_FolL, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K2_DMSO_THFL, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K2_DMSO_SAML, var.equal = TRUE) # no 

# high

# check normality of groups -> WT
group_WT_DMSO_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "WT" & Metabolite == "DMSO")
group_WT_form_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "WT" & Metabolite == "Formate")
group_WT_fol_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "WT" & Metabolite == "Folate")
group_WT_THF_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "WT" & Metabolite == "meTHF")
group_WT_SAM_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "WT" & Metabolite == "SAM")
shapiro.test(group_WT_DMSO_H$perc_mitotic)
shapiro.test(group_WT_form_H$perc_mitotic)
shapiro.test(group_WT_fol_H$perc_mitotic)
shapiro.test(group_WT_THF_H$perc_mitotic)
shapiro.test(group_WT_SAM_H$perc_mitotic)

# check homogeneity of variance -> yes
leveneTest(perc_mitotic ~ Metabolite, data=mitotic_perc_info_sc[mitotic_perc_info_sc$Condition == "WT"
                                                                & mitotic_perc_info_sc$Level == "High",])

# t-test between metabolites and DMSO
group_WT_DMSO_FormH <- mitotic_perc_info_sc %>%
    filter(Condition == "WT" & Level == "High" & (Metabolite == "DMSO" | Metabolite == "Formate"))
group_WT_DMSO_FolH <- mitotic_perc_info_sc %>%
    filter(Condition == "WT" & Level == "High" & (Metabolite == "DMSO" | Metabolite == "Folate"))
group_WT_DMSO_THFH <- mitotic_perc_info_sc %>%
    filter(Condition == "WT" & Level == "High" & (Metabolite == "DMSO" | Metabolite == "meTHF"))
group_WT_DMSO_SAMH <- mitotic_perc_info_sc %>%
    filter(Condition == "WT" & Level == "High" & (Metabolite == "DMSO" | Metabolite == "SAM"))

t.test(perc_mitotic ~ Metabolite, data=group_WT_DMSO_FormH, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_WT_DMSO_FolH, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_WT_DMSO_THFH, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_WT_DMSO_SAMH, var.equal = TRUE) # no 

# check normality of groups -> KO1
group_K1_DMSO_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO1" & Metabolite == "DMSO")
group_K1_form_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO1" & Metabolite == "Formate")
group_K1_fol_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO1" & Metabolite == "Folate")
group_K1_THF_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO1" & Metabolite == "meTHF")
group_K1_SAM_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO1" & Metabolite == "SAM")
shapiro.test(group_K1_DMSO_H$perc_mitotic)
shapiro.test(group_K1_form_H$perc_mitotic)
shapiro.test(group_K1_fol_H$perc_mitotic)
shapiro.test(group_K1_THF_H$perc_mitotic)
shapiro.test(group_K1_SAM_H$perc_mitotic) 

# check homogeneity of variance -> yes
leveneTest(perc_mitotic ~ Metabolite, data=mitotic_perc_info_sc[mitotic_perc_info_sc$Condition == "KO1" & mitotic_perc_info_sc$Level == "High",])

# t-test between metabolites and DMSO
group_K1_DMSO_FormH <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO1" & (Metabolite == "DMSO" | Metabolite == "Formate"))
group_K1_DMSO_FolH <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO1" & (Metabolite == "DMSO" | Metabolite == "Folate"))
group_K1_DMSO_THFH <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO1" & (Metabolite == "DMSO" | Metabolite == "meTHF"))
group_K1_DMSO_SAMH <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO1" & (Metabolite == "DMSO" | Metabolite == "SAM"))

t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_FormH, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_FolH, var.equal = TRUE) # no
t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_THFH, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_SAMH, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K1_DMSO_thyH, var.equal = TRUE) # no 

# check normality of groups -> KO2
group_K2_DMSO_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO2" & Metabolite == "DMSO")
group_K2_form_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO2" & Metabolite == "Formate")
group_K2_fol_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO2" & Metabolite == "Folate")
group_K2_THF_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO2" & Metabolite == "meTHF")
group_K2_SAM_H <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO2" & Metabolite == "SAM")
shapiro.test(group_K2_DMSO_H$perc_mitotic)
shapiro.test(group_K2_form_H$perc_mitotic) 
shapiro.test(group_K2_fol_H$perc_mitotic)
shapiro.test(group_K2_THF_H$perc_mitotic)
shapiro.test(group_K2_SAM_H$perc_mitotic)

# check homogeneity of variance -> yes
leveneTest(perc_mitotic ~ Metabolite, data=mitotic_perc_info_sc[mitotic_perc_info_sc$Condition == "KO2" & mitotic_perc_info_sc$Level == "High",])

# t-test between metabolites and DMSO
group_K2_DMSO_FormH <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO2" & (Metabolite == "DMSO" | Metabolite == "Formate"))
group_K2_DMSO_FolH <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO2" & (Metabolite == "DMSO" | Metabolite == "Folate"))
group_K2_DMSO_THFH <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO2" & (Metabolite == "DMSO" | Metabolite == "meTHF"))
group_K2_DMSO_SAMH <- mitotic_perc_info_sc %>%
    filter(Level == "High" & Condition == "KO2" & (Metabolite == "DMSO" | Metabolite == "SAM"))

t.test(perc_mitotic ~ Metabolite, data=group_K2_DMSO_FormH, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K2_DMSO_FolH, var.equal = TRUE) # no 
t.test(perc_mitotic ~ Metabolite, data=group_K2_DMSO_THFH, var.equal = TRUE) # yes
t.test(perc_mitotic ~ Metabolite, data=group_K2_DMSO_SAMH, var.equal = TRUE) # yes

# 5. Plot normalized with DMSO --------------------------------------------------------------

# normalized with DMSO
WT_DMSO_low_1 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "WT" &
                                                       mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                       mitotic_perc_info_sc$Level == "Low" &
                                                       mitotic_perc_info_sc$Replicate == 1]
WT_DMSO_low_2 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "WT" &
                                                       mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                       mitotic_perc_info_sc$Level == "Low"&
                                                       mitotic_perc_info_sc$Replicate == 2]
WT_DMSO_low_3 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "WT" &
                                                       mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                       mitotic_perc_info_sc$Level == "Low"&
                                                       mitotic_perc_info_sc$Replicate == 3]

KO1_DMSO_low_1 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO1" &
                                                        mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                        mitotic_perc_info_sc$Level == "Low" &
                                                        mitotic_perc_info_sc$Replicate == 1]
KO1_DMSO_low_2 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO1" &
                                                        mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                        mitotic_perc_info_sc$Level == "Low"&
                                                        mitotic_perc_info_sc$Replicate == 2]
KO1_DMSO_low_3 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO1" &
                                                        mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                        mitotic_perc_info_sc$Level == "Low"&
                                                        mitotic_perc_info_sc$Replicate == 3]

KO2_DMSO_low_1 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO2" &
                                                        mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                        mitotic_perc_info_sc$Level == "Low" &
                                                        mitotic_perc_info_sc$Replicate == 1]
KO2_DMSO_low_2 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO2" &
                                                        mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                        mitotic_perc_info_sc$Level == "Low"&
                                                        mitotic_perc_info_sc$Replicate == 2]
KO2_DMSO_low_3 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO2" &
                                                        mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                        mitotic_perc_info_sc$Level == "Low"&
                                                        mitotic_perc_info_sc$Replicate == 3]

WT_DMSO_high_1 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "WT" &
                                                        mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                        mitotic_perc_info_sc$Level == "High" &
                                                        mitotic_perc_info_sc$Replicate == 1]
WT_DMSO_high_2 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "WT" &
                                                        mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                        mitotic_perc_info_sc$Level == "High"&
                                                        mitotic_perc_info_sc$Replicate == 2]
WT_DMSO_high_3 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "WT" &
                                                        mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                        mitotic_perc_info_sc$Level == "High"&
                                                        mitotic_perc_info_sc$Replicate == 3]

KO1_DMSO_high_1 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO1" &
                                                         mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                         mitotic_perc_info_sc$Level == "High" &
                                                         mitotic_perc_info_sc$Replicate == 1]
KO1_DMSO_high_2 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO1" &
                                                         mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                         mitotic_perc_info_sc$Level == "High"&
                                                         mitotic_perc_info_sc$Replicate == 2]
KO1_DMSO_high_3 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO1" &
                                                         mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                         mitotic_perc_info_sc$Level == "High"&
                                                         mitotic_perc_info_sc$Replicate == 3]

KO2_DMSO_high_1 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO2" &
                                                         mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                         mitotic_perc_info_sc$Level == "High" &
                                                         mitotic_perc_info_sc$Replicate == 1]
KO2_DMSO_high_2 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO2" &
                                                         mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                         mitotic_perc_info_sc$Level == "High"&
                                                         mitotic_perc_info_sc$Replicate == 2]
KO2_DMSO_high_3 <- mitotic_perc_info_sc$perc_mitotic[mitotic_perc_info_sc$Condition == "KO2" &
                                                         mitotic_perc_info_sc$Metabolite == "DMSO" &
                                                         mitotic_perc_info_sc$Level == "High"&
                                                         mitotic_perc_info_sc$Replicate == 3]

mitotic_perc_info_sc_norm_byrep <- mitotic_perc_info_sc %>%
    mutate(perc_mitotic_norm = case_when(
        Condition == "WT" & Level == "Low" & Replicate == 1 ~ perc_mitotic/WT_DMSO_low_1,
        Condition == "WT" & Level == "Low" & Replicate == 2 ~ perc_mitotic/WT_DMSO_low_2,
        Condition == "WT" & Level == "Low" & Replicate == 3 ~ perc_mitotic/WT_DMSO_low_3,
        Condition == "KO1" & Level == "Low" & Replicate == 1 ~ perc_mitotic/KO1_DMSO_low_1,
        Condition == "KO1" & Level == "Low" & Replicate == 2 ~ perc_mitotic/KO1_DMSO_low_2,
        Condition == "KO1" & Level == "Low" & Replicate == 3 ~ perc_mitotic/KO1_DMSO_low_3,
        Condition == "KO2" & Level == "Low" & Replicate == 1 ~ perc_mitotic/KO2_DMSO_low_1,
        Condition == "KO2" & Level == "Low" & Replicate == 2 ~ perc_mitotic/KO2_DMSO_low_2,
        Condition == "KO2" & Level == "Low" & Replicate == 3 ~ perc_mitotic/KO2_DMSO_low_3,
        
        Condition == "WT" & Level == "High" & Replicate == 1 ~ perc_mitotic/WT_DMSO_high_1,
        Condition == "WT" & Level == "High" & Replicate == 2 ~ perc_mitotic/WT_DMSO_high_2,
        Condition == "WT" & Level == "High" & Replicate == 3 ~ perc_mitotic/WT_DMSO_high_3,
        Condition == "KO1" & Level == "High" & Replicate == 1 ~ perc_mitotic/KO1_DMSO_high_1,
        Condition == "KO1" & Level == "High" & Replicate == 2 ~ perc_mitotic/KO1_DMSO_high_2,
        Condition == "KO1" & Level == "High" & Replicate == 3 ~ perc_mitotic/KO1_DMSO_high_3,
        Condition == "KO2" & Level == "High" & Replicate == 1 ~ perc_mitotic/KO2_DMSO_high_1,
        Condition == "KO2" & Level == "High" & Replicate == 2 ~ perc_mitotic/KO2_DMSO_high_2,
        Condition == "KO2" & Level == "High" & Replicate == 3 ~ perc_mitotic/KO2_DMSO_high_3),
        MetaboliteLevel = ifelse(Level == "Low", paste0(Metabolite,"L"), paste0(Metabolite,"H")))

# statistic summary
mitotic_perc_stats_sc_norm_byrep <- mitotic_perc_info_sc_norm_byrep %>%
    group_by(Condition, MetaboliteLevel) %>%
    summarise(
        count = n(),
        mean_mit = mean(perc_mitotic_norm),
        median_mit= median(perc_mitotic_norm),
        sd_mit = sd(perc_mitotic_norm),
        se_mit = sd_mit/sqrt(count))

# ggplot data
my_comparisons <- list(c("DMSOL", "FormateL"), c("DMSOL", "FolateL"), c("DMSOL", "meTHFL"),
                       c("DMSOL", "SAML"),
                       c("DMSOH", "FormateH"), c("DMSOH", "FolateH"), c("DMSOH", "meTHFH"),
                       c("DMSOH", "SAMH"))

mitotic_perc_stats_sc_norm_byrep$MetaboliteLevel <- factor(mitotic_perc_stats_sc_norm_byrep$MetaboliteLevel, 
                                                           levels = c("DMSOL", "DMSOH", "FormateL", "FormateH",
                                                                      "FolateL", "FolateH", "meTHFL", "meTHFH",
                                                                      "SAML", "SAMH"))
mitotic_perc_info_sc_norm_byrep$MetaboliteLevel <- factor(mitotic_perc_info_sc_norm_byrep$MetaboliteLevel, 
                                                          levels = c("DMSOL", "DMSOH", "FormateL", "FormateH",
                                                                     "FolateL", "FolateH", "meTHFL", "meTHFH",
                                                                     "SAML", "SAMH"))

ggplot() +
    geom_errorbar(data = mitotic_perc_stats_sc_norm_byrep, aes(x=MetaboliteLevel, ymin=mean_mit-sd_mit, 
                                                               ymax=mean_mit+sd_mit), 
                  width=.2,
                  position=position_dodge(0.05)) +
    geom_col(data = mitotic_perc_stats_sc_norm_byrep, aes(x=MetaboliteLevel, y=mean_mit, fill = Condition), 
             width = 0.5)+
    geom_jitter(data = mitotic_perc_info_sc_norm_byrep, aes(x=MetaboliteLevel, y=perc_mitotic_norm), 
                height = 0, width = 0.1, alpha = 1, size = 1.5)+
    stat_compare_means(comparisons = my_comparisons, method = "t.test") +
    facet_wrap(~Condition, nrow = 3) +
    xlab("") +
    ylab("Percentage of mitotic cells (%) \n normalized to DMSO") +
    scale_y_continuous(expand=c(0,0), limits = c(0,2)) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58", "KO2" = "#f9beb3")) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("plots/perc_mitotic_cells_bar_pval_scaled_low_high_byrep.pdf", width = 6, height = 5)

# statistics one sample t test low

# t-test low WT
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "WT" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Formate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "WT" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Folate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "WT" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "meTHF"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "WT" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "SAM"], 
       mu = 1, alternative = "two.sided") 


# t-test low KO1
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO1" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Formate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO1" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Folate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO1" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "meTHF"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO1" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "SAM"], 
       mu = 1, alternative = "two.sided") 

# t-test low KO2
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO2" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Formate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO2" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Folate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO2" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "meTHF"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO2" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "Low" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "SAM"], 
       mu = 1, alternative = "two.sided") 

# t-test high WT
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "WT" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Formate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "WT" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Folate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "WT" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "meTHF"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "WT" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "SAM"], 
       mu = 1, alternative = "two.sided") 

# t-test high KO1
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO1" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Formate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO1" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Folate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO1" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "meTHF"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO1" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "SAM"], 
       mu = 1, alternative = "two.sided") 

# t-test high KO2
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO2" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Formate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO2" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "Folate"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO2" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "meTHF"], 
       mu = 1, alternative = "two.sided") 
t.test(mitotic_perc_info_sc_norm_byrep$perc_mitotic_norm[mitotic_perc_info_sc_norm_byrep$Condition == "KO2" 
                                                         & mitotic_perc_info_sc_norm_byrep$Level == "High" 
                                                         & mitotic_perc_info_sc_norm_byrep$Metabolite == "SAM"], 
       mu = 1, alternative = "two.sided") 
