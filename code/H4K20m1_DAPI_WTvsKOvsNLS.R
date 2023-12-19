# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("nortest")
library("ggrepel")

# Objective: Analyse H4K20m1 and DAPI between MTHFD2 WT, KO and NLS cells

# 1. Get input data ----------------------------------------------------

# get data
H4K20_rep1 <- read.delim("rep1/H4K20me1 WTvsNLS 40X 010323__2023-03-01T18_16_52-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep1)

H4K20_rep2 <- read.delim("rep2 old/H4K20me1 WTvsNLS 40X old 010323__2023-03-01T18_26_11-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep2)

H4K20_rep3 <- read.delim("rep3/H4K20me1 WTvsNLS 40X 080323__2023-03-08T15_57_38-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep3)

# check colnames are the same
table(colnames(H4K20_rep1) == colnames(H4K20_rep2))
table(colnames(H4K20_rep1) == colnames(H4K20_rep3))

# add replicate
H4K20_rep1 <- H4K20_rep1 %>%
    mutate(Replicate = rep(1, nrow(H4K20_rep1)))
H4K20_rep2 <- H4K20_rep2 %>%
    mutate(Replicate = rep(2, nrow(H4K20_rep2)))
H4K20_rep3 <- H4K20_rep3 %>%
    mutate(Replicate = rep(3, nrow(H4K20_rep3)))

# put all data together -> since colnames are the same, do rbind
H4K20_alldata <- rbind(H4K20_rep1, H4K20_rep2, H4K20_rep3)

# check
table(H4K20_alldata$Replicate, H4K20_alldata$Row) 
table(H4K20_alldata$Replicate, H4K20_alldata$Column)  

# 2. Process data ----------------------------------------------

# add condition
H4K20_alldata_int <- H4K20_alldata %>%
    mutate(Condition = case_when(
        Row == 2  ~ "WT",
        Row == 3  ~ "KO1",
        Row == 4  ~ "KO2",
        Row == 5  ~ "NLS"),
        H4K20m1_log = log2(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.Mean),
        log2CV = log2(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.CV....),
        logDAPI = log(All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean),
        log2DAPI = log2(All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean),
        log2CVDAPI = log2(All.Nuclei.Selected...Intensity.Nucleus.DAPI.CV....)) %>%
    filter(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.CV.... != 0)

# check
table(H4K20_alldata_int$Replicate, H4K20_alldata_int$Condition) 
table(H4K20_alldata_int$Row, H4K20_alldata_int$Condition) 
table(H4K20_alldata_int$Column, H4K20_alldata_int$Condition)  

# change colname
colnames(H4K20_alldata_int)[18] <- "H4K20m1_CV"

# filter single nuclei and round
H4K20_alldata_singl <- H4K20_alldata_int %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200 & `All.Nuclei.Selected...Area.Nuclei.Area..µm..` > 60,
           `All.Nuclei.Selected...Area.Nuclei.Roundness` > 0.9)
nrow(H4K20_alldata_singl) # 19191
nrow(H4K20_alldata_int) # 26772
table(H4K20_alldata_singl$Condition) 

# remove mitotic cells with DAPI 
# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_DAPI_data <- function(n) {
    repdata <- H4K20_alldata_singl %>%
        filter(Replicate == n) 
    repdata_scale <- repdata %>%
        mutate(logDAPI_sc = range01(logDAPI))
}

reps <- 1:3
H4K20_data_scaled_DAPI <- map(reps, scale_DAPI_data)

# build data frame
H4K20_alldata_scaled_DAPI <- bind_rows(H4K20_data_scaled_DAPI)
table(H4K20_alldata_scaled_DAPI$Condition)

# mitotic cells have higher DAPI signal
hist(H4K20_alldata_scaled_DAPI$All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)
hist(H4K20_alldata_scaled_DAPI$logDAPI)
hist(H4K20_alldata_scaled_DAPI$logDAPI_sc)
# we could try 0.7 as threshold

# remove mitotic cells
H4K20_alldata_scaled_DAPI_nomit <- H4K20_alldata_scaled_DAPI %>%
    filter(logDAPI_sc < 0.7)
nrow(H4K20_alldata_scaled_DAPI_nomit) # 17448
nrow(H4K20_alldata_scaled_DAPI) # 19191

# export data
save(H4K20_alldata_scaled_DAPI_nomit, file="H4K20_alldata_nomit.RData")

# 3. Statistics and plots ---------------------------------------------------------------

# stats
H4K20_alldata_scaled_DAPI_mit_stas <- H4K20_alldata_scaled_DAPI_mit %>%
    group_by(Condition) %>%
    summarise(mean_log2H4K20m1 = mean(H4K20m1_log),
              median_log2H4K20m1 = median(H4K20m1_log),
              sd_log2H4K20m1 = sd(H4K20m1_log),
              mean_CV_log2H4K20m1 = mean(log2CV),
              median_CV_log2H4K20m1 = median(log2CV),
              sd_CV_log2H4K20m1 = sd(log2CV),
              mean_DAPI = mean(log2DAPI),
              median_DAPI = median(log2DAPI),
              sd_DAPI = sd(log2DAPI),
              mean_CV_DAPI = mean(log2CVDAPI),
              median_CV_DAPI = median(log2CVDAPI),
              sd_CV_DAPI = sd(log2CVDAPI))

# plot H4K20me1

my_comp <- list(c("WT", "KO1"), c("WT", "KO2"), c("WT", "NLS"))

# add H4K20m1 facet
H4K20_alldata_scaled_DAPI_nomit$Mark <- rep("H4K20m1", nrow(H4K20_alldata_scaled_DAPI_nomit)) 
H4K20_alldata_scaled_DAPI_nomit$Mark <- factor(H4K20_alldata_scaled_DAPI_nomit$Mark)

# plot violin
ggplot(H4K20_alldata_scaled_DAPI_nomit, aes(x = Condition, y = H4K20m1_log)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO2" = "#f9beb3", "KO1" = "#e56e58",
                                 "NLS" = "#890808")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(Mean Intensity)", fill = "") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("reps_plots/log2H4K20m1_allreps_violin.pdf",  width = 4.5, height = 3)

ggplot(H4K20_alldata_scaled_DAPI_nomit, aes(x = Condition, y = log2CV)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO2" = "#f9beb3", "KO1" = "#e56e58",
                                 "NLS" = "#890808")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(Coefficient of Variation)", fill = "") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("reps_plots/log2H4K20m1CV_allreps_violin.pdf", width = 4.5, height = 3)

# plot DAPI

# add DAPI facet
H4K20_alldata_scaled_DAPI_nomit$Mark <- rep("DAPI", nrow(H4K20_alldata_scaled_DAPI_nomit)) 
H4K20_alldata_scaled_DAPI_nomit$Mark <- factor(H4K20_alldata_scaled_DAPI_nomit$Mark)

# plot violin
ggplot(H4K20_alldata_scaled_DAPI_nomit, aes(x = Condition, y = log2DAPI)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO2" = "#f9beb3", "KO1" = "#e56e58",
                                 "NLS" = "#890808")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(Mean Intensity)", fill = "") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("reps_plots/DAPI_allreps_violin.pdf",  width = 4.5, height = 3)

ggplot(H4K20_alldata_scaled_DAPI_nomit, aes(x = Condition, y = log2CVDAPI)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO2" = "#f9beb3", "KO1" = "#e56e58",
                                 "NLS" = "#890808")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(Coefficient of Variation)", fill = "") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("reps_plots/DAPICV_allreps_violin.pdf", width = 4.5, height = 3)
