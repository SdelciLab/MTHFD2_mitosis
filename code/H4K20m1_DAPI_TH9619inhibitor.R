# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("nortest")
library("ggrepel")

# Objective: Analyse H4K20me1 and DAPI levels of HCT116 treated with Th9619 inhibitor

# 1. Get input data

# get data
H4K20_rep1 <- read.delim("rep2/H4K20me1 HAP1 HCT MTHinh 230223__2023-02-23T12_09_57-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep1a)

H4K20_rep2 <- read.delim("rep3/H4K20me1 HCT116 MTHinh 010323__2023-03-01T17_31_02-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep2b)

H4K20_rep3 <- read.delim("rep4/H4K20me1 HCT116 MTHinh 070323__2023-03-07T19_41_52-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep3b)

# check colnames are the same
table(colnames(H4K20_rep1) == colnames(H4K20_rep2))
table(colnames(H4K20_rep1) == colnames(H4K20_rep3))

# add replicate
H4K20_rep1 <- H4K20_rep1 %>%
    mutate(Replicate = rep(1, nrow(H4K20_rep1)),
           Plate = rep("HAP1_HCT116", nrow(H4K20_rep1)))
H4K20_rep2 <- H4K20_rep2 %>%
    mutate(Replicate = rep(2, nrow(H4K20_rep2)),
           Plate = rep("HCT116", nrow(H4K20_rep2)))
H4K20_rep3 <- H4K20_rep3 %>%
    mutate(Replicate = rep(3, nrow(H4K20_rep3)),
           Plate = rep("HCT116", nrow(H4K20_rep3)))

# put all data together -> since colnames are the same, do rbind
H4K20_alldata <- rbind(H4K20_rep1, H4K20_rep2, H4K20_rep3)

# check
table(H4K20_alldata$Replicate, H4K20_alldata$Row) 
table(H4K20_alldata$Replicate, H4K20_alldata$Column)  

# 2. Process data -------------------------------------------------------------------------------

# add condition
H4K20_alldata_int <- H4K20_alldata %>%
    mutate(Cell_Line = case_when(
        (Row == 4 | Row == 5) & Plate == "HCT116" ~ "HCT116",
        (Row == 2 | Row == 3) & Plate == "HAP1_HCT116" ~ "HAP1-WT",
        (Row == 4 | Row == 5) & Plate == "HAP1_HCT116" ~ "HAP1-mut",
        (Row == 6 | Row == 7) & Plate == "HAP1_HCT116" ~ "HCT116"),
        Drug = case_when(
            (Column == 2) ~ 0,
            (Column == 4 & Replicate != 1) ~ 10,
            (Column == 4 & Replicate == 1) ~ 25,
            (Column == 6 & Replicate != 1) ~ 25,
            (Column == 6 & Replicate == 1) ~ 63,
            (Column == 8 & Replicate != 1) ~ 63,
            (Column == 8 & Replicate == 1) ~ 156,
            (Column == 10 & Replicate != 1) ~ 156,
            (Column == 10 & Replicate == 1) ~ 391),
        log2H4K20me1 = log2(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.Mean),
        IntDen = All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.Mean * All.Nuclei.Selected...Area.Nuclei.Area..µm..,
        log2IntDen = log2(IntDen),
        log2CV = log2(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.CV....),
        logDAPI = log(All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean),
        log2DAPI = log2(All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean),
        log2CVDAPI = log2(All.Nuclei.Selected...Intensity.Nucleus.DAPI.CV....)) %>% 
    filter(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.CV.... != 0)

# check
table(H4K20_alldata_int$Replicate, H4K20_alldata_int$Cell_Line) 
table(H4K20_alldata_int$Row, H4K20_alldata_int$Cell_Line) 
table(H4K20_alldata_int$Replicate, H4K20_alldata_int$Drug)  
table(H4K20_alldata_int$Column, H4K20_alldata_int$Drug)  

# remove first and last concentration because no 3 replicates
H4K20_alldata_int2 <- H4K20_alldata_int %>%
    filter(Drug != 391 & Drug != 10 )
table(H4K20_alldata_int2$Column, H4K20_alldata_int2$Drug)  

# change colname
colnames(H4K20_alldata_int2)[18] <- "H4K20m1_CV"

# filter single nuclei and round
H4K20_alldata_singl <- H4K20_alldata_int2 %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200 & `All.Nuclei.Selected...Area.Nuclei.Area..µm..` > 50,
           `All.Nuclei.Selected...Area.Nuclei.Roundness` > 0.8)
nrow(H4K20_alldata_singl) 
nrow(H4K20_alldata_int2) 
table(H4K20_alldata_singl$Cell_Line) 

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
table(H4K20_alldata_scaled_DAPI$Cell_Line)

# mitotic cells have higher DAPI signal
hist(H4K20_alldata_scaled_DAPI$All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)
hist(H4K20_alldata_scaled_DAPI$logDAPI)
hist(H4K20_alldata_scaled_DAPI$logDAPI_sc)
# we could try 0.7 as threshold

# remove mitotic cells
H4K20_alldata_scaled_DAPI_nomit <- H4K20_alldata_scaled_DAPI %>%
    filter(logDAPI_sc < 0.7)
nrow(H4K20_alldata_scaled_DAPI_nomit) 
nrow(H4K20_alldata_scaled_DAPI) 

# export data
save(H4K20_alldata_scaled_DAPI_nomit, file="H4K20_alldata_nomit.RData")

# 3. Statistics and plots --------------------------------------------------------------------

# load data
load("H4K20_alldata_nomit.RData")

# Analysis of HCT116
H4K20_alldata_HCT116 <- H4K20_alldata_scaled_DAPI_nomit %>%
    filter(Cell_Line == "HCT116")
H4K20_alldata_HCT116$Drug <- factor(H4K20_alldata_HCT116$Drug)
H4K20_alldata_HCT116$Replicate <- factor(H4K20_alldata_HCT116$Replicate)

# add facet
H4K20_alldata_HCT116$Mark <- rep("H4K20m1", nrow(H4K20_alldata_HCT116)) 
H4K20_alldata_HCT116$Mark <- factor(H4K20_alldata_HCT116$Mark)

# plot 
my_comp_drug <- list(c("0", "25"), c("0", "63"), c("0", "156"), c("25", "63"), c("25", "156"), c("63", "156"))

# plot with violin mean intensity
ggplot(H4K20_alldata_HCT116, aes(x = Drug, y = log2H4K20me1)) +
    geom_violin(fill = "#aaafb0",  col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    stat_compare_means(comparisons = my_comp_drug, size = 2, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "TH9619 Inhibitor (nM)", y = "log2(Mean Intensity)", fill = "") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("reps_plots/log2MeanIntH4K20m1_allreps_violin.pdf", width = 4.5, height = 3)

# plot with violin coefficient of variation
ggplot(H4K20_alldata_HCT116, aes(x = Drug, y = log2CV)) +
    geom_violin(fill = "#aaafb0", col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    stat_compare_means(comparisons = my_comp_drug, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "TH9619 Inhibitor (nM)", y = "log2(Coefficient of Variation)", fill = "") +
    theme_classic() +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("reps_plots/log2CVallreps_violin.pdf", width = 4.5, height = 3)

# plot DAPI

# add facet
H4K20_alldata_HCT116$Mark <- rep("DAPI", nrow(H4K20_alldata_HCT116)) 
H4K20_alldata_HCT116$Mark <- factor(H4K20_alldata_HCT116$Mark)

# plot log2DAPI
ggplot(H4K20_alldata_HCT116, aes(x = Drug, y = log2DAPI)) +
    geom_violin(fill = "#aaafb0",  col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    stat_compare_means(comparisons = my_comp_drug, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "TH9619 Inhibitor (nM)", y = "log2(Mean Intensity)", fill = "") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("reps_plots/log2DAPI_allreps_violin.pdf", width = 4.5, height = 3)

# plot CV DAPI
ggplot(H4K20_alldata_HCT116, aes(x = Drug, y = log2CVDAPI)) +
    geom_violin(fill = "#aaafb0", col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    stat_compare_means(comparisons = my_comp_drug, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "TH9619 Inhibitor (nM)", y = "log2(Coefficient of Variation)", fill = "") +
    theme_classic() +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("reps_plots/log2CVDAPIallreps_violin.pdf", width = 4.5, height = 3)

