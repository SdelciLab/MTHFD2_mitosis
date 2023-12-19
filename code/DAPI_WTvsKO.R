# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("ggrepel")
library("nortest")

# 1. Get input data ---------------------------------------------------------------------------

# get data H3K9
H3K9_rep1 <- read.delim("analyse_H3K9me3/H3K9me3 atub WT KO KD 40X 030322__2022-03-03T11_06_11-Measurement 1/Evaluation11/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep1)

H3K9_rep2 <- read.delim("analyse_H3K9me3/H3K9me3 atub WT KO KD 40X 070422__2022-04-07T15_19_27-Measurement 1/Evaluation4/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep2)

H3K9_rep3 <- read.delim("analyse_H3K9me3/H3K9me3 WT KO KD 40X 210422__2022-04-21T13_44_19-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep3)

H3K9_rep4 <- read.delim("analyse_H3K9me3/H3K9me3 WT KO KD 40X110522__2022-05-11T14_58_00-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep4)

# check colnames are the same
table(colnames(H3K9_rep1) == colnames(H3K9_rep2))
table(colnames(H3K9_rep1) == colnames(H3K9_rep3))
table(colnames(H3K9_rep1) == colnames(H3K9_rep4))

# add replicate
H3K9_rep1 <- H3K9_rep1 %>%
    mutate(Replicate = rep("rep1", nrow(H3K9_rep1)))
H3K9_rep2 <- H3K9_rep2 %>%
    mutate(Replicate = rep("rep2", nrow(H3K9_rep2)))
H3K9_rep3 <- H3K9_rep3 %>%
    mutate(Replicate = rep("rep3", nrow(H3K9_rep3)))
H3K9_rep4 <- H3K9_rep4 %>%
    mutate(Replicate = rep("rep4", nrow(H3K9_rep4)))

# put all data together -> since colnames are the same, do rbind
H3K9_alldata <- rbind(H3K9_rep1, H3K9_rep2, H3K9_rep3, H3K9_rep4)

# check
table(H3K9_alldata$Replicate, H3K9_alldata$Row) 
table(H3K9_alldata$Replicate, H3K9_alldata$Column)  

# check colnames
colnames(H3K9_alldata)

# select DAPI info
H3K9_alldata_DAPI <- H3K9_alldata %>%
    dplyr::select(Row, Column, Replicate,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.StdDev,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.CV....,
                  `All.Nuclei.Selected...Area.Nuclei.Area..µm..`,
                  `All.Nuclei.Selected...Area.Nuclei.Roundness`) %>%
    mutate(Data = rep("H3K9", nrow(H3K9_alldata)))

colnames(H3K9_alldata_DAPI) <- c("Row", "Column", "Replicate", "DAPI.Mean",
                                 "DAPI.StDev", "DAPI.CV", "Area", "Roundness",
                                 "Data")
head(H3K9_alldata_DAPI)

# get data H3K27me3
H3K27_rep1 <- read.delim("analyse_H3K27me3/H3K27me3 atub WT KO KD 40X 030322__2022-03-03T11_30_02-Measurement 1/Evaluation3/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep1)

H3K27_rep2 <- read.delim("analyse_H3K27me3/H3K27me3 atub WT KO KD 40X 070422__2022-04-07T15_36_17-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep2)

H3K27_rep3 <- read.delim("analyse_H3K27me3/H3K27me3 WT KO KD 40X 210422__2022-04-21T15_22_27-Measurement 1/Evaluation3/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep3)

H3K27_rep4 <- read.delim("analyse_H3K27me3/H3K27me3 WT KO KD 40X 110522__2022-05-11T15_33_25-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep4)

# check colnames are the same
table(colnames(H3K27_rep1) == colnames(H3K27_rep2))
table(colnames(H3K27_rep1) == colnames(H3K27_rep3))
table(colnames(H3K27_rep1) == colnames(H3K27_rep4))

# add replicate
H3K27_rep1 <- H3K27_rep1 %>%
    mutate(Replicate = rep("rep1", nrow(H3K27_rep1)))
H3K27_rep2 <- H3K27_rep2 %>%
    mutate(Replicate = rep("rep2", nrow(H3K27_rep2)))
H3K27_rep3 <- H3K27_rep3 %>%
    mutate(Replicate = rep("rep3", nrow(H3K27_rep3)))
H3K27_rep4 <- H3K27_rep4 %>%
    mutate(Replicate = rep("rep4", nrow(H3K27_rep4)))

# put all data together -> since colnames are the same, do rbind
H3K27_alldata <- rbind(H3K27_rep1, H3K27_rep2, H3K27_rep3, H3K27_rep4)

# check
table(H3K27_alldata$Replicate, H3K27_alldata$Row) 
table(H3K27_alldata$Replicate, H3K27_alldata$Column)  

# select DAPI info
H3K27_alldata_DAPI <- H3K27_alldata %>%
    dplyr::select(Row, Column, Replicate,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.StdDev,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.CV....,
                  `All.Nuclei.Selected...Area.Nuclei.Area..µm..`,
                  `All.Nuclei.Selected...Area.Nuclei.Roundness`) %>%
    mutate(Data = rep("H3K27", nrow(H3K27_alldata)))

colnames(H3K27_alldata_DAPI) <- c("Row", "Column", "Replicate", "DAPI.Mean",
                                  "DAPI.StDev", "DAPI.CV", "Area", "Roundness",
                                  "Data")                               
head(H3K27_alldata_DAPI)

# get data H4K20me1
H4K20_rep1 <- read.delim("analyse_H4K20me1/H4K20me1 atub WT KO KD 070422__2022-04-07T15_53_56-Measurement 1/Evaluation4/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep1)

H4K20_rep2 <- read.delim("analyse_H4K20me1/H4K20me1 WT KO KD 210422__2022-04-21T15_47_58-Measurement 1/Evaluation3/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep2)

H4K20_rep3 <- read.delim("analyse_H4K20me1/H4K20me1 WT KO KD 110522__2022-05-11T15_47_27-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep3)

# check colnames are the same
table(colnames(H4K20_rep1) == colnames(H4K20_rep2))
table(colnames(H4K20_rep1) == colnames(H4K20_rep3))

# add replicate
H4K20_rep1 <- H4K20_rep1 %>%
    mutate(Replicate = rep("rep1", nrow(H4K20_rep1))) 
H4K20_rep2 <- H4K20_rep2 %>%
    mutate(Replicate = rep("rep2", nrow(H4K20_rep2))) 
H4K20_rep3 <- H4K20_rep3 %>%
    mutate(Replicate = rep("rep3", nrow(H4K20_rep3))) 

# put all data together -> since colnames are the same, do rbind
H4K20_alldata <- rbind(H4K20_rep1, H4K20_rep2, H4K20_rep3)

# check
table(H4K20_alldata$Replicate, H4K20_alldata$Row) 
table(H4K20_alldata$Replicate, H4K20_alldata$Column) 

# select DAPI info
H4K20_alldata_DAPI <- H4K20_alldata %>%
    dplyr::select(Row, Column, Replicate,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.StdDev,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.CV....,
                  `All.Nuclei.Selected...Area.Nuclei.Area..µm..`,
                  `All.Nuclei.Selected...Area.Nuclei.Roundness`) %>%
    mutate(Data = rep("H4K20", nrow(H4K20_alldata)))

colnames(H4K20_alldata_DAPI) <- c("Row", "Column", "Replicate", "DAPI.Mean",
                                  "DAPI.StDev", "DAPI.CV", "Area", "Roundness",
                                  "Data")
head(H4K20_alldata_DAPI)

# 2. Process input data ---------------------------------------------------------------------------

# put together all data
DAPI_all_data <- H3K9_alldata_DAPI %>%
    rbind(H3K27_alldata_DAPI) %>%
    rbind(H4K20_alldata_DAPI)

# check
table(DAPI_all_data$Column, DAPI_all_data$Row)
table(DAPI_all_data$Data, DAPI_all_data$Replicate)

# add condition
DAPI_all_data_cond <- DAPI_all_data %>%
    mutate(Condition = case_when(
        Row == 2 ~ "WT",
        Row == 4 ~ "KO1"),
        Replicate2 = case_when(
            Data == "H3K9" & Replicate == "rep1" ~ 1,
            Data == "H3K9" & Replicate == "rep2" ~ 2,
            Data == "H3K9" & Replicate == "rep3" ~ 3,
            Data == "H3K9" & Replicate == "rep4" ~ 4,
            Data == "H3K27" & Replicate == "rep1" ~ 5,
            Data == "H3K27" & Replicate == "rep2" ~ 6,
            Data == "H3K27" & Replicate == "rep3" ~ 7,
            Data == "H3K27" & Replicate == "rep4" ~ 8,
            Data == "H4K20" & Replicate == "rep1" ~ 9,
            Data == "H4K20" & Replicate == "rep2" ~ 10,
            Data == "H4K20" & Replicate == "rep3" ~ 11))

# check
table(DAPI_all_data_cond$Replicate, DAPI_all_data_cond$Row) 
table(DAPI_all_data_cond$Replicate, DAPI_all_data_cond$Column)  
table(DAPI_all_data_cond$Condition, DAPI_all_data_cond$Row)  
table(DAPI_all_data_cond$Data, DAPI_all_data_cond$Column)
table(DAPI_all_data_cond$Data, DAPI_all_data_cond$Replicate2)
table(DAPI_all_data_cond$Condition, DAPI_all_data_cond$Replicate2)

# filter single nuclei and round
DAPI_alldata_singl <- DAPI_all_data_cond %>%
    filter(Area < 200 & Area > 60, Roundness > 0.9)
nrow(DAPI_all_data_cond) 
nrow(DAPI_alldata_singl) 
table(DAPI_alldata_singl$Condition) 

# remove mitotic cells with DAPI 
# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_DAPI_data <- function(n) {
    repdata <- DAPI_alldata_singl %>%
        filter(Replicate2 == n) 
    repdata_scale <- repdata %>%
        mutate(logDAPI = log(DAPI.Mean),
               logDAPI_sc = range01(logDAPI))
}

reps <- 1:11
DAPI_data_scaled <- map(reps, scale_DAPI_data)

# build data frame
DAPI_alldata_scaled <- bind_rows(DAPI_data_scaled)
table(DAPI_alldata_scaled$Condition)

# mitotic cells have higher DAPI signal
hist(DAPI_alldata_scaled$DAPI.Mean)
hist(DAPI_alldata_scaled$logDAPI)
hist(DAPI_alldata_scaled$logDAPI_sc)
# we could try 0.7 as threshold

# remove mitotic cells
DAPI_alldata_scaled_nomit <- DAPI_alldata_scaled %>%
    filter(logDAPI_sc < 0.7)
nrow(DAPI_alldata_scaled_nomit) 
nrow(DAPI_alldata_scaled) 

# get WTvsKO dats
DAPI_KO_scaled_nomit <- DAPI_alldata_scaled_nomit %>%
    filter(Condition %in% c("WT", "KO1"))
table(DAPI_KO_scaled_nomit$Condition)

# convert condition to factor
DAPI_KO_scaled_nomit$Condition <- factor(DAPI_KO_scaled_nomit$Condition, levels = c("WT", "KO1"))

# create logCV
DAPI_KO_scaled_nomit <- DAPI_KO_scaled_nomit %>%
    mutate(log2CV = log2(DAPI.CV),
           log2DAPI = log2(DAPI.Mean))

# 3. Analyse and plot data ---------------------------------------------------------------------------

# statistics for each sample
DAPI_KO_stats <- DAPI_KO_scaled_nomit %>%
    group_by(Condition, Replicate2) %>%
    summarise(
        count = n(),
        mean_DAPI_mean = mean(DAPI.Mean),
        median_DAPI_mean = median(DAPI.Mean),
        sd_DAPI_mean = sd(DAPI.Mean),
        se_DAPI_mean = sd_DAPI_mean/sqrt(count),
        mean_DAPI_CV = mean(DAPI.CV),
        median_DAPI_CV = median(DAPI.CV),
        sd_DAPI_CV = sd(DAPI.CV),
        se_DAPI_CV = sd_DAPI_CV/sqrt(count))

# add DAPI facet
DAPI_KO_scaled_nomit$Mark <- rep("DAPI", nrow(DAPI_KO_scaled_nomit)) 
DAPI_KO_scaled_nomit$Mark <- factor(DAPI_KO_scaled_nomit$Mark)

# plot only KO2
my_comp<- list(c("WT", "KO1"))
ggplot(DAPI_KO_scaled_nomit, aes(x = Condition, y = log2DAPI)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58")) +
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
ggsave("plots/DAPI_allreps_KO_violin.pdf", width = 2, height = 2.2)

ggplot(DAPI_KO_scaled_nomit, aes(x = Condition, y = log2CV)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58")) +
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
ggsave("plots/DAPICV_allreps_KO_violin.pdf", width = 2, height = 2.2)