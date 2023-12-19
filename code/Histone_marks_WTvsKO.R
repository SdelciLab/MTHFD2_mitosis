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

H3K9_rep5a <- read.delim("analyse_H3K9me3/H3K9me3 RO3306 WT KO CD 40X 240522__2022-05-24T12_18_53-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep5a)

H3K9_rep5b <- read.delim("analyse_H3K9me3/H3K9me3 thym WT KO CD 40X 240522__2022-05-24T11_37_36-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep5b)

# check colnames are the same
table(colnames(H3K9_rep1) == colnames(H3K9_rep2))
table(colnames(H3K9_rep1) == colnames(H3K9_rep3))
table(colnames(H3K9_rep1) == colnames(H3K9_rep4))
table(colnames(H3K9_rep1) == colnames(H3K9_rep5a))
table(colnames(H3K9_rep1) == colnames(H3K9_rep5b))

# add replicate
H3K9_rep1 <- H3K9_rep1 %>%
    mutate(Replicate = rep(1, nrow(H3K9_rep1)))
H3K9_rep2 <- H3K9_rep2 %>%
    mutate(Replicate = rep(2, nrow(H3K9_rep2)))
H3K9_rep3 <- H3K9_rep3 %>%
    mutate(Replicate = rep(3, nrow(H3K9_rep3)))
H3K9_rep4 <- H3K9_rep4 %>%
    mutate(Replicate = rep(4, nrow(H3K9_rep4)))
H3K9_rep5a <- H3K9_rep5a %>%
    mutate(Replicate = rep(5, nrow(H3K9_rep5a)))
H3K9_rep5b <- H3K9_rep5b %>%
    mutate(Replicate = rep(5, nrow(H3K9_rep5b)))

# put all data together -> since colnames are the same, do rbind
H3K9_alldata <- rbind(H3K9_rep1, H3K9_rep2, H3K9_rep3, H3K9_rep4, H3K9_rep5a, H3K9_rep5b)

# check
table(H3K9_alldata$Replicate, H3K9_alldata$Row) # rows 1-WT 2-KO1 3-KO2 4-shC 5-shMTHFD2 or rows 1-WT 2-KO1 3-KO2 4-CD
table(H3K9_alldata$Replicate, H3K9_alldata$Column)  

# add condition
# select only WT and KO1 
# remove weird points with StDev 0
H3K9_alldata_int <- H3K9_alldata %>%
    mutate(Condition = case_when(
        Row == 2 ~ "WT",
        Row == 4 ~ "KO1"),
        log2H3K9me3 = log2(All.Nuclei.Selected...Intensity.Nucleus.H3K9m3.Mean),
        IntDen = All.Nuclei.Selected...Intensity.Nucleus.H3K9m3.Mean * All.Nuclei.Selected...Area.Nuclei.Area..µm..,
        log2IntDen = log2(IntDen),
        log2CV = log2(All.Nuclei.Selected...Intensity.Nucleus.H3K9m3.CV....),
        logDAPI = log(All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)) %>%
    filter(Condition %in% c("WT", "KO1")) %>%
    filter(!(Replicate == 5 & (Column == 6 | Column == 7))) %>%
    filter(All.Nuclei.Selected...Intensity.Nucleus.H3K9m3.CV.... != 0)
colnames(H3K9_alldata_int)[18] <- "H3K9m3_CV"

# check
table(H3K9_alldata_int$Replicate, H3K9_alldata_int$Row) 
table(H3K9_alldata_int$Replicate, H3K9_alldata_int$Column)  

# filter single nuclei and round
H3K9_alldata_singl <- H3K9_alldata_int %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200 & `All.Nuclei.Selected...Area.Nuclei.Area..µm..` > 60,
           `All.Nuclei.Selected...Area.Nuclei.Roundness` > 0.9)
nrow(H3K9_alldata_singl) 
nrow(H3K9_alldata_int) 
table(H3K9_alldata_singl$Condition) 

# remove mitotic cells with DAPI 
# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_DAPI_data <- function(n) {
    repdata <- H3K9_alldata_singl %>%
        filter(Replicate == n) 
    repdata_scale <- repdata %>%
        mutate(logDAPI_sc = range01(logDAPI))
}

reps <- 1:5
H3K9_data_scaled_DAPI <- map(reps, scale_DAPI_data)

# build data frame
H3K9_alldata_scaled_DAPI <- bind_rows(H3K9_data_scaled_DAPI)
table(H3K9_alldata_scaled_DAPI$Condition)

# mitotic cells have higher DAPI signal
hist(H3K9_alldata_scaled_DAPI$All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)
hist(H3K9_alldata_scaled_DAPI$logDAPI)
hist(H3K9_alldata_scaled_DAPI$logDAPI_sc)
# we could try 0.7 as threshold

# remove mitotic cells
H3K9_alldata_scaled_DAPI_nomit <- H3K9_alldata_scaled_DAPI %>%
    filter(logDAPI_sc < 0.7)
nrow(H3K9_alldata_scaled_DAPI_nomit) 
nrow(H3K9_alldata_scaled_DAPI) 
# export data
save(H3K9_alldata_scaled_DAPI_nomit, file="analyse_H3K9me3/H3K9_alldata_nomit.RData")

# get data H3K27
H3K27_rep1 <- read.delim("analyse_H3K27me3/H3K27me3 atub WT KO KD 40X 030322__2022-03-03T11_30_02-Measurement 1/Evaluation3/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep1)

H3K27_rep2 <- read.delim("analyse_H3K27me3/H3K27me3 atub WT KO KD 40X 070422__2022-04-07T15_36_17-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep2)

H3K27_rep3 <- read.delim("analyse_H3K27me3/H3K27me3 WT KO KD 40X 210422__2022-04-21T15_22_27-Measurement 1/Evaluation3/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep3)

H3K27_rep4 <- read.delim("analyse_H3K27me3/H3K27me3 WT KO KD 40X 110522__2022-05-11T15_33_25-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep4)

H3K27_rep5a <- read.delim("analyse_H3K27me3/H3K27me3 RO3306 WT KO CD 40X__2022-05-24T12_29_09-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep5a)

H3K27_rep5b <- read.delim("analyse_H3K27me3/H3K27me3 thym WT KO CD 40X 240522__2022-05-24T12_55_54-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep5b)

# check colnames are the same
table(colnames(H3K27_rep1) == colnames(H3K27_rep2))
table(colnames(H3K27_rep1) == colnames(H3K27_rep3))
table(colnames(H3K27_rep1) == colnames(H3K27_rep4))
table(colnames(H3K27_rep1) == colnames(H3K27_rep5a))
table(colnames(H3K27_rep1) == colnames(H3K27_rep5b))

# add replicate
H3K27_rep1 <- H3K27_rep1 %>%
    mutate(Replicate = rep(1, nrow(H3K27_rep1)))
H3K27_rep2 <- H3K27_rep2 %>%
    mutate(Replicate = rep(2, nrow(H3K27_rep2)))
H3K27_rep3 <- H3K27_rep3 %>%
    mutate(Replicate = rep(3, nrow(H3K27_rep3)))
H3K27_rep4 <- H3K27_rep4 %>%
    mutate(Replicate = rep(4, nrow(H3K27_rep4)))
H3K27_rep5a <- H3K27_rep5a %>%
    mutate(Replicate = rep(5, nrow(H3K27_rep5a)))
H3K27_rep5b <- H3K27_rep5b %>%
    mutate(Replicate = rep(5, nrow(H3K27_rep5b)))

# put all data together -> since colnames are the same, do rbind
H3K27_alldata <- rbind(H3K27_rep1, H3K27_rep2, H3K27_rep3, H3K27_rep4, H3K27_rep5a, H3K27_rep5b)

# check
table(H3K27_alldata$Replicate, H3K27_alldata$Row) 
table(H3K27_alldata$Replicate, H3K27_alldata$Column)  

# add condition
# select only WT, KO1 
H3K27_alldata_int <- H3K27_alldata %>%
    mutate(Condition = case_when(
        Row == 2 ~ "WT",
        Row == 4 ~ "KO1"),
        log2H3K27me3 = log2(All.Nuclei.Selected...Intensity.Nucleus.H3K27m3.Mean),
        IntDen = All.Nuclei.Selected...Intensity.Nucleus.H3K27m3.Mean * All.Nuclei.Selected...Area.Nuclei.Area..µm..,
        log2IntDen = log2(IntDen),
        log2CV = log2(All.Nuclei.Selected...Intensity.Nucleus.H3K27m3.CV....),
        logDAPI = log(All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)) %>%
    filter(Condition %in% c("WT", "KO1")) %>%
    filter(!(Replicate == "rep5" & (Column == 8 | Column == 9))) %>%
    filter(All.Nuclei.Selected...Intensity.Nucleus.H3K27m3.CV.... != 0)
colnames(H3K27_alldata_int)[18] <- "H3K27m3_CV"

# check
table(H3K27_alldata_int$Replicate, H3K27_alldata_int$Row) 
table(H3K27_alldata_int$Replicate, H3K27_alldata_int$Column)  

# filter single nuclei and round
H3K27_alldata_singl <- H3K27_alldata_int %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200 & `All.Nuclei.Selected...Area.Nuclei.Area..µm..` > 60,
           `All.Nuclei.Selected...Area.Nuclei.Roundness` > 0.9)
nrow(H3K27_alldata_singl) 
nrow(H3K27_alldata_int) 
table(H3K27_alldata_singl$Condition) 

# remove mitotic cells with DAPI 
# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_DAPI_data <- function(n) {
    repdata <- H3K27_alldata_singl %>%
        filter(Replicate == n) 
    repdata_scale <- repdata %>%
        mutate(logDAPI_sc = range01(logDAPI))
}

reps <- 1:5
H3K27_data_scaled_DAPI <- map(reps, scale_DAPI_data)

# build data frame
H3K27_alldata_scaled_DAPI <- bind_rows(H3K27_data_scaled_DAPI)
table(H3K27_alldata_scaled_DAPI$Condition)

# mitotic cells have higher DAPI signal
hist(H3K27_alldata_scaled_DAPI$All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)
hist(H3K27_alldata_scaled_DAPI$logDAPI)
hist(H3K27_alldata_scaled_DAPI$logDAPI_sc)
# we could try 0.7 as threshold

# remove mitotic cells
H3K27_alldata_scaled_DAPI_nomit <- H3K27_alldata_scaled_DAPI %>%
    filter(logDAPI_sc < 0.7)
nrow(H3K27_alldata_scaled_DAPI_nomit) 
nrow(H3K27_alldata_scaled_DAPI)
# export data
save(H3K27_alldata_scaled_DAPI_nomit, file="analyse_H3K27me3/H3K27_alldata_nomit.RData")

# get data H4K20
H4K20_rep1 <- read.delim("analyse_H4K20me1/H4K20me1 atub WT KO KD 070422__2022-04-07T15_53_56-Measurement 1/Evaluation4/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep1)

H4K20_rep2 <- read.delim("analyse_H4K20me1/H4K20me1 WT KO KD 210422__2022-04-21T15_47_58-Measurement 1/Evaluation3/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep2)

H4K20_rep3 <- read.delim("analyse_H4K20me1/H4K20me1 WT KO KD 110522__2022-05-11T15_47_27-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep3)

H4K20_rep4a <- read.delim("analyse_H4K20me1/H4K20me1 RO3306 WT KO CD 40X 240522__2022-05-24T12_41_34-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep4a)

H4K20_rep4b <- read.delim("analyse_H4K20me1/H4K20me1 thym WT KO CD 40X 240522__2022-05-24T12_03_40-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep4b)

# check colnames are the same
table(colnames(H4K20_rep1) == colnames(H4K20_rep2))
table(colnames(H4K20_rep1) == colnames(H4K20_rep3))
table(colnames(H4K20_rep1) == colnames(H4K20_rep4a))
table(colnames(H4K20_rep1) == colnames(H4K20_rep4b))

# add replicate
H4K20_rep1 <- H4K20_rep1 %>%
    mutate(Replicate = rep(1, nrow(H4K20_rep1))) 
H4K20_rep2 <- H4K20_rep2 %>%
    mutate(Replicate = rep(2, nrow(H4K20_rep2))) 
H4K20_rep3 <- H4K20_rep3 %>%
    mutate(Replicate = rep(3, nrow(H4K20_rep3))) 
H4K20_rep4a <- H4K20_rep4a %>%
    mutate(Replicate = rep(4, nrow(H4K20_rep4a)))
H4K20_rep4b <- H4K20_rep4b %>%
    mutate(Replicate = rep(4, nrow(H4K20_rep4b))) 

# put all data together -> since colnames are the same, do rbind
H4K20_alldata <- rbind(H4K20_rep1, H4K20_rep2, H4K20_rep3, H4K20_rep4a, H4K20_rep4b)

# check
table(H4K20_alldata$Replicate, H4K20_alldata$Row) 
table(H4K20_alldata$Replicate, H4K20_alldata$Column)  

# add condition
# select only WT, KO1 
H4K20_alldata_int <- H4K20_alldata %>%
    mutate(Condition = case_when(
        Row == 2 ~ "WT",
        Row == 4 ~ "KO1"),
        log2H4K20me1 = log2(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.Mean),
        IntDen = All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.Mean * All.Nuclei.Selected...Area.Nuclei.Area..µm..,
        log2IntDen = log2(IntDen),
        log2CV = log2(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.CV....),
        logDAPI = log(All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)) %>%
    filter(Condition %in% c("WT", "KO1")) %>%
    filter(!(Replicate == 4 & (Column == 10 | Column == 11))) %>%
    filter(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.CV.... != 0)
colnames(H4K20_alldata_int)[18] <- "H4K20m1_CV"

# check
table(H4K20_alldata_int$Replicate, H4K20_alldata_int$Row) 
table(H4K20_alldata_int$Replicate, H4K20_alldata_int$Column)  

# filter single nuclei and round
H4K20_alldata_singl <- H4K20_alldata_int %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200 & `All.Nuclei.Selected...Area.Nuclei.Area..µm..` > 60,
           `All.Nuclei.Selected...Area.Nuclei.Roundness` > 0.9)
nrow(H4K20_alldata_singl) 
nrow(H4K20_alldata_int) 
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

reps <- 1:4
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
nrow(H4K20_alldata_scaled_DAPI_nomit) # 
nrow(H4K20_alldata_scaled_DAPI) 
# export data
save(H4K20_alldata_scaled_DAPI_nomit, file="analyse_H4K20me1/H4K20_alldata_nomit.RData")

# 2. Analyse and plot data ---------------------------------------------------------------------------

# load data
load("analyse_H3K9me3/H3K9_alldata_nomit.RData")
load("analyse_H3K27me3/H3K27_alldata_nomit.RData")
load("analyse_H4K20me1/H4K20_alldata_nomit.RData")

# select useful data
H4K20_alldata_int <- H4K20_alldata_scaled_DAPI_nomit %>%
    dplyr::select(Condition, Replicate, log2H4K20me1, log2IntDen, H4K20m1_CV, log2CV) %>%
    mutate(Mark = rep("H4K20m1", nrow(H4K20_alldata_scaled_DAPI_nomit))) 
nrow(H4K20_alldata_int) 
H4K20_alldata_int$Condition <- factor(H4K20_alldata_int$Condition, levels = c("WT", "KO1"))
H4K20_alldata_int$Replicate <- factor(H4K20_alldata_int$Replicate, levels = c(1,2,3,4))

H3K9_alldata_int <- H3K9_alldata_scaled_DAPI_nomit %>%
    dplyr::select(Condition, Replicate, log2H3K9me3, log2IntDen, H3K9m3_CV, log2CV) %>%
    mutate(Mark = rep("H3K9m3", nrow(H3K9_alldata_scaled_DAPI_nomit))) 
nrow(H3K9_alldata_int) 
H3K9_alldata_int$Condition <- factor(H3K9_alldata_int$Condition, levels = c("WT", "KO1"))
H3K9_alldata_int$Replicate <- factor(H3K9_alldata_int$Replicate, levels = c(1,2,3,4,5))

H3K27_alldata_int <- H3K27_alldata_scaled_DAPI_nomit %>%
    dplyr::select(Condition, Replicate, log2H3K27me3, log2IntDen, H3K27m3_CV, log2CV) %>%
    mutate(Mark = rep("H3K27m3", nrow(H3K27_alldata_scaled_DAPI_nomit))) 
nrow(H3K27_alldata_int) 
H3K27_alldata_int$Condition <- factor(H3K27_alldata_int$Condition, levels = c("WT", "KO1"))
H3K27_alldata_int$Replicate <- factor(H3K27_alldata_int$Replicate, levels = c(1,2,3,4,5))

# put together all data
colnames(H4K20_alldata_int) <- c("Condition", "Replicate", "log2MeanInt", "log2IntDen", "CV", "log2CV", "Mark")
colnames(H3K9_alldata_int) <- c("Condition", "Replicate", "log2MeanInt", "log2IntDen", "CV", "log2CV", "Mark")
colnames(H3K27_alldata_int) <- c("Condition", "Replicate", "log2MeanInt", "log2IntDen", "CV", "log2CV", "Mark")

# put everything together
all_marks_data <- H4K20_alldata_int %>%
    rbind(H3K9_alldata_int) %>%
    rbind(H3K27_alldata_int)
all_marks_data$Mark <- factor(all_marks_data$Mark, levels = c("H4K20m1", "H3K9m3", "H3K27m3"))


# plot without correction KO1 and KO2
my_comp_KO <- list(c("WT", "KO1"))
all_marks_data$Condition <- factor(all_marks_data$Condition, levels = c("WT", "KO1"))

ggplot(all_marks_data, aes(x = Condition, y = log2MeanInt)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58")) +
    stat_compare_means(comparisons = my_comp_KO, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox") +
    labs(x = "", y = "log2(Mean Intensity)", fill = "") +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none")
ggsave("global_plots/log2MeanInt_allreps_violin_KO.pdf", width = 5.5, height = 2.7)

ggplot(all_marks_data, aes(x = Condition, y = log2CV)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_wrap(~Mark) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58")) +
    stat_compare_means(comparisons = my_comp_KO, size = 3, tip.length = 0.01, braket.size = 0.01,
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
ggsave("global_plots/log2CVallreps_violin_KO.pdf", device = "pdf", width = 5.5, height = 2.7)

