# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("ggrepel")
library("nortest")

# get data H3K9
H3K9_rep1_thym <- read.delim("prev repeat with new evaluations/H3K9me3 thym WT KO CD 40X 240522__2022-05-24T11_37_36-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep1_thym)

H3K9_rep1_RO <- read.delim("prev repeat with new evaluations/H3K9me3 RO3306 WT KO CD 40X 240522__2022-05-24T12_18_53-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep1_RO)

H3K9_rep2_thym <- read.delim("rep2/H3K9me3 wt-ko-nls 40X thym rep2__2023-09-13T17_10_32-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep2_thym)

H3K9_rep2_RO <- read.delim("rep2/H3K9me3 wt-ko-nls 40X RO rep2__2023-09-21T17_19_59-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9_rep2_RO)

# check colnames are the same
table(colnames(H3K9_rep1_thym) == colnames(H3K9_rep1_RO))
table(colnames(H3K9_rep1_thym) == colnames(H3K9_rep2_thym))
table(colnames(H3K9_rep1_thym) == colnames(H3K9_rep2_RO))

# check
table(H3K9_rep1_thym$Row, H3K9_rep1_thym$Column)
table(H3K9_rep1_RO$Row, H3K9_rep1_RO$Column)

table(H3K9_rep2_thym$Row, H3K9_rep2_thym$Column)
table(H3K9_rep2_RO$Row, H3K9_rep2_RO$Column)

# add replicate and treatment
H3K9_rep1_thym <- H3K9_rep1_thym %>%
    mutate(Replicate = rep(1, nrow(H3K9_rep1_thym)),
           Treatment = ifelse(Column == 3, "DMSO", "thymidine"))
H3K9_rep1_RO <- H3K9_rep1_RO %>%
    mutate(Replicate = rep(1, nrow(H3K9_rep1_RO)),
           Treatment = ifelse(Column == 3, "DMSO", "RO3306"))

H3K9_rep2_thym <- H3K9_rep2_thym %>%
    mutate(Replicate = rep(2, nrow(H3K9_rep2_thym)),
           Treatment = ifelse(Column == 2, "DMSO", "thymidine"))
H3K9_rep2_RO <- H3K9_rep2_RO %>%
    mutate(Replicate = rep(2, nrow(H3K9_rep2_RO)),
           Treatment = ifelse(Column == 2, "DMSO", "RO3306"))

# put all data together -> since colnames are the same, do rbind
H3K9_alldata <- rbind(H3K9_rep1_thym, H3K9_rep1_RO, H3K9_rep2_thym, H3K9_rep2_RO)

# check
table(H3K9_alldata$Replicate, H3K9_alldata$Row)
table(H3K9_alldata$Replicate, H3K9_alldata$Column)  
table(H3K9_alldata$Treatment, H3K9_alldata$Column) 

# check colnames
colnames(H3K9_alldata)

# select DAPI info
H3K9_alldata_DAPI <- H3K9_alldata %>%
    dplyr::select(Row, Column, Treatment, Replicate,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.StdDev,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.CV....,
                  `All.Nuclei.Selected...Area.Nuclei.Area..µm..`,
                  `All.Nuclei.Selected...Area.Nuclei.Roundness`) %>%
    mutate(Data = rep("H3K9", nrow(H3K9_alldata)))

colnames(H3K9_alldata_DAPI) <- c("Row", "Column", "Treatment", "Replicate",
                                 "DAPI.Mean", "DAPI.StDev", "DAPI.CV", 
                                 "Area", "Roundness", "Data")
head(H3K9_alldata_DAPI)

# get data H3K27me3
H3K27_rep1_thym <- read.delim("prev repeat with new evaluations/H3K27me3 thym WT KO CD 40X 240522__2022-05-24T12_55_54-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep1_thym)

H3K27_rep1_RO <- read.delim("prev repeat with new evaluations/H3K27me3 RO3306 WT KO CD 40X__2022-05-24T12_29_09-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep1_RO)

H3K27_rep2_thym <- read.delim("rep2/H3K27me3 wt-ko-nls 40X thym 2nd rep__2023-09-13T17_19_11-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep2_thym)

H3K27_rep2_RO <- read.delim("rep2/H3K27me3 wt-ko-nls 40X RO 2ndrep__2023-09-21T17_10_12-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27_rep2_RO)

# check colnames are the same
table(colnames(H3K27_rep1_thym) == colnames(H3K27_rep1_RO))
table(colnames(H3K27_rep1_thym) == colnames(H3K27_rep2_thym))
table(colnames(H3K27_rep1_thym) == colnames(H3K27_rep2_RO))

# check
table(H3K27_rep1_thym$Row, H3K27_rep1_thym$Column)
table(H3K27_rep1_RO$Row, H3K27_rep1_RO$Column)

table(H3K27_rep2_thym$Row, H3K27_rep2_thym$Column)
table(H3K27_rep2_RO$Row, H3K27_rep2_RO$Column)

# add replicate and treatment
H3K27_rep1_thym <- H3K27_rep1_thym %>%
    mutate(Replicate = rep(1, nrow(H3K27_rep1_thym)),
           Treatment = ifelse(Column == 4, "DMSO", "thymidine"))
H3K27_rep1_RO <- H3K27_rep1_RO %>%
    mutate(Replicate = rep(1, nrow(H3K27_rep1_RO)),
           Treatment = ifelse(Column == 4, "DMSO", "RO3306"))

H3K27_rep2_thym <- H3K27_rep2_thym %>%
    mutate(Replicate = rep(2, nrow(H3K27_rep2_thym)),
           Treatment = ifelse(Column == 3, "DMSO", "thymidine"))
H3K27_rep2_RO <- H3K27_rep2_RO %>%
    mutate(Replicate = rep(2, nrow(H3K27_rep2_RO)),
           Treatment = ifelse(Column == 3, "DMSO", "RO3306"))

# put all data together -> since colnames are the same, do rbind
H3K27_alldata <- rbind(H3K27_rep1_thym, H3K27_rep1_RO, H3K27_rep2_thym, H3K27_rep2_RO)

# check
table(H3K27_alldata$Replicate, H3K27_alldata$Row)
table(H3K27_alldata$Replicate, H3K27_alldata$Column)  
table(H3K27_alldata$Treatment, H3K27_alldata$Column)  

# select DAPI info
H3K27_alldata_DAPI <- H3K27_alldata %>%
    dplyr::select(Row, Column, Treatment, Replicate,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.StdDev,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.CV....,
                  `All.Nuclei.Selected...Area.Nuclei.Area..µm..`,
                  `All.Nuclei.Selected...Area.Nuclei.Roundness`) %>%
    mutate(Data = rep("H3K27", nrow(H3K27_alldata)))

colnames(H3K27_alldata_DAPI) <- c("Row", "Column", "Treatment", "Replicate",
                                  "DAPI.Mean","DAPI.StDev", "DAPI.CV", 
                                  "Area", "Roundness","Data")                               
head(H3K27_alldata_DAPI)

# get data H4K20me1
H4K20_rep1_thym <- read.delim("prev repeat with new evaluations/H4K20me1 thym WT KO CD 40X 240522__2022-05-24T12_03_40-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep1_thym)

H4K20_rep1_RO <- read.delim("prev repeat with new evaluations/H4K20me1 RO3306 WT KO CD 40X 240522__2022-05-24T12_41_34-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep1_RO)

H4K20_rep2_thym <- read.delim("rep2/H4K20me1 WTvsKOvsNLS 40X thymidine 2nd rep__2023-09-13T16_54_10-Measurement 2/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep2_thym)

H4K20_rep2_RO <- read.delim("rep2/H4K20me1 WTvsKOvsNLS 40X Ro 2ndrep__2023-09-21T16_56_43-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20_rep2_RO)

# check colnames are the same
table(colnames(H4K20_rep1_thym) == colnames(H4K20_rep1_RO))
table(colnames(H4K20_rep1_thym) == colnames(H4K20_rep2_thym))
table(colnames(H4K20_rep1_thym) == colnames(H4K20_rep2_RO))

# check
table(H4K20_rep1_thym$Row, H4K20_rep1_thym$Column)
table(H4K20_rep1_RO$Row, H4K20_rep1_RO$Column)

table(H4K20_rep2_thym$Row, H4K20_rep2_thym$Column)
table(H4K20_rep2_RO$Row, H4K20_rep2_RO$Column)

# add replicate and treatment
H4K20_rep1_thym <- H4K20_rep1_thym %>%
    mutate(Replicate = rep(1, nrow(H4K20_rep1_thym)),
           Treatment = ifelse(Column == 5, "DMSO", "thymidine"))
H4K20_rep1_RO <- H4K20_rep1_RO %>%
    mutate(Replicate = rep(1, nrow(H4K20_rep1_RO)),
           Treatment = ifelse(Column == 5, "DMSO", "RO3306"))

H4K20_rep2_thym <- H4K20_rep2_thym %>%
    mutate(Replicate = rep(2, nrow(H4K20_rep2_thym)),
           Treatment = ifelse(Column == 4, "DMSO", "thymidine"))
H4K20_rep2_RO <- H4K20_rep2_RO %>%
    mutate(Replicate = rep(2, nrow(H4K20_rep2_RO)),
           Treatment = ifelse(Column == 4, "DMSO", "RO3306"))

# put all data together -> since colnames are the same, do rbind
H4K20_alldata <- rbind(H4K20_rep1_thym, H4K20_rep1_RO, H4K20_rep2_thym, H4K20_rep2_RO)

# check
table(H4K20_alldata$Replicate, H4K20_alldata$Row)
table(H4K20_alldata$Replicate, H4K20_alldata$Column) 
table(H4K20_alldata$Treatment, H4K20_alldata$Column) 

# select DAPI info
H4K20_alldata_DAPI <- H4K20_alldata %>%
    dplyr::select(Row, Column, Treatment, Replicate,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.StdDev,
                  All.Nuclei.Selected...Intensity.Nucleus.DAPI.CV....,
                  `All.Nuclei.Selected...Area.Nuclei.Area..µm..`,
                  `All.Nuclei.Selected...Area.Nuclei.Roundness`) %>%
    mutate(Data = rep("H4K20", nrow(H4K20_alldata)))

colnames(H4K20_alldata_DAPI) <- c("Row", "Column", "Treatment", "Replicate",
                                  "DAPI.Mean","DAPI.StDev", "DAPI.CV",
                                  "Area","Roundness","Data")
head(H4K20_alldata_DAPI)

# put together all data
DAPI_all_data <- H3K9_alldata_DAPI %>%
    rbind(H3K27_alldata_DAPI) %>%
    rbind(H4K20_alldata_DAPI)

# check
table(DAPI_all_data$Column, DAPI_all_data$Row, DAPI_all_data$Treatment)
table(DAPI_all_data$Data, DAPI_all_data$Replicate)
table(DAPI_all_data$Row, DAPI_all_data$Data, DAPI_all_data$Replicate)

# add condition
DAPI_all_data_cond <- DAPI_all_data %>%
    mutate(Condition = case_when(
        Row == 2 ~ "WT",
        Row == 4 ~ "KO1"),
        Replicate2 = case_when(
            Data == "H3K9" & Replicate == 1 ~ 1,
            Data == "H3K9" & Replicate == 2 ~ 2,
            Data == "H3K27" & Replicate == 1 ~ 3,
            Data == "H3K27" & Replicate == 2 ~ 4,
            Data == "H4K20" & Replicate == 1 ~ 5,
            Data == "H4K20" & Replicate == 2 ~ 6))

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
nrow(DAPI_all_data_cond) # 183172
nrow(DAPI_alldata_singl) # 104085
table(DAPI_alldata_singl$Condition) 

# remove mitotic cells with DAPI 
# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_DAPI_data <- function(n, treat) {
    repdata_treat <- DAPI_alldata_singl %>%
        filter(Replicate2 == n & Treatment == treat) 
    repdata_scale <- repdata_treat %>%
        mutate(logDAPI = log(DAPI.Mean),
               logDAPI_sc = range01(logDAPI))
}
reps <- list(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6)
treats <- list("DMSO", "thymidine", "RO3306", "DMSO", "thymidine", "RO3306",
               "DMSO", "thymidine", "RO3306", "DMSO", "thymidine", "RO3306",
               "DMSO", "thymidine", "RO3306", "DMSO", "thymidine", "RO3306")
DAPI_data_scaled <- map2(reps, treats, scale_DAPI_data)

# build data frame
DAPI_alldata_scaled <- bind_rows(DAPI_data_scaled)
table(DAPI_alldata_scaled$Condition, DAPI_alldata_scaled$Treatment)

# mitotic cells have higher DAPI signal
hist(DAPI_alldata_scaled$DAPI.Mean)
hist(DAPI_alldata_scaled$logDAPI)
hist(DAPI_alldata_scaled$logDAPI_sc)
# we could try 0.7 as threshold

# remove mitotic cells
DAPI_alldata_scaled_nomit <- DAPI_alldata_scaled %>%
    filter(logDAPI_sc < 0.7)
nrow(DAPI_alldata_scaled_nomit) # 98592
nrow(DAPI_alldata_scaled) # 104085

# observe mitotic DAPI values
DAPI_alldata_scaled_mit <- DAPI_alldata_scaled %>%
    filter(logDAPI_sc > 0.7)

# filter WT and KO1
DAPI_WT_KO_scaled_nomit <- DAPI_alldata_scaled_nomit %>%
    filter(Condition %in% c("WT", "KO1"))
table(DAPI_WT_KO_scaled_nomit$Condition, DAPI_WT_KO_scaled_nomit$Treatment)

# convert condition to factor
DAPI_WT_KO_scaled_nomit$Condition <- factor(DAPI_WT_KO_scaled_nomit$Condition, levels = c("WT", "KO1"))

DAPI_WT_KO_scaled_nomit$Treatment <- factor(DAPI_WT_KO_scaled_nomit$Treatment, levels = c("DMSO", "thymidine", "RO3306"))

# create logCV
DAPI_WT_KO_scaled_nomit <- DAPI_WT_KO_scaled_nomit %>%
    mutate(log2CV = log2(DAPI.CV),
           log2DAPI = log2(DAPI.Mean))

# add DAPI facet
DAPI_WT_KO_scaled_nomit$Mark <- rep("DAPI", nrow(DAPI_WT_KO_scaled_nomit)) 
DAPI_WT_KO_scaled_nomit$Mark <- factor(DAPI_WT_KO_scaled_nomit$Mark)

# plot violin
my_comp<- list(c("WT", "KO1"))

ggplot(DAPI_WT_KO_scaled_nomit, aes(x = Condition, y = log2CV)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_grid(Mark ~ Treatment) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58")) +
    stat_compare_means(comparisons = my_comp, size = 3, tip.length = 0.01, braket.size = 0.01,
                       method = "wilcox", label.y = 5.2) +
    labs(x = "", y = "log2(Coefficient of variation)", fill = "") +
    ylim(3.5,5.5) +
    theme_bw() +
    
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("DAPI_plots/DAPICV_allreps_WTKO1_treats_violin_hor.pdf",  width = 4.5, height = 2.2)