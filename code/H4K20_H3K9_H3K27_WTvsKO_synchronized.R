# libraries
library("tidyverse")
library("ggpubr")
library("car")
library("ggsignif")
library("ggrepel")
library("nortest")

### H3K9me3 data

# get data
H3K9m3_rep1_thym <- read.delim("prev repeat with new evaluations/H3K9me3 thym WT KO CD 40X 240522__2022-05-24T11_37_36-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9m3_rep1_thym)

H3K9m3_rep1_RO <- read.delim("prev repeat with new evaluations/H3K9me3 RO3306 WT KO CD 40X 240522__2022-05-24T12_18_53-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9m3_rep1_RO)

H3K9m3_rep2_thym <- read.delim("rep2/H3K9me3 wt-ko-nls 40X thym rep2__2023-09-13T17_10_32-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9m3_rep2_thym)

H3K9m3_rep2_RO <- read.delim("rep2/H3K9me3 wt-ko-nls 40X RO rep2__2023-09-21T17_19_59-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K9m3_rep2_RO)

# check colnames are the same
table(colnames(H3K9m3_rep1_thym) == colnames(H3K9m3_rep1_RO))
table(colnames(H3K9m3_rep1_thym) == colnames(H3K9m3_rep2_thym))
table(colnames(H3K9m3_rep1_thym) == colnames(H3K9m3_rep2_RO))

# check
table(H3K9m3_rep1_thym$Row, H3K9m3_rep1_thym$Column)
table(H3K9m3_rep1_RO$Row, H3K9m3_rep1_RO$Column)

table(H3K9m3_rep2_thym$Row, H3K9m3_rep2_thym$Column)
table(H3K9m3_rep2_RO$Row, H3K9m3_rep2_RO$Column)

# add replicate and treatment
H3K9m3_rep1_thym <- H3K9m3_rep1_thym %>%
    mutate(Replicate = rep(1, nrow(H3K9m3_rep1_thym)),
           Treatment = ifelse(Column == 3, "DMSO", "thymidine"))
H3K9m3_rep1_RO <- H3K9m3_rep1_RO %>%
    mutate(Replicate = rep(1, nrow(H3K9m3_rep1_RO)),
           Treatment = ifelse(Column == 3, "DMSO", "RO3306"))

H3K9m3_rep2_thym <- H3K9m3_rep2_thym %>%
    mutate(Replicate = rep(2, nrow(H3K9m3_rep2_thym)),
           Treatment = ifelse(Column == 2, "DMSO", "thymidine"))
H3K9m3_rep2_RO <- H3K9m3_rep2_RO %>%
    mutate(Replicate = rep(2, nrow(H3K9m3_rep2_RO)),
           Treatment = ifelse(Column == 2, "DMSO", "RO3306"))

# put all data together -> since colnames are the same, do rbind
H3K9_alldata <- rbind(H3K9m3_rep1_thym, H3K9m3_rep1_RO, 
                      H3K9m3_rep2_thym, H3K9m3_rep2_RO)

# check
table(H3K9_alldata$Treatment, H3K9_alldata$Column) 
table(H3K9_alldata$Replicate, H3K9_alldata$Row) 
table(H3K9_alldata$Replicate, H3K9_alldata$Column) 

# add condition
# select only WT and KO
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
    filter(All.Nuclei.Selected...Intensity.Nucleus.H3K9m3.CV.... != 0)
colnames(H3K9_alldata_int)[18] <- "H3K9m3_CV"

# check
table(H3K9_alldata_int$Replicate, H3K9_alldata_int$Row) 
table(H3K9_alldata_int$Replicate, H3K9_alldata_int$Column)  

table(H3K9_alldata_int$Condition, H3K9_alldata_int$Row) 
table(H3K9_alldata_int$Treatment, H3K9_alldata_int$Column)  

# filter single nuclei and round
H3K9_alldata_singl <- H3K9_alldata_int %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200 & `All.Nuclei.Selected...Area.Nuclei.Area..µm..` > 60,
           `All.Nuclei.Selected...Area.Nuclei.Roundness` > 0.9)
nrow(H3K9_alldata_singl) # 20116
nrow(H3K9_alldata_int) # 36592
table(H3K9_alldata_singl$Condition) 

# remove mitotic cells with DAPI 
# only for DMSO

# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_DAPI_data <- function(n, treat) {
    repdata_treat <- H3K9_alldata_singl %>%
        filter(Replicate == n & Treatment == treat) 
    repdata_scale <- repdata_treat %>%
        mutate(logDAPI_sc = range01(logDAPI))
}
reps <- list(1,1,1,2,2,2)
treats <- list("DMSO", "thymidine", "RO3306", "DMSO", "thymidine", "RO3306")
H3K9_data_scaled_DAPI <- map2(reps, treats, scale_DAPI_data)

# check
table(H3K9_data_scaled_DAPI[[1]]["Replicate"])
table(H3K9_data_scaled_DAPI[[2]]["Replicate"])
table(H3K9_data_scaled_DAPI[[3]]["Replicate"])
table(H3K9_data_scaled_DAPI[[4]]["Replicate"])
table(H3K9_data_scaled_DAPI[[5]]["Replicate"])
table(H3K9_data_scaled_DAPI[[6]]["Replicate"])

table(H3K9_data_scaled_DAPI[[1]]["Treatment"])
table(H3K9_data_scaled_DAPI[[2]]["Treatment"])
table(H3K9_data_scaled_DAPI[[3]]["Treatment"])
table(H3K9_data_scaled_DAPI[[4]]["Treatment"])
table(H3K9_data_scaled_DAPI[[5]]["Treatment"])
table(H3K9_data_scaled_DAPI[[6]]["Treatment"])

# build data frame
H3K9_alldata_scaled_DAPI <- bind_rows(H3K9_data_scaled_DAPI)
table(H3K9_alldata_scaled_DAPI$Condition)
table(H3K9_alldata_scaled_DAPI$Treatment)

# mitotic cells have higher DAPI signal
hist(H3K9_alldata_scaled_DAPI$All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)
hist(H3K9_alldata_scaled_DAPI$logDAPI)
hist(H3K9_alldata_scaled_DAPI$logDAPI_sc)
# we could try 0.7 as threshold

# remove mitotic cells only in DMSO
H3K9_alldata_scaled_DAPI_nomit <- H3K9_alldata_scaled_DAPI %>%
    filter(Treatment != "DMSO" | (Treatment == "DMSO" & logDAPI_sc < 0.7))
nrow(H3K9_alldata_scaled_DAPI_nomit) # 19281
nrow(H3K9_alldata_scaled_DAPI) # 20116

# check
table(H3K9_alldata_scaled_DAPI$Treatment)
table(H3K9_alldata_scaled_DAPI_nomit$Treatment)

# export data
save(H3K9_alldata_scaled_DAPI_nomit, file="H3K9_alldata_nomit.RData")

### H3K27me3 data

# get data
H3K27m3_rep1_thym <- read.delim("prev repeat with new evaluations/H3K27me3 thym WT KO CD 40X 240522__2022-05-24T12_55_54-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27m3_rep1_thym)

H3K27m3_rep1_RO <- read.delim("prev repeat with new evaluations/H3K27me3 RO3306 WT KO CD 40X__2022-05-24T12_29_09-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27m3_rep1_RO)

H3K27m3_rep2_thym <- read.delim("rep2/H3K27me3 wt-ko-nls 40X thym 2nd rep__2023-09-13T17_19_11-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27m3_rep2_thym)

H3K27m3_rep2_RO <- read.delim("rep2/H3K27me3 wt-ko-nls 40X RO 2ndrep__2023-09-21T17_10_12-Measurement 1/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H3K27m3_rep2_RO)

# check colnames are the same
table(colnames(H3K27m3_rep1_thym) == colnames(H3K27m3_rep1_RO))
table(colnames(H3K27m3_rep1_thym) == colnames(H3K27m3_rep2_thym))
table(colnames(H3K27m3_rep1_thym) == colnames(H3K27m3_rep2_RO))

# check
table(H3K27m3_rep1_thym$Row, H3K27m3_rep1_thym$Column)
table(H3K27m3_rep1_RO$Row, H3K27m3_rep1_RO$Column)

table(H3K27m3_rep2_thym$Row, H3K27m3_rep2_thym$Column)
table(H3K27m3_rep2_RO$Row, H3K27m3_rep2_RO$Column)

# add replicate and treatment
H3K27m3_rep1_thym <- H3K27m3_rep1_thym %>%
    mutate(Replicate = rep(1, nrow(H3K27m3_rep1_thym)),
           Treatment = ifelse(Column == 4, "DMSO", "thymidine"))
H3K27m3_rep1_RO <- H3K27m3_rep1_RO %>%
    mutate(Replicate = rep(1, nrow(H3K27m3_rep1_RO)),
           Treatment = ifelse(Column == 4, "DMSO", "RO3306"))

H3K27m3_rep2_thym <- H3K27m3_rep2_thym %>%
    mutate(Replicate = rep(2, nrow(H3K27m3_rep2_thym)),
           Treatment = ifelse(Column == 3, "DMSO", "thymidine"))
H3K27m3_rep2_RO <- H3K27m3_rep2_RO %>%
    mutate(Replicate = rep(2, nrow(H3K27m3_rep2_RO)),
           Treatment = ifelse(Column == 3, "DMSO", "RO3306"))

# put all data together -> since colnames are the same, do rbind
H3K27_alldata <- rbind(H3K27m3_rep1_thym, H3K27m3_rep1_RO, 
                       H3K27m3_rep2_thym, H3K27m3_rep2_RO)

# check
table(H3K27_alldata$Treatment, H3K27_alldata$Column) 
table(H3K27_alldata$Replicate, H3K27_alldata$Row) 
table(H3K27_alldata$Replicate, H3K27_alldata$Column)  

# add condition
# select only WT and KO
# remove weird points with StDev 0
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
    filter(All.Nuclei.Selected...Intensity.Nucleus.H3K27m3.CV.... != 0)
colnames(H3K27_alldata_int)[18] <- "H3K27m3_CV"

# check
table(H3K27_alldata_int$Replicate, H3K27_alldata_int$Row) 
table(H3K27_alldata_int$Replicate, H3K27_alldata_int$Column)  

table(H3K27_alldata_int$Condition, H3K27_alldata_int$Row) 
table(H3K27_alldata_int$Treatment, H3K27_alldata_int$Column)  

# filter single nuclei and round
H3K27_alldata_singl <- H3K27_alldata_int %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200 & `All.Nuclei.Selected...Area.Nuclei.Area..µm..` > 60,
           `All.Nuclei.Selected...Area.Nuclei.Roundness` > 0.9)
nrow(H3K27_alldata_singl) # 21942
nrow(H3K27_alldata_int) # 37209
table(H3K27_alldata_singl$Condition) 

# remove mitotic cells with DAPI 
# only for DMSO

# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_DAPI_data <- function(n, treat) {
    repdata_treat <- H3K27_alldata_singl %>%
        filter(Replicate == n & Treatment == treat) 
    repdata_scale <- repdata_treat %>%
        mutate(logDAPI_sc = range01(logDAPI))
}
reps <- list(1,1,1,2,2,2)
treats <- list("DMSO", "thymidine", "RO3306", "DMSO", "thymidine", "RO3306")
H3K27_data_scaled_DAPI <- map2(reps, treats, scale_DAPI_data)

# check
table(H3K27_data_scaled_DAPI[[1]]["Replicate"])
table(H3K27_data_scaled_DAPI[[2]]["Replicate"])
table(H3K27_data_scaled_DAPI[[3]]["Replicate"])
table(H3K27_data_scaled_DAPI[[4]]["Replicate"])
table(H3K27_data_scaled_DAPI[[5]]["Replicate"])
table(H3K27_data_scaled_DAPI[[6]]["Replicate"])

table(H3K27_data_scaled_DAPI[[1]]["Treatment"])
table(H3K27_data_scaled_DAPI[[2]]["Treatment"])
table(H3K27_data_scaled_DAPI[[3]]["Treatment"])
table(H3K27_data_scaled_DAPI[[4]]["Treatment"])
table(H3K27_data_scaled_DAPI[[5]]["Treatment"])
table(H3K27_data_scaled_DAPI[[6]]["Treatment"])

# build data frame
H3K27_alldata_scaled_DAPI <- bind_rows(H3K27_data_scaled_DAPI)
table(H3K27_alldata_scaled_DAPI$Condition)
table(H3K27_alldata_scaled_DAPI$Treatment)

# mitotic cells have higher DAPI signal
hist(H3K27_alldata_scaled_DAPI$All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)
hist(H3K27_alldata_scaled_DAPI$logDAPI)
hist(H3K27_alldata_scaled_DAPI$logDAPI_sc)
# we could try 0.7 as threshold

# remove mitotic cells only in DMSO
H3K27_alldata_scaled_DAPI_nomit <- H3K27_alldata_scaled_DAPI %>%
    filter(Treatment != "DMSO" | (Treatment == "DMSO" & logDAPI_sc < 0.7))
nrow(H3K27_alldata_scaled_DAPI_nomit) # 21742
nrow(H3K27_alldata_scaled_DAPI) # 21942

# check
table(H3K27_alldata_scaled_DAPI$Treatment)
table(H3K27_alldata_scaled_DAPI_nomit$Treatment)

# export data
save(H3K27_alldata_scaled_DAPI_nomit, file="H3K27_alldata_nomit.RData")

### H4K20me1 data

# get data
H4K20m1_rep1_thym <- read.delim("prev repeat with new evaluations/H4K20me1 thym WT KO CD 40X 240522__2022-05-24T12_03_40-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20m1_rep1_thym)

H4K20m1_rep1_RO <- read.delim("prev repeat with new evaluations/H4K20me1 RO3306 WT KO CD 40X 240522__2022-05-24T12_41_34-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20m1_rep1_RO)

H4K20m1_rep2_thym <- read.delim("rep2/H4K20me1 WTvsKOvsNLS 40X thymidine 2nd rep__2023-09-13T16_54_10-Measurement 2/Evaluation1/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20m1_rep2_thym)

H4K20m1_rep2_RO <- read.delim("rep2/H4K20me1 WTvsKOvsNLS 40X Ro 2ndrep__2023-09-21T16_56_43-Measurement 1/Evaluation2/Objects_Population - All Nuclei Selected.txt", skip = 9)
head(H4K20m1_rep2_RO)

# check colnames are the same
table(colnames(H4K20m1_rep1_thym) == colnames(H4K20m1_rep1_RO))
table(colnames(H4K20m1_rep1_thym) == colnames(H4K20m1_rep2_thym))
table(colnames(H4K20m1_rep1_thym) == colnames(H4K20m1_rep2_RO))

# check
table(H4K20m1_rep1_thym$Row, H4K20m1_rep1_thym$Column)
table(H4K20m1_rep1_RO$Row, H4K20m1_rep1_RO$Column)

table(H4K20m1_rep2_thym$Row, H4K20m1_rep2_thym$Column)
table(H4K20m1_rep2_RO$Row, H4K20m1_rep2_RO$Column)

# add replicate and treatment
H4K20m1_rep1_thym <- H4K20m1_rep1_thym %>%
    mutate(Replicate = rep(1, nrow(H4K20m1_rep1_thym)),
           Treatment = ifelse(Column == 5, "DMSO", "thymidine"))
H4K20m1_rep1_RO <- H4K20m1_rep1_RO %>%
    mutate(Replicate = rep(1, nrow(H4K20m1_rep1_RO)),
           Treatment = ifelse(Column == 5, "DMSO", "RO3306"))

H4K20m1_rep2_thym <- H4K20m1_rep2_thym %>%
    mutate(Replicate = rep(2, nrow(H4K20m1_rep2_thym)),
           Treatment = ifelse(Column == 4, "DMSO", "thymidine"))
H4K20m1_rep2_RO <- H4K20m1_rep2_RO %>%
    mutate(Replicate = rep(2, nrow(H4K20m1_rep2_RO)),
           Treatment = ifelse(Column == 4, "DMSO", "RO3306"))

# put all data together -> since colnames are the same, do rbind
H4K20_alldata <- rbind(H4K20m1_rep1_thym, H4K20m1_rep1_RO, 
                       H4K20m1_rep2_thym, H4K20m1_rep2_RO)

# check
table(H4K20_alldata$Treatment, H4K20_alldata$Column) 
table(H4K20_alldata$Replicate, H4K20_alldata$Row) 
table(H4K20_alldata$Replicate, H4K20_alldata$Column)  

# add condition
# select only WT and KO
# remove weird points with StDev 0
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
    filter(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.CV.... != 0)
colnames(H4K20_alldata_int)[18] <- "H4K20m1_CV"

# check
table(H4K20_alldata_int$Replicate, H4K20_alldata_int$Row) 
table(H4K20_alldata_int$Replicate, H4K20_alldata_int$Column)  

table(H4K20_alldata_int$Condition, H4K20_alldata_int$Row) 
table(H4K20_alldata_int$Treatment, H4K20_alldata_int$Column)  

# filter single nuclei and round
H4K20_alldata_singl <- H4K20_alldata_int %>%
    filter(`All.Nuclei.Selected...Area.Nuclei.Area..µm..` < 200 & `All.Nuclei.Selected...Area.Nuclei.Area..µm..` > 60,
           `All.Nuclei.Selected...Area.Nuclei.Roundness` > 0.9)
nrow(H4K20_alldata_singl) # 21558
nrow(H4K20_alldata_int) # 36641
table(H4K20_alldata_singl$Condition) 

# remove mitotic cells with DAPI 
# only for DMSO

# scale each replicate separately between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

scale_DAPI_data <- function(n, treat) {
    repdata_treat <- H4K20_alldata_singl %>%
        filter(Replicate == n & Treatment == treat) 
    repdata_scale <- repdata_treat %>%
        mutate(logDAPI_sc = range01(logDAPI))
}
reps <- list(1,1,1,2,2,2)
treats <- list("DMSO", "thymidine", "RO3306", "DMSO", "thymidine", "RO3306")
H4K20_data_scaled_DAPI <- map2(reps, treats, scale_DAPI_data)

# check
table(H4K20_data_scaled_DAPI[[1]]["Replicate"])
table(H4K20_data_scaled_DAPI[[2]]["Replicate"])
table(H4K20_data_scaled_DAPI[[3]]["Replicate"])
table(H4K20_data_scaled_DAPI[[4]]["Replicate"])
table(H4K20_data_scaled_DAPI[[5]]["Replicate"])
table(H4K20_data_scaled_DAPI[[6]]["Replicate"])

table(H4K20_data_scaled_DAPI[[1]]["Treatment"])
table(H4K20_data_scaled_DAPI[[2]]["Treatment"])
table(H4K20_data_scaled_DAPI[[3]]["Treatment"])
table(H4K20_data_scaled_DAPI[[4]]["Treatment"])
table(H4K20_data_scaled_DAPI[[5]]["Treatment"])
table(H4K20_data_scaled_DAPI[[6]]["Treatment"])

# build data frame
H4K20_alldata_scaled_DAPI <- bind_rows(H4K20_data_scaled_DAPI)
table(H4K20_alldata_scaled_DAPI$Condition)
table(H4K20_alldata_scaled_DAPI$Treatment)

# mitotic cells have higher DAPI signal
hist(H4K20_alldata_scaled_DAPI$All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean)
hist(H4K20_alldata_scaled_DAPI$logDAPI)
hist(H4K20_alldata_scaled_DAPI$logDAPI_sc)
# we could try 0.7 as threshold

# remove mitotic cells only in DMSO
H4K20_alldata_scaled_DAPI_nomit <- H4K20_alldata_scaled_DAPI %>%
    filter(Treatment != "DMSO" | (Treatment == "DMSO" & logDAPI_sc < 0.7))
nrow(H4K20_alldata_scaled_DAPI_nomit) # 19986
nrow(H4K20_alldata_scaled_DAPI) # 21558

# check
table(H4K20_alldata_scaled_DAPI_nomit$Treatment)
table(H4K20_alldata_scaled_DAPI$Treatment)

# export data
save(H4K20_alldata_scaled_DAPI_nomit, file="H4K20_alldata_nomit.RData")

### Plot

# load data
load("H3K9_alldata_nomit.RData")
load("H3K27_alldata_nomit.RData")
load("H4K20_alldata_nomit.RData")

# select useful data
H4K20_alldata_int <- H4K20_alldata_scaled_DAPI_nomit %>%
    dplyr::select(Condition, Treatment, Replicate, log2H4K20me1, log2IntDen, H4K20m1_CV, log2CV) %>%
    mutate(Mark = rep("H4K20m1", nrow(H4K20_alldata_scaled_DAPI_nomit))) 
nrow(H4K20_alldata_int) # 19986

H4K20_alldata_int$Condition <- factor(H4K20_alldata_int$Condition, levels = c("WT", "KO1"))
H4K20_alldata_int$Treatment <- factor(H4K20_alldata_int$Treatment, levels = c("DMSO", "thymidine", "RO3306"))
H4K20_alldata_int$Replicate <- factor(H4K20_alldata_int$Replicate)

H3K9_alldata_int <- H3K9_alldata_scaled_DAPI_nomit %>%
    dplyr::select(Condition, Treatment, Replicate, log2H3K9me3, log2IntDen, H3K9m3_CV, log2CV) %>%
    mutate(Mark = rep("H3K9m3", nrow(H3K9_alldata_scaled_DAPI_nomit))) 
nrow(H3K9_alldata_int) # 26480
H3K9_alldata_int$Condition <- factor(H3K9_alldata_int$Condition, levels = c("WT", "KO1"))
H3K9_alldata_int$Treatment <- factor(H3K9_alldata_int$Treatment, levels = c("DMSO", "thymidine", "RO3306"))
H3K9_alldata_int$Replicate <- factor(H3K9_alldata_int$Replicate)

H3K27_alldata_int <- H3K27_alldata_scaled_DAPI_nomit %>%
    dplyr::select(Condition, Treatment, Replicate, log2H3K27me3, log2IntDen, H3K27m3_CV, log2CV) %>%
    mutate(Mark = rep("H3K27m3", nrow(H3K27_alldata_scaled_DAPI_nomit))) 
nrow(H3K27_alldata_int) # 38138
H3K27_alldata_int$Condition <- factor(H3K27_alldata_int$Condition, levels = c("WT", "KO1"))
H3K27_alldata_int$Treatment <- factor(H3K27_alldata_int$Treatment, levels = c("DMSO", "thymidine", "RO3306"))
H3K27_alldata_int$Replicate <- factor(H3K27_alldata_int$Replicate)

# put together all data
colnames(H4K20_alldata_int) <- c("Condition", "Treatment", "Replicate", "log2MeanInt", "log2IntDen", "CV", "log2CV", "Mark")
colnames(H3K9_alldata_int) <- c("Condition", "Treatment", "Replicate", "log2MeanInt", "log2IntDen", "CV", "log2CV", "Mark")
colnames(H3K27_alldata_int) <- c("Condition", "Treatment", "Replicate", "log2MeanInt", "log2IntDen", "CV", "log2CV", "Mark")

# put everything together
all_marks_together <- H4K20_alldata_int %>%
    rbind(H3K9_alldata_int) %>%
    rbind(H3K27_alldata_int)
all_marks_together$Mark <- factor(all_marks_together$Mark, levels = c("H4K20m1", "H3K9m3", "H3K27m3"))

# plot WT and KO
my_comp_KO1 <- list(c("WT", "KO1"))
all_marks_together$Condition <- factor(all_marks_together$Condition, levels = c("WT", "KO1"))

table(all_marks_together$Mark, all_marks_together$Condition)

# plot with violin
ggplot(all_marks_together, aes(x = Condition, y = log2MeanInt)) +
    geom_violin(aes(fill = Condition), col = NA, show.legend = T) +
    geom_boxplot(width = 0.25, alpha = 0, fatten = 2) +
    facet_grid(Treatment~Mark) +
    scale_fill_manual(values = c("WT" = "#aaafb0", "KO1" = "#e56e58")) +
    stat_compare_means(comparisons = my_comp_KO1, size = 3, tip.length = 0.01, braket.size = 0.01, label.y = 15,
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
ggsave("global_plots/log2MeanInt_allreps_violin_WTKO1.pdf", width = 4.5, height = 4.5)
