# libraries
library(dplyr)
library(ggplot2)
library(ggpubr) 
library(plotly) 
library(ggforce)
library(viridis)
library(stringr)  
library(rstatix)
library(data.table)

# get files
CENPA_nuclei <-  read.csv("new pipeline_CREST_488int_2/230926_1 Natalia_s__2023-10-16T13_13_52-Measurement 8 (CENPA)/Evaluation8/Objects_Population - Nuclei C Selected Selected.txt", 
                          header = TRUE, sep="\t",skip=9) %>%
    select(where(function(x) any(!is.na(x))))

CENPA_spots <- read.csv("new pipeline_CREST_488int_2/230926_1 Natalia_s__2023-10-16T13_13_52-Measurement 8 (CENPA)/Evaluation8/Objects_Population - Spots_488.txt", 
                        header = TRUE, sep="\t",skip=9) %>%
    select(where(function(x) any(!is.na(x))))

CREST_CENPA_spots <- read.csv("new pipeline_CREST_488int_2/230926_1 Natalia_s__2023-10-16T13_13_52-Measurement 8 (CENPA)/Evaluation8/Objects_Population - Spots_647.txt", 
                              header = TRUE, sep="\t",skip=9) %>%
    select(where(function(x) any(!is.na(x))))

######### prepare nuclei data ########

# leave the interesting columns only
nuclei_data = CENPA_nuclei[,c(1,2,4,5,12:14,16:36)]

### Change column names to make them more handy
{
    names(nuclei_data)[1]='Row'
    names(nuclei_data)[2]='Column'
    names(nuclei_data)[3]='Field'
    names(nuclei_data)[4]='Nuclei_Object_no'
    
    names(nuclei_data)[5]='Nuclei_Area'
    names(nuclei_data)[6]='Nuclei_Roundness'
    names(nuclei_data)[7]='Nuclei_Ratio_wxL'
    
    names(nuclei_data)[8]='Nuclei_Intensity_Nuc_488'
    names(nuclei_data)[9]='Nuclei_Ring_Intensity_Nuc_488'
    
    names(nuclei_data)[10]='Nuclei_Intensity_Nuc_DAPI'
    names(nuclei_data)[11]='Nuclei_Ring_Intensity_Nuc_DAPI'
    
    names(nuclei_data)[12]='Nuclei_Intensity_Nuc_555'
    names(nuclei_data)[13]='Nuclei_Ring_Intensity_Nuc_555'
    
    names(nuclei_data)[14]='Nuclei_Intensity_Nuc_647'
    names(nuclei_data)[15]='Nuclei_Ring_Intensity_Nuc_647'
    
    names(nuclei_data)[16]='Nuclei_Total_Spot_Area_488'
    names(nuclei_data)[17]='Nuclei_Rel_Spot_Intensity_488'
    names(nuclei_data)[18]='Nuclei_Number_of_Spots_488'
    names(nuclei_data)[19]='Nuclei_Number_of_Spots_x_Area_488'
    
    names(nuclei_data)[20]='Nuclei_Total_Spot_Area_555'
    names(nuclei_data)[21]='Nuclei_Rel_Spot_Intensity_555'
    names(nuclei_data)[22]='Nuclei_Number_of_Spots_555'
    names(nuclei_data)[23]='Nuclei_Number_of_Spots_x_Area_555'
    
    names(nuclei_data)[24]='Nuclei_Total_Spot_Area_647'
    names(nuclei_data)[25]='Nuclei_Rel_Spot_Intensity_647'
    names(nuclei_data)[26]='Nuclei_Number_of_Spots_647'
    names(nuclei_data)[27]='Nuclei_Number_of_Spots_x_Area_647'
    
    names(nuclei_data)[28]='Intensity_488_Spots_647'
}

######### prepare 488 spots data ########

# leave the interesting columns only
CENPA_spots_int = CENPA_spots[,c(1,2,4,5,11:13,16,19,21:23)]

{
    names(CENPA_spots_int)[1]='Row'
    names(CENPA_spots_int)[2]='Column'
    names(CENPA_spots_int)[3]='Field'
    names(CENPA_spots_int)[4]='Spot_Object_no'
    names(CENPA_spots_int)[5]='Spot_Rel_Intensity'
    names(CENPA_spots_int)[6]='Spot_Corr_Intensity'
    names(CENPA_spots_int)[7]='Spot_UnCorr_Intensity'
    names(CENPA_spots_int)[8]='Spot_Area'
    names(CENPA_spots_int)[9]='Nuclei_Object_no'
    names(CENPA_spots_int)[10]='Spot_488_Nearest-647_distance'
    names(CENPA_spots_int)[11]='Spot_488_Nearest-647_Object_no'
    names(CENPA_spots_int)[12]='Spot_488_Nearest-647_Overlap'
}

spots_data = CENPA_spots_int

######### prepare CREST spots data ########
CREST_spots_int = CREST_CENPA_spots[,c(1,2,4,5,11:13,16,19,21:24)]

{
    names(CREST_spots_int)[1]='Row'
    names(CREST_spots_int)[2]='Column'
    names(CREST_spots_int)[3]='Field'
    names(CREST_spots_int)[4]='Spot_Object_no'
    names(CREST_spots_int)[5]='Spot_Rel_Intensity'
    names(CREST_spots_int)[6]='Spot_Corr_Intensity'
    names(CREST_spots_int)[7]='Spot_UnCorr_Intensity'
    names(CREST_spots_int)[8]='Spot_Area'
    names(CREST_spots_int)[9]='Nuclei_Object_no'
    names(CREST_spots_int)[10]='Spot_647_Nearest-488_distance'
    names(CREST_spots_int)[11]='Spot_647_Nearest-488_Object_no'
    names(CREST_spots_int)[12]='Spot_647_Nearest-488_Overlap'
    names(CREST_spots_int)[13]='Spot_647_488_Mean_Intensity'
}

### Remove background
nuclei_data$dif_488 = nuclei_data$Nuclei_Intensity_Nuc_488 - nuclei_data$Nuclei_Ring_Intensity_Nuc_488
nuclei_data$dif_DAPI = nuclei_data$Nuclei_Intensity_Nuc_DAPI - nuclei_data$Nuclei_Ring_Intensity_Nuc_DAPI
nuclei_data$dif_555 = nuclei_data$Nuclei_Intensity_Nuc_555 - nuclei_data$Nuclei_Ring_Intensity_Nuc_555
nuclei_data$dif_647 = nuclei_data$Nuclei_Intensity_Nuc_647 - nuclei_data$Nuclei_Ring_Intensity_Nuc_647

### Calculate integrated 
nuclei_data$integrated_488 = nuclei_data$Nuclei_Area * nuclei_data$dif_488/1000000
nuclei_data$integrated_DAPI = nuclei_data$Nuclei_Area * nuclei_data$dif_DAPI/750000
nuclei_data$integrated_555 = nuclei_data$Nuclei_Area * nuclei_data$dif_555/1000000
nuclei_data$integrated_647 = nuclei_data$Nuclei_Area * nuclei_data$dif_647/1000000

### nuclei: Assign a name to the groups
nuclei_data_cond <- nuclei_data %>%
    mutate(Condition=
               case_when((Row == 2 & (between(Column,2,11))) ~ "WT",
                         (Row == 3 & (between(Column,2,11))) ~ "KO1",
               ),
           Antibodies_488=
               case_when((between(Row, 2,5) & (between(Column,2,3))) ~ "CENPA_488",
                         (between(Row, 2,5) & (between(Column,6,7))) ~ "CENPA_488",
                         (between(Row, 2,5) & (between(Column,10,11))) ~ "no_antibody"
               ),
           Antibodies_647=
               case_when((between(Row, 2,5) & (between(Column,2,9))) ~ "CREST_647", 
                         (between(Row, 2,5) & (between(Column,10,11))) ~ "no_antibody"
               ),
           Nuclei_ID=paste0(Row,Column,Field,Nuclei_Object_no)
    )%>%
    filter(Condition %in% c("WT","KO1"),
           Antibodies_488 %in% c("CENPA_488"),
           Antibodies_647 == "CREST_647") 

### 488 spots: Assign a name to the groups
spots_data_cond <- spots_data %>%
    mutate(Condition=
               case_when((Row == 2 & (between(Column,2,11))) ~ "WT",
                         (Row == 3 & (between(Column,2,11))) ~ "KO1"
               ),
           Antibodies_488=
               case_when((between(Row, 2,5) & (between(Column,2,3))) ~ "CENPA_488", 
                         (between(Row, 2,5) & (between(Column,6,7))) ~ "CENPA_488",
                         (between(Row, 2,5) & (between(Column,10,11))) ~ "no_antibody"
               ),
           Antibodies_647=
               case_when((between(Row, 2,5) & (between(Column,2,9))) ~ "CREST_647", 
                         (between(Row, 2,5) & (between(Column,10,11))) ~ "no_antibody"
               ),
           Nuclei_ID=paste0(Row,Column,Field,Nuclei_Object_no)
    )%>%
    filter(Condition %in% c("WT","KO1"),
           Antibodies_488 %in% c("CENPA_488"),
           Antibodies_647 == "CREST_647") 

### 647 spots: Assign a name to the groups
CREST_spots_cond <- CREST_spots_int %>%
    mutate(Condition=
               case_when((Row == 2 & (between(Column,2,11))) ~ "WT",
                         (Row == 3 & (between(Column,2,11))) ~ "KO1"
               ),
           Antibodies_488=
               case_when((between(Row, 2,5) & (between(Column,2,3))) ~ "CENPA_488",
                         (between(Row, 2,5) & (between(Column,6,7))) ~ "CENPA_488",
                         (between(Row, 2,5) & (between(Column,10,11))) ~ "no_antibody"
               ),
           Antibodies_647=
               case_when((between(Row, 2,5) & (between(Column,2,9))) ~ "CREST_647", 
                         (between(Row, 2,5) & (between(Column,10,11))) ~ "no_antibody"
               ),
           Nuclei_ID=paste0(Row,Column,Field,Nuclei_Object_no)
    )%>%
    filter(Condition %in% c("WT","KO1"),
           Antibodies_488 %in% c("CENPA_488"),
           Antibodies_647 == "CREST_647") 

########### TRY TO SEPARATE PHASES

### DAPI test plot 
nuclei_data_test  = nuclei_data_cond %>% as.data.table()
nuclei_data_test[, median_dapi := median(integrated_DAPI), by = Condition]
nuclei_data_test[, median_norm_dapi := integrated_DAPI - median_dapi, by = .(Row, Column,Field,Nuclei_Object_no)]

ggplot(nuclei_data, aes(x=integrated_DAPI)) +
    geom_histogram(aes(y=..density..),binwidth=0.02, position = "identity", alpha=0.1)+ 
    geom_density(alpha=0, size=.2, adjust= 0.5)+
    ggtitle(paste("integrated_DAPI - all wells"))+
    theme_bw()
ggsave("plots/testplot - Integrated_DAPI - all wells.png")

### DAPI test plot Integrated signal from nuclei
ggplot(nuclei_data_test, aes(x=median_norm_dapi)) +
    geom_density(aes(color =  Condition), alpha=0, size=.2, adjust= 0.5)+
    ggtitle(paste("Integrated_DAPI - Treatments?"))+
    theme_bw()
ggsave("plots/plot_1_raw - Integrated_DAPI - Condition.png")

### Normalize profiles
all_norm=1.5
nuclei_data_cond = nuclei_data_cond %>%
    mutate(norm_integrated_DAPI=
               case_when((Condition == "WT") ~ integrated_DAPI*1*all_norm,
                         (Condition == "KO1") ~ integrated_DAPI*1*all_norm)
    )

### Set the gates (adjust according to the plots below)
move_all=0
{
    G1_gate_start =0.6+move_all
    G1_gate_end = 0.95+move_all
    
    G2M_gate_start = 1.25+move_all
    G2M_gate_end = 1.6+move_all
    
    S_gate_start = G1_gate_end+0.01
    S_gate_end = G2M_gate_start-0.01
    
    
    
    custom_gates=c(geom_segment(aes(x = G1_gate_start, xend = G1_gate_end, y = 1, yend = 1)),
                   geom_segment(aes(x = G1_gate_start, xend = G1_gate_start, y = 1.05, yend = .95)),
                   geom_segment(aes(x = G1_gate_end, xend = G1_gate_end, y = 1.05, yend = .95)),
                   
                   geom_segment(aes(x = S_gate_start, xend = S_gate_end, y = 1, yend = 1)),
                   geom_segment(aes(x = S_gate_start, xend = S_gate_start, y = 1.05, yend = .95)),
                   geom_segment(aes(x = S_gate_end, xend = S_gate_end, y = 1.05, yend = .95)),
                   
                   geom_segment(aes(x = G2M_gate_start, xend = G2M_gate_end, y = 1, yend = 1)),
                   geom_segment(aes(x = G2M_gate_start, xend = G2M_gate_start, y = 1.05, yend = .95)),
                   geom_segment(aes(x = G2M_gate_end, xend = G2M_gate_end, y = 1.05, yend = .95)))
    
    custom_gates_vline = geom_vline(xintercept = c(G1_gate_start,G1_gate_end,S_gate_start,S_gate_end,G2M_gate_start,G2M_gate_end),linetype="dotted")
    
}

### DAPI test plot normalized DAPI Integrated signal from nuclei
ggplot(nuclei_data_cond, aes(x=norm_integrated_DAPI)) +
    coord_cartesian(xlim = c(.25, 2))+
    geom_density(aes(color =  Condition), alpha=0, size=.2, adjust= 0.5)+
    custom_gates+
    scale_x_continuous(breaks=seq(.25,2,by=0.1))+
    ggtitle(paste("Integrated_DAPI - Condition"))+
    theme_bw()
ggsave("plots/plot_1 - Integrated_DAPI - Conditions.png")

###  DAPI test plot normalized DAPI Integrated signal from nuclei facets
ggplot(nuclei_data_cond, aes(x=norm_integrated_DAPI, fill = Condition)) +
    coord_cartesian(xlim = c(0.25, 2))+
    geom_histogram(aes(y=..density..),binwidth=0.02, position = "identity", alpha=0.1)+
    geom_density(aes(color =  Condition), alpha=0, size=.2, adjust= 0.55)+
    custom_gates+
    facet_wrap(Antibodies_488~Condition)+
    ggtitle(paste("Integrated_DAPI - Condition - Facets"))+
    theme_bw()
ggsave("plots/plot_2 - Integrated_DAPI - Conditions - Facets.png")

### add the groups to the data table
nuclei_data_cond_phases = nuclei_data_cond %>%
    mutate(Phase_Group=                                                 
               case_when((between(norm_integrated_DAPI,G1_gate_start,G1_gate_end)) ~ "G1",
                         (between(norm_integrated_DAPI,S_gate_start,S_gate_end)) ~ "S",
                         (between(norm_integrated_DAPI,G2M_gate_start,G2M_gate_end)) ~ "G2M"))%>%
    mutate(Condition_Phase_Group=paste0(Condition,"\n",Phase_Group))

# add phases to CREST spots
CREST_spots_cond_phases  <- merge(CREST_spots_cond,nuclei_data_cond_phases)

# focus only WT and KO
CREST_spots_cond_phases_WTKO <- CREST_spots_cond_phases %>%
    filter(Condition %in% c("WT", "KO1")) %>%
    filter(!is.na(Phase_Group))

# factor
CREST_spots_cond_phases_WTKO$Condition <- factor(CREST_spots_cond_phases_WTKO$Condition, 
                                                 levels = c("WT", "KO1"))
CREST_spots_cond_phases_WTKO$Condition_Phase_Group <- factor(CREST_spots_cond_phases_WTKO$Condition_Phase_Group, 
                                                             levels = c("WT\nG1", "KO1\nG1",
                                                                        "WT\nS", "KO1\nS",
                                                                        "WT\nG2M", "KO1\nG2M"))
# check numbers
table(CREST_spots_cond_phases_WTKO$Condition_Phase_Group, 
      CREST_spots_cond_phases_WTKO$Antibodies_488)

# plot CENPA intensity in CREST spots
my_comp <- list(c("WT\nG1", "KO1\nG1"), c("WT\nS", "KO1\nS"), c("WT\nG2M", "KO1\nG2M"))
ggplot(CREST_spots_cond_phases_WTKO, aes(x=Condition_Phase_Group, y=log2(Intensity_488_Spots_647), fill = Condition)) +
    geom_boxplot(width = 0.75) +
    facet_grid(~Antibodies_488) +
    labs(x="", y="log2(CENP mean intensity) \nin CREST spots)") +
    scale_fill_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    stat_compare_means(
        comparisons = my_comp, 
        label = "p.format", size = 2, method = "wilcox"
    ) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("plots/CREST_spots_488int_individual_phases_nice.pdf", device = "pdf", width = 5, height = 3)

# add phases to CENPA spots
spots_data_cond_phases  <- merge(spots_data_cond,nuclei_data_cond_phases)

# focus only WT and KO
spots_data_cond_phases_WTKO <- spots_data_cond_phases %>%
    filter(Condition %in% c("WT", "KO1")) %>%
    filter(!is.na(Phase_Group))

# factor
spots_data_cond_phases_WTKO$Condition <- factor(spots_data_cond_phases_WTKO$Condition, 
                                                levels = c("WT", "KO1"))
spots_data_cond_phases_WTKO$Condition_Phase_Group <- factor(spots_data_cond_phases_WTKO$Condition_Phase_Group, 
                                                            levels = c("WT\nG1", "KO1\nG1",
                                                                       "WT\nS", "KO1\nS",
                                                                       "WT\nG2M", "KO1\nG2M"))
# check numbers
table(spots_data_cond_phases_WTKO$Condition_Phase_Group,
      spots_data_cond_phases_WTKO$Antibodies_488)

# keep only CENPA spots close to CREST (distance = 0)
spots_data_cond_phases_WTKO_dist <- spots_data_cond_phases_WTKO %>%
    filter(`Spot_488_Nearest-647_distance` == 0)

# plot overlap between CENPA and CREST in spots with distance = 0
my_comp <- list(c("WT\nG1", "KO1\nG1"), c("WT\nS", "KO1\nS"), c("WT\nG2M", "KO1\nG2M"))
ggplot(spots_data_cond_phases_WTKO_dist, aes(x=Condition_Phase_Group, y=`Spot_488_Nearest-647_Overlap`, fill = Condition)) +
    geom_boxplot(width = 0.75) +
    facet_grid(~Antibodies_488) +
    labs(x="", y="Overlap of CREST spots \n within CENP spots (%)") +
    scale_fill_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    stat_compare_means(
        comparisons = my_comp, 
        label = "p.format", size = 2, method = "wilcox"
    ) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none"
    )
ggsave("plots/CENPA_spots_overlapCREST_phases_nice.pdf", device = "pdf", width = 5, height = 3)
