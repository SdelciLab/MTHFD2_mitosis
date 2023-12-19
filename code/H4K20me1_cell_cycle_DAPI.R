# libraries
library(dplyr)
library(ggplot2)
library(ggpubr) 
library(plotly)
library(ggforce)
library(viridis)
library(stringr)  
library(rstatix)

### get the data

load("../analyse_H4K20me1/H4K20_alldata_nomit.RData")

# leave the interesting columns only
colnames(H4K20_alldata_scaled_DAPI_nomit)

H4K20_alldata <- H4K20_alldata_scaled_DAPI_nomit[,c("Row", "Column", "Condition", "Replicate",
                                                    "All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.Mean",
                                                    "All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean",
                                                    "All.Nuclei.Selected...Area.Nuclei.Area..µm..")]

### Calculate integrated 
H4K20_alldata$integrated_488 = H4K20_alldata$All.Nuclei.Selected...Area.Nuclei.Area..µm.. * H4K20_alldata$All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.Mean
H4K20_alldata$integrated_DAPI = H4K20_alldata$All.Nuclei.Selected...Area.Nuclei.Area..µm.. * H4K20_alldata$All.Nuclei.Selected...Intensity.Nucleus.DAPI.Mean/500000

### Select interesting groups
H4K20_alldata_int <- H4K20_alldata %>%
    filter(Condition %in% c("WT", "KO1"))

### Set the gates (adjust according to the plots below)

move_all=0.25
{
    G1_gate_start =0.55+move_all
    G1_gate_end = 1.0+move_all
    
    G2M_gate_start = 1.65+move_all
    G2M_gate_end = 1.95+move_all
    
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
    
}

### DAPI test plot 
ggplot(H4K20_alldata_int, aes(x=integrated_DAPI)) +
    geom_histogram(aes(y=..density..), position = "identity", alpha=0.1) + 
    geom_density(alpha=0, size=.2, adjust= 0.5) +
    custom_gates+
    ggtitle(paste("integrated_DAPI - all wells"))+
    theme_bw()
ggsave("plots/H4K20me1_Integrated_DAPI_all_wells.pdf", device = "pdf")

## check replicates separated

H4K20_alldata_int_WT <- H4K20_alldata_int[H4K20_alldata_int$Condition == "WT",]
H4K20_alldata_int_KO1 <- H4K20_alldata_int[H4K20_alldata_int$Condition == "KO1",]

ggplot(H4K20_alldata_int_WT, aes(x=integrated_DAPI)) +
    geom_histogram(aes(y=..density..), position = "identity", alpha=0.1) + 
    facet_wrap(~factor(Replicate))+
    geom_density(alpha=0, size=.2, adjust= 0.5) +
    custom_gates+
    ggtitle(paste("integrated_DAPI - WT by rep"))+
    theme_bw()
ggsave("plots/H4K20me1_Integrated_DAPI_all_wells_WTbyrep.pdf", device = "pdf")

ggplot(H4K20_alldata_int_KO1, aes(x=integrated_DAPI)) +
    geom_histogram(aes(y=..density..), position = "identity", alpha=0.1) + 
    facet_wrap(~factor(Replicate))+
    geom_density(alpha=0, size=.2, adjust= 0.5) +
    custom_gates+
    ggtitle(paste("integrated_DAPI - KO1 by rep"))+
    theme_bw()
ggsave("plots/H4K20me1_Integrated_DAPI_all_wells_KO1byrep.pdf", device = "pdf")

## select rep 2 and 3 which show good cell cycle profile

H4K20_alldata_int_rep2 <- H4K20_alldata_int[H4K20_alldata_int$Replicate == 2,]
H4K20_alldata_int_rep3 <- H4K20_alldata_int[H4K20_alldata_int$Replicate == 3,]

## replicate 2

### Histogram DAPI Integrated signal from nuclei 

H4K20_alldata_int_rep2$Condition <- factor(H4K20_alldata_int_rep2$Condition,
                                           levels = c("WT", "KO1"))

ggplot(H4K20_alldata_int_rep2, aes(x=integrated_DAPI)) +
    coord_cartesian(xlim = c(.25, 3))+
    geom_density(aes(color =  Condition), alpha=0, size=.2, adjust= 0.35)+
    custom_gates+
    scale_color_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    ggtitle(paste("Integrated_DAPI - Treatments"))+
    theme_bw()
ggsave("plots/H4K20me1_rep2_Integrated_DAPI_Conditions.pdf", device = "pdf")

### Normalize profiles
H4K20_alldata_int_rep2 = H4K20_alldata_int_rep2 %>%
    mutate(norm_integrated_DAPI=
               case_when((Condition == "WT") ~ integrated_DAPI*1,
                         (Condition == "KO1") ~ integrated_DAPI*1.01
               )
    )

ggplot(H4K20_alldata_int_rep2, aes(x=norm_integrated_DAPI)) +
    coord_cartesian(xlim = c(.25, 3))+
    geom_density(aes(color =  Condition), alpha=0, size=.2, adjust= 0.35)+
    custom_gates+
    scale_color_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    ggtitle(paste("Integrated_DAPI - Treatments"))+
    theme_bw()
ggsave("plots/H4K20me1_rep2_Integrated_DAPI_Conditions_norm.pdf", device = "pdf")

### Histogram DAPI Integrated signal from nuclei facets

ggplot(H4K20_alldata_int_rep2, aes(x=norm_integrated_DAPI, fill = Condition)) +
    coord_cartesian(xlim = c(0.25, 3))+
    geom_histogram(aes(y=..density..),binwidth=0.02, position = "identity", alpha=0.1)+
    geom_density(aes(color =  Condition), alpha=0, size=.2, adjust= 0.35)+
    custom_gates+
    scale_color_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    scale_fill_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    facet_wrap(~factor(Condition), nrow = 2)+
    theme_bw() +
    theme(axis.title.y = element_text(size = 12), axis.text = element_text(size = 10), 
          axis.text.x = element_text(size = 12),
          strip.text = element_text(size = 12), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none")
ggsave("plots/H4K20me1_rep2_Integrated_DAPI_Conditions.pdf", device = "pdf",
       width = 2, height = 4)

### add the groups to the data table
H4K20_alldata_int_rep2_cycle = H4K20_alldata_int_rep2 %>%
    mutate(Phase_Group=                                                 
               case_when((between(norm_integrated_DAPI,G1_gate_start,G1_gate_end)) ~ "G1",
                         (between(norm_integrated_DAPI,S_gate_start,S_gate_end)) ~ "S",
                         (between(norm_integrated_DAPI,G2M_gate_start,G2M_gate_end)) ~ "G2M"))%>%
    mutate(Condition_Phase_Group=paste0(Condition,"\n",Phase_Group))

# select gates
move_all=0.15
{
    G1_gate_start =0.55+move_all
    G1_gate_end = 1.0+move_all
    
    G2M_gate_start = 1.65+move_all
    G2M_gate_end = 1.95+move_all
    
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
    
}

### Histogram DAPI Integrated signal from nuclei 
H4K20_alldata_int_rep3$Condition <- factor(H4K20_alldata_int_rep3$Condition,
                                           levels = c("WT", "KO1"))

ggplot(H4K20_alldata_int_rep3, aes(x=integrated_DAPI)) +
    coord_cartesian(xlim = c(.25, 3))+
    geom_density(aes(color =  Condition), alpha=0, size=.2, adjust= 0.35)+
    custom_gates+
    scale_color_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    ggtitle(paste("Integrated_DAPI - Treatments"))+
    theme_bw()
ggsave("plots/H4K20me1_rep3_Integrated_DAPI_Conditions.pdf", device = "pdf")

### Normalize profiles
H4K20_alldata_int_rep3 = H4K20_alldata_int_rep3 %>%
    mutate(norm_integrated_DAPI=
               case_when((Condition == "WT") ~ integrated_DAPI*1,
                         (Condition == "KO1") ~ integrated_DAPI*1.15
               )
    )

ggplot(H4K20_alldata_int_rep3, aes(x=norm_integrated_DAPI)) +
    coord_cartesian(xlim = c(.25, 3))+
    geom_density(aes(color =  Condition), alpha=0, size=.2, adjust= 0.35)+
    custom_gates+
    scale_color_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    ggtitle(paste("Integrated_DAPI - Treatments"))+
    theme_bw()
ggsave("plots/H4K20me1_rep3_Integrated_DAPI_Conditions_norm.pdf", device = "pdf")

### Histogram DAPI Integrated signal from nuclei facets

ggplot(H4K20_alldata_int_rep3, aes(x=norm_integrated_DAPI, fill = Condition)) +
    coord_cartesian(xlim = c(0.25, 3))+
    geom_histogram(aes(y=..density..),binwidth=0.02, position = "identity", alpha=0.1)+
    geom_density(aes(color =  Condition), alpha=0, size=.2, adjust= 0.35)+
    custom_gates+
    scale_color_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    scale_fill_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    facet_wrap(~factor(Condition))+#, scales = "free" )
    ggtitle(paste("Integrated_DAPI - Treatments - Facets"))+
    theme_bw()
ggsave("plots/H4K20me1_rep3_Integrated_DAPI_Conditions.pdf", device = "pdf")

### add the groups to the data table
H4K20_alldata_int_rep3_cycle = H4K20_alldata_int_rep3 %>%
    mutate(Phase_Group=                                                 
               case_when((between(norm_integrated_DAPI,G1_gate_start,G1_gate_end)) ~ "G1",
                         (between(norm_integrated_DAPI,S_gate_start,S_gate_end)) ~ "S",
                         (between(norm_integrated_DAPI,G2M_gate_start,G2M_gate_end)) ~ "G2M"))%>%
    mutate(Condition_Phase_Group=paste0(Condition,"\n",Phase_Group))

## add 2 replicates together
H4K20_alldata_int_both_rep <- H4K20_alldata_int_rep2_cycle %>%
    rbind(H4K20_alldata_int_rep3_cycle)

# remove NA
H4K20_alldata_int_cycle_noNA <- na.omit(H4K20_alldata_int_both_rep)

# factor
H4K20_alldata_int_cycle_noNA$Condition_Phase_Group <- factor(H4K20_alldata_int_cycle_noNA$Condition_Phase_Group,
                                                             levels = c("WT\nG1", "KO1\nG1",
                                                                        "WT\nS", "KO1\nS",
                                                                        "WT\nG2M", "KO1\nG2M"))

### Boxplot 488 mean intensity by Phases and Treatments - Treatments grouped
my_comp <- list(c("WT\nG1", "KO1\nG1"), c("WT\nS", "KO1\nS"), c("WT\nG2M", "KO1\nG2M"))

# plot mean intensity
ggplot(H4K20_alldata_int_cycle_noNA, aes(x=Condition_Phase_Group, y =log2(All.Nuclei.Selected...Intensity.Nucleus.H4K20m1.Mean),
                                         fill = Condition)) +
    geom_boxplot(width=0.75)+
    scale_color_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    scale_fill_manual(values = c("WT" = "#aaafb0",  "KO1" = "#e56e58")) +
    stat_compare_means(
        comparisons = my_comp, 
        label = "p.format", size = 2, method = "wilcox"
    ) +
    labs(x = "", y="log2(H4K20me1 mean intensity)") +
    theme_classic() +
    theme(legend.position="none") 
ggsave("plots/H4K20me1_boxplot_conditions_cellphases_rep23_meanint_nice.pdf", device = "pdf", width = 4, height = 3)

# numbers
table(H4K20_alldata_int_cycle_noNA$Condition_Phase_Group)