# libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr) # for ggarrange

# 1. Get input data U2OS ------------------------------------------------------------------------

# open file
fuccidata=read.csv("U2OS/FUCCI MTHFD2 647_220603 Ev3.txt", header = TRUE,sep="\t",skip=9)

# get the interesting columns
fuccidata= fuccidata[,20:39] # 

# change the name of the columns to make them more handy
#colnames(fuccidata) to see a list of the columns names

names(fuccidata)[1]='Object n sel sel'
names(fuccidata)[2]='Nuc_Area'
names(fuccidata)[3]='Nuc_Roundness'
names(fuccidata)[4]='Cyt_Area'
names(fuccidata)[5]='Cyt_Roundness'
names(fuccidata)[6]='Cell_Area'
names(fuccidata)[7]='Cell_Roundness'
names(fuccidata)[8]='NAdivCAxCR'
names(fuccidata)[9]='nuc_int_488'
names(fuccidata)[10]='nuc_int_Turq'
names(fuccidata)[11]='nuc_int_546'
names(fuccidata)[12]='nuc_int_647'
names(fuccidata)[13]='cyt_int_488'
names(fuccidata)[14]='cyt_int_Turq'
names(fuccidata)[15]='cyt_int_546'
names(fuccidata)[16]='cyt_int_647'
names(fuccidata)[17]='nuc_intR_488'
names(fuccidata)[18]='nuc_intR_Turq'
names(fuccidata)[19]='nuc_intR_546'
names(fuccidata)[20]='nuc_intR_647'
            
#substract green signal from around the nucleus 
fuccidata$dif488=fuccidata$nuc_int_488 - fuccidata$nuc_intR_488
fuccidata$difTurq=fuccidata$nuc_int_Turq - fuccidata$nuc_intR_Turq
fuccidata$dif546=fuccidata$nuc_int_546 - fuccidata$nuc_intR_546

# compensate bleed-through
fuccidata$dif488_int_btcomp = fuccidata$dif488-(0.571*fuccidata$difTurq-56.59) 

# 2. Obtain plots U2OS -------------------------------------------------------------------------------------------

# remove negative values, get logs and scale
fuccidata2 <- fuccidata %>%
  filter(difTurq > 0 & dif488_int_btcomp > 0 & dif546 > 0,
         nuc_int_647 > 0 & cyt_int_647 > 0) %>%
  mutate(logdifTurq = log10(difTurq),
         logdif488_comp = log10(dif488_int_btcomp),
         logdif546 = log10(dif546),
         lognuc_647 = log10(nuc_int_647),
         logcyt_647 = log10(cyt_int_647)) %>% 
  mutate(diflog_Turq_sc = (logdifTurq - min(logdifTurq)) / (max(logdifTurq) - min(logdifTurq)),
         diflog_488_sc = (logdif488_comp - min(logdif488_comp)) / (max(logdif488_comp) - min(logdif488_comp)),
         diflog_546_sc = (logdif546 - min(logdif546)) / (max(logdif546) - min(logdif546)),
         diflog_nuc647_sc = (lognuc_647 - min(lognuc_647)) / (max(lognuc_647) - min(lognuc_647)),
         diflog_cyt647_sc = (logcyt_647 - min(logcyt_647)) / (max(logcyt_647) - min(logcyt_647)))

# plot MTHFD2 with red - nucleus
ggplot(fuccidata2, aes(x=difTurq, y=dif488_int_btcomp)) + 
  geom_point(aes(colour = nuc_int_647), size = 1.5) +
  ggtitle(paste0("MTHFD2 - Nucleus"))+
  scale_colour_gradient(low = "#FFFFFF", high = "#54002A", na.value = NA)+ 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic() +
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  scale_x_continuous(trans='log10', limits = c(30,14000)) +
  scale_y_continuous(trans='log10', limits = c(1, 14000))
ggsave("U2OS/plots/MTHFD2_nuc_plot_NPL.pdf", device = "pdf", width = 4.5, height = 4.5)


# plot MTHFD2 with red - cytosol
ggplot(fuccidata2, aes(x=difTurq, y=dif488_int_btcomp)) + 
  geom_point(aes(colour = cyt_int_647), size = 1.5) +
  theme_bw()+
  ggtitle(paste0("MTHFD2 - Cytosol"))+
  scale_colour_gradient(low = "#FFFFFF", high = "#54002A", na.value = NA)+ 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic() +
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  scale_x_continuous(trans='log10', limits = c(30,14000)) +
  scale_y_continuous(trans='log10', limits = c(1, 14000))
ggsave("U2OS/plots/MTHFD2_cyt_plot_NPL.pdf", device = "pdf", width = 4.5, height = 4.5)

# plot FUCCI
ggplot(fuccidata2, aes(x=difTurq, y=dif488_int_btcomp)) + 
  geom_point(aes(colour = dif546), size = 1.5, alpha = 1) +
  ggtitle(paste0("FUCCIplot"))+
  scale_color_gradient(low = "#DAD5D0", high = "#9E4F00", na.value = NA,
                       limits= c(0.457,600), oob = scales::squish) + 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "mKO2")+
  theme_classic() +
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
  ) +
  scale_x_continuous(trans='log10', limits = c(30,14000)) +
  scale_y_continuous(trans='log10', limits = c(1, 14000))
ggsave("U2OS/plots/FUCCIplot_NPL.pdf", device = "pdf", width = 4.6, height = 3.8)

# plot highest values nucleus
fuccidata2_levels <- fuccidata2 %>%
  mutate(level_nuc = ifelse(nuc_int_647 > quantile(nuc_int_647, 0.95), "high",
                        ifelse(nuc_int_647 < quantile(nuc_int_647, 0.05), "low", "medium")))

ggplot(fuccidata2_levels, aes(color = nuc_int_647)) + 
  geom_point(data = fuccidata2_levels %>% subset(level_nuc == "high"), aes(x=difTurq, y=dif488_int_btcomp),
             size = 1.5) +
  ggtitle(paste0("MTHFD2 - Nucleus"))+
  scale_colour_gradient(low = "#FFFFFF", high = "#54002A", na.value = NA,
                        limits = c(845.701,17961.300))+ 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic()+
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  scale_x_continuous(trans='log10', limits = c(30,14000)) +
  scale_y_continuous(trans='log10', limits = c(1, 14000))
ggsave("U2OS/plots/MTHFD2_nuc_plot_high_NPL.pdf", device = "pdf", width = 4.5, height = 4.5)

# plot highest values nucleus density
ggplot(data = fuccidata2_levels %>% subset(level_nuc == "high")) + 
  stat_density_2d(aes(x=difTurq, y=dif488_int_btcomp, fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  ggtitle(paste0("MTHFD2 - Nucleus"))+
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic()+
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  ) +
  scale_x_continuous(trans='log10', limits = c(30,14000)) +
  scale_y_continuous(trans='log10', limits = c(1, 14000))
ggsave("U2OS/plots/MTHFD2_nuc_plot_high_NPL_5top_density.pdf", device = "pdf", width = 4.4, height = 3.65)

# plot highest values cytosol
fuccidata2_levels <- fuccidata2 %>%
  mutate(level_cyt = ifelse(cyt_int_647 > quantile(cyt_int_647, 0.95), "high",
                        ifelse(cyt_int_647 < quantile(cyt_int_647, 0.05), "low", "medium")))

ggplot(fuccidata2_levels, aes(color = cyt_int_647)) + 
  geom_point(data = fuccidata2_levels %>% subset(level_cyt == "high"), aes(x=difTurq, y=dif488_int_btcomp),
             size = 1.5) +
  ggtitle(paste0("MTHFD2 - Cytosol"))+
  scale_colour_gradient(low = "#FFFFFF", high = "#54002A", na.value = NA,
                        limits = c(816.344,10613.300))+ 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic()+
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  scale_x_continuous(trans='log10', limits = c(30,14000)) +
  scale_y_continuous(trans='log10', limits = c(1, 14000))
ggsave("U2OS/plots/MTHFD2_cyt_plot_high_NPL.pdf", device = "pdf", width = 4.5, height = 4.5)

# plot highest values cytosol density
ggplot(data = fuccidata2_levels %>% subset(level_cyt == "high")) + 
  stat_density_2d(aes(x=difTurq, y=dif488_int_btcomp, fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  ggtitle(paste0("MTHFD2 - Cytosol"))+
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic()+
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  ) +
  scale_x_continuous(trans='log10', limits = c(30,14000)) +
  scale_y_continuous(trans='log10', limits = c(1, 14000))
ggsave("U2OS/plots/MTHFD2_cyt_plot_high_NPL_5top_density.pdf", device = "pdf", width = 4.4, height = 3.65)

# 3. Get input data MCF7-------------------------------------------------------------------------------------

# open file
fuccidata=read.table("MCF7/MTHFD2 - 230328_1 Meas 3 Ev1 - other.txt", header = TRUE,sep="\t")

# this file was already processed externally to get and modify the interesting columns

# 4. Obtain plots MCF7 -------------------------------------------------------------------------------------------

# remove negative values, get logs and scale
fuccidata2 <- fuccidata %>%
  filter(difTurq > 0 & dif488_int_btcomp > 0 & dif546 > 0,
         nuc_int_647 > 0 & cyt_int_647 > 0) %>%
  mutate(logdifTurq = log10(difTurq),
         logdif488_comp = log10(dif488_int_btcomp),
         logdif546 = log10(dif546),
         lognuc_647 = log10(nuc_int_647),
         logcyt_647 = log10(cyt_int_647)) %>% 
  mutate(diflog_Turq_sc = (logdifTurq - min(logdifTurq)) / (max(logdifTurq) - min(logdifTurq)),
         diflog_488_sc = (logdif488_comp - min(logdif488_comp)) / (max(logdif488_comp) - min(logdif488_comp)),
         diflog_546_sc = (logdif546 - min(logdif546)) / (max(logdif546) - min(logdif546)),
         diflog_nuc647_sc = (lognuc_647 - min(lognuc_647)) / (max(lognuc_647) - min(lognuc_647)),
         diflog_cyt647_sc = (logcyt_647 - min(logcyt_647)) / (max(logcyt_647) - min(logcyt_647)))

# plot MTHFD2 with red - nucleus
ggplot(fuccidata2, aes(x=difTurq, y=dif488_int_btcomp)) + 
  geom_point(aes(colour = nuc_int_647), size = 1.5) +
  ggtitle(paste0("MTHFD2 - Nucleus"))+
  scale_colour_gradient(low = "#FFFFFF", high = "#54002A", na.value = NA)+ 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic() +
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  scale_x_continuous(trans='log10', limits = c(30,21000)) +
  scale_y_continuous(trans='log10', limits = c(1,21000))
ggsave("MCF7/plots/MTHFD2_nuc_plot_NPL.pdf", device = "pdf", width = 4.5, height = 4.5)


# plot MTHFD2 with red - cytosol
ggplot(fuccidata2, aes(x=difTurq, y=dif488_int_btcomp)) + 
  geom_point(aes(colour = cyt_int_647), size = 1.5) +
  theme_bw()+
  ggtitle(paste0("MTHFD2 - Cytosol"))+
  scale_colour_gradient(low = "#FFFFFF", high = "#54002A", na.value = NA)+ 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic() +
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  scale_x_continuous(trans='log10', limits = c(30,21000)) +
  scale_y_continuous(trans='log10', limits = c(1,21000))
ggsave("MCF7/plots/MTHFD2_cyt_plot_NPL.pdf", device = "pdf", width = 4.5, height = 4.5)

# plot FUCCI
ggplot(fuccidata2, aes(x=difTurq, y=dif488_int_btcomp)) + 
  geom_point(aes(colour = dif546), size = 1.5, alpha = 1) +
  ggtitle(paste0("FUCCIplot"))+
  scale_color_gradient(low = "#DAD5D0", high = "#9E4F00", na.value = NA,
                       limits= c(0.08,1800), oob = scales::squish)+ 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "mKO2")+
  theme_classic() +
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    #legend.position = "none"
  ) +
  scale_x_continuous(trans='log10', limits = c(30,21000)) +
  scale_y_continuous(trans='log10', limits = c(1,21000))
ggsave("MCF7/plots/FUCCIplot_NPL.pdf", device = "pdf", width = 4.6, height = 3.8)

# plot highest values nucleus
fuccidata2_levels <- fuccidata2 %>%
  mutate(level_nuc = ifelse(nuc_int_647 > quantile(nuc_int_647, 0.95), "high",
                            ifelse(nuc_int_647 < quantile(nuc_int_647, 0.05), "low", "medium")))

ggplot(fuccidata2_levels, aes(color = nuc_int_647)) + 
  geom_point(data = fuccidata2_levels %>% subset(level_nuc == "high"), aes(x=difTurq, y=dif488_int_btcomp),
             size = 1.5) +
  ggtitle(paste0("MTHFD2 - Nucleus"))+
  scale_colour_gradient(low = "#FFFFFF", high = "#54002A", na.value = NA,
                        limits = c(845.701,17961.300))+ 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic()+
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  scale_x_continuous(trans='log10', limits = c(30,21000)) +
  scale_y_continuous(trans='log10', limits = c(1,21000))
ggsave("MCF7/plots/MTHFD2_nuc_plot_high_NPL.pdf", device = "pdf", width = 4.5, height = 4.5)

# plot highest values nucleus density
ggplot(data = fuccidata2_levels %>% subset(level_nuc == "high")) + 
  stat_density_2d(aes(x=difTurq, y=dif488_int_btcomp, fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  ggtitle(paste0("MTHFD2 - Nucleus"))+
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic()+
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  ) +
  scale_x_continuous(trans='log10', limits = c(30,21000)) +
  scale_y_continuous(trans='log10', limits = c(1,21000))
ggsave("MCF7/plots/MTHFD2_nuc_plot_high_NPL_5top_density.pdf", device = "pdf", width = 4.4, height = 3.65)

# plot highest values cytosol
fuccidata2_levels <- fuccidata2 %>%
  mutate(level_cyt = ifelse(cyt_int_647 > quantile(cyt_int_647, 0.95), "high",
                            ifelse(cyt_int_647 < quantile(cyt_int_647, 0.05), "low", "medium")))

ggplot(fuccidata2_levels, aes(color = cyt_int_647)) + 
  geom_point(data = fuccidata2_levels %>% subset(level_cyt == "high"), aes(x=difTurq, y=dif488_int_btcomp),
             size = 1.5) +
  ggtitle(paste0("MTHFD2 - Cytosol"))+
  scale_colour_gradient(low = "#FFFFFF", high = "#54002A", na.value = NA,
                        limits = c(816.344,10613.300))+ 
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic()+
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  scale_x_continuous(trans='log10', limits = c(30,21000)) +
  scale_y_continuous(trans='log10', limits = c(1,21000))
ggsave("MCF7/plots/MTHFD2_cyt_plot_high_NPL.pdf", device = "pdf", width = 4.5, height = 4.5)

# plot highest values cytosol density
ggplot(data = fuccidata2_levels %>% subset(level_cyt == "high")) + 
  stat_density_2d(aes(x=difTurq, y=dif488_int_btcomp, fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  ggtitle(paste0("MTHFD2 - Cytosol"))+
  xlab("Turquoise2")+
  ylab("Clover")+
  labs(color= "MTHFD2 647")+
  theme_classic()+
  theme(
    plot.title = element_text(color="black", size=14,hjust = 0.5),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  ) +
  scale_x_continuous(trans='log10', limits = c(30,21000)) +
  scale_y_continuous(trans='log10', limits = c(1,21000))
ggsave("MCF7/plots/MTHFD2_cyt_plot_high_NPL_5top_density.pdf", device = "pdf", width = 4.4, height = 3.65)



