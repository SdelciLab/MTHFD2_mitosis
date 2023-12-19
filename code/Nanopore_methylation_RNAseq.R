# libraries
library(tidyverse)
library(nortest)
library(biomaRt)

# Objective: cross cpg methylation with RNA-seq

# 1. Get input data -------------------------------------------------------------------------------
library(tidyverse)
library(nortest)
library(biomaRt)

# read RNA-seq
KO1_t07_common <- read.csv("DEG_KO1_vs_WT_t0t7_common_annotated_shrink.csv")

# n
nrow(KO1_t07_common) # 1192

# read methylation data
hypomethyl_data <-  read.delim("neg.anno", sep = "\t")
hypermethyl_data <-  read.delim("pos.anno", sep = "\t")

# get location annotation
hypomethyl_data_annot <- hypomethyl_data %>%
    mutate(annot_int = str_remove(Annotation, " .*"))
table(hypomethyl_data_annot$annot_int)

hypermethyl_data_annot <- hypermethyl_data %>%
    mutate(annot_int = str_remove(Annotation, " .*"))
table(hypermethyl_data_annot$annot_int)

# divide gene up and down
KO1_t07_common_up <- KO1_t07_common[KO1_t07_common$log2FoldChange_t0 > 0 
                                    & KO1_t07_common$log2FoldChange_t7 > 0,]
nrow(KO1_t07_common_up) #  469

KO1_t07_common_down <- KO1_t07_common[KO1_t07_common$log2FoldChange_t0 < 0
                                      & KO1_t07_common$log2FoldChange_t7 < 0,]
nrow(KO1_t07_common_down) #  723

# get methylation info of promoter, exon, intron and TTS
hypomethyl_data_annot_prom <- hypomethyl_data_annot %>%
    filter(annot_int == "promoter-TSS")
hypomethyl_data_annot_exon <- hypomethyl_data_annot %>%
    filter(annot_int == "exon")
hypomethyl_data_annot_intron <- hypomethyl_data_annot %>%
    filter(annot_int == "intron")
hypomethyl_data_annot_TTS <- hypomethyl_data_annot %>%
    filter(annot_int == "TTS")

hypermethyl_data_annot_prom <- hypermethyl_data_annot %>%
    filter(annot_int == "promoter-TSS")
hypermethyl_data_annot_exon <- hypermethyl_data_annot %>%
    filter(annot_int == "exon")
hypermethyl_data_annot_intron <- hypermethyl_data_annot %>%
    filter(annot_int == "intron")
hypermethyl_data_annot_TTS <- hypermethyl_data_annot %>%
    filter(annot_int == "TTS")

# analysis % genes differentially methylated from the total differentially expressed
KO1_t07_common_up_methyl <- data.frame(DE_total = nrow(KO1_t07_common_up),
                                       hypo_prom = length(intersect(KO1_t07_common_up$HGNC_symbol, hypomethyl_data_annot_prom$Entrez.ID)),
                                       hypo_exon = length(intersect(KO1_t07_common_up$HGNC_symbol, hypomethyl_data_annot_exon$Entrez.ID)),
                                       hypo_intron = length(intersect(KO1_t07_common_up$HGNC_symbol, hypomethyl_data_annot_intron$Entrez.ID)),
                                       hypo_TTS = length(intersect(KO1_t07_common_up$HGNC_symbol, hypomethyl_data_annot_TTS$Entrez.ID)),
                                       hyper_prom = length(intersect(KO1_t07_common_up$HGNC_symbol, hypermethyl_data_annot_prom$Entrez.ID)),
                                       hyper_exon = length(intersect(KO1_t07_common_up$HGNC_symbol, hypermethyl_data_annot_exon$Entrez.ID)),
                                       hyper_intron = length(intersect(KO1_t07_common_up$HGNC_symbol, hypermethyl_data_annot_intron$Entrez.ID)),
                                       hyper_TTS = length(intersect(KO1_t07_common_up$HGNC_symbol, hypermethyl_data_annot_TTS$Entrez.ID)))

KO1_t07_common_down_methyl <- data.frame(DE_total = nrow(KO1_t07_common_down),
                                         hypo_prom = length(intersect(KO1_t07_common_down$HGNC_symbol, hypomethyl_data_annot_prom$Entrez.ID)),
                                         hypo_exon = length(intersect(KO1_t07_common_down$HGNC_symbol, hypomethyl_data_annot_exon$Entrez.ID)),
                                         hypo_intron = length(intersect(KO1_t07_common_down$HGNC_symbol, hypomethyl_data_annot_intron$Entrez.ID)),
                                         hypo_TTS = length(intersect(KO1_t07_common_down$HGNC_symbol, hypomethyl_data_annot_TTS$Entrez.ID)),
                                         hyper_prom = length(intersect(KO1_t07_common_down$HGNC_symbol, hypermethyl_data_annot_prom$Entrez.ID)),
                                         hyper_exon = length(intersect(KO1_t07_common_down$HGNC_symbol, hypermethyl_data_annot_exon$Entrez.ID)),
                                         hyper_intron = length(intersect(KO1_t07_common_down$HGNC_symbol, hypermethyl_data_annot_intron$Entrez.ID)),
                                         hyper_TTS = length(intersect(KO1_t07_common_down$HGNC_symbol, hypermethyl_data_annot_TTS$Entrez.ID)))


# add DE
KO1_t07_common_up_methyl <- KO1_t07_common_up_methyl %>%
    mutate(DE = rep("up", nrow(KO1_t07_common_up_methyl)))

KO1_t07_common_down_methyl <- KO1_t07_common_down_methyl %>%
    mutate(DE = rep("down", nrow(KO1_t07_common_down_methyl)))

# put together
table(colnames(KO1_t07_common_up_methyl) == colnames(KO1_t07_common_down_methyl))
KO1_t07_common_methyl_data <- rbind(KO1_t07_common_up_methyl, KO1_t07_common_down_methyl)

# modify data frame
KO1_t07_common_methyl_long <- KO1_t07_common_methyl_data %>%
    pivot_longer(-c(DE_total,DE), names_to = "variables", values_to = "values") %>%
    mutate(methylation = case_when(
        str_detect(variables, "hypo") ~ "hypo",
        str_detect(variables, "hyper") ~ "hyper"),
        region = case_when(
            str_detect(variables, "prom") ~ "promoter-TSS",
            str_detect(variables, "exon") ~ "exon",
            str_detect(variables, "intron") ~ "intron",
            str_detect(variables, "TTS") ~ "TTS",
        ), perc_DE = values/DE_total*100)

# factor
KO1_t07_common_methyl_long$DE <- factor(KO1_t07_common_methyl_long$DE, levels = c("up", "down"))
KO1_t07_common_methyl_long$methylation <- factor(KO1_t07_common_methyl_long$methylation, levels = c("hyper", "hypo"))
KO1_t07_common_methyl_long$region <- factor(KO1_t07_common_methyl_long$region, 
                                            levels = c("promoter-TSS", "exon", "intron", "TTS"))


# plot
ggplot(KO1_t07_common_methyl_long, aes(x=region, y = perc_DE, fill = methylation)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~DE) +
    ylab("Percentage (%)") +
    theme_bw() +
    coord_flip() +
    scale_fill_manual(values = c("hyper" = "#EF3D3A", "hypo" = "#33669B")) +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 11), 
          strip.text = element_text(size = 11), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots_all/KO1_t07_common_methylation.pdf", device = "pdf", width = 6, height = 2.5)

# intersect DE genes with methylation data
hypo_methyl_intersect_t07 <- data.frame(prom_total = length(unique(hypomethyl_data_annot_prom$Entrez.ID)),
                                        exon_total = length(unique(hypomethyl_data_annot_exon$Entrez.ID)),
                                        intron_total = length(unique(hypomethyl_data_annot_intron$Entrez.ID)),
                                        TTS_total = length(unique(hypomethyl_data_annot_TTS$Entrez.ID)),
                                        prom_up = length(intersect(hypomethyl_data_annot_prom$Entrez.ID, KO1_t07_common_up$HGNC_symbol)),
                                        exon_up = length(intersect(hypomethyl_data_annot_exon$Entrez.ID, KO1_t07_common_up$HGNC_symbol)),
                                        intron_up = length(intersect(hypomethyl_data_annot_intron$Entrez.ID, KO1_t07_common_up$HGNC_symbol)),
                                        TTS_up = length(intersect(hypomethyl_data_annot_TTS$Entrez.ID, KO1_t07_common_up$HGNC_symbol)),
                                        prom_down = length(intersect(hypomethyl_data_annot_prom$Entrez.ID, KO1_t07_common_down$HGNC_symbol)),
                                        exon_down = length(intersect(hypomethyl_data_annot_exon$Entrez.ID, KO1_t07_common_down$HGNC_symbol)),
                                        intron_down = length(intersect(hypomethyl_data_annot_intron$Entrez.ID, KO1_t07_common_down$HGNC_symbol)),
                                        TTS_down = length(intersect(hypomethyl_data_annot_TTS$Entrez.ID, KO1_t07_common_down$HGNC_symbol)))

hyper_methyl_intersect_t07 <- data.frame(prom_total = length(unique(hypermethyl_data_annot_prom$Entrez.ID)),
                                         exon_total = length(unique(hypermethyl_data_annot_exon$Entrez.ID)),
                                         intron_total = length(unique(hypermethyl_data_annot_intron$Entrez.ID)),
                                         TTS_total = length(unique(hypermethyl_data_annot_TTS$Entrez.ID)),
                                         prom_up = length(intersect(hypermethyl_data_annot_prom$Entrez.ID, KO1_t07_common_up$HGNC_symbol)),
                                         exon_up = length(intersect(hypermethyl_data_annot_exon$Entrez.ID, KO1_t07_common_up$HGNC_symbol)),
                                         intron_up = length(intersect(hypermethyl_data_annot_intron$Entrez.ID, KO1_t07_common_up$HGNC_symbol)),
                                         TTS_up = length(intersect(hypermethyl_data_annot_TTS$Entrez.ID, KO1_t07_common_up$HGNC_symbol)),
                                         prom_down = length(intersect(hypermethyl_data_annot_prom$Entrez.ID, KO1_t07_common_down$HGNC_symbol)),
                                         exon_down = length(intersect(hypermethyl_data_annot_exon$Entrez.ID, KO1_t07_common_down$HGNC_symbol)),
                                         intron_down = length(intersect(hypermethyl_data_annot_intron$Entrez.ID, KO1_t07_common_down$HGNC_symbol)),
                                         TTS_down = length(intersect(hypermethyl_data_annot_TTS$Entrez.ID, KO1_t07_common_down$HGNC_symbol)))

# add methyl
hypo_methyl_intersect_t07 <- hypo_methyl_intersect_t07 %>%
    mutate(methylation = rep("hypo", nrow(hypo_methyl_intersect_t07)))

hyper_methyl_intersect_t07 <- hyper_methyl_intersect_t07 %>%
    mutate(methylation = rep("hyper", nrow(hyper_methyl_intersect_t07)))

# put together
table(colnames(hypo_methyl_intersect_t07) == colnames(hyper_methyl_intersect_t07))
hyper_hypo_KO1_t07_data <- rbind(hypo_methyl_intersect_t07, hyper_methyl_intersect_t07)

# modify data frame
hyper_hypo_KO1_t07_long <- hyper_hypo_KO1_t07_data %>%
    pivot_longer(-c(methylation, prom_total, exon_total,
                    intron_total, TTS_total),
                 names_to = "variables", values_to = "values") %>%
    mutate(DE = case_when(
        str_detect(variables, "up") ~ "up",
        str_detect(variables, "down") ~ "down"),
        region = case_when(
            str_detect(variables, "prom") ~ "promoter-TSS",
            str_detect(variables, "exon") ~ "exon",
            str_detect(variables, "intron") ~ "intron",
            str_detect(variables, "TTS") ~ "TTS"),
        perc = case_when(
            region == "promoter-TSS" ~ values/prom_total*100,
            region == "exon" ~ values/exon_total*100,
            region == "intron" ~ values/intron_total*100,
            region == "TTS" ~ values/TTS_total*100))

# factor
hyper_hypo_KO1_t07_long$DE <- factor(hyper_hypo_KO1_t07_long$DE, levels = c("up", "down"))
hyper_hypo_KO1_t07_long$methylation <- factor(hyper_hypo_KO1_t07_long$methylation, levels = c("hyper", "hypo"))
hyper_hypo_KO1_t07_long$region <- factor(hyper_hypo_KO1_t07_long$region, 
                                         levels = c("promoter-TSS", "exon", "intron", "TTS"))


# plot
ggplot(hyper_hypo_KO1_t07_long, aes(x=region, y = perc, fill = DE)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~methylation) +
    ylab("Percentage (%)") +
    theme_bw() +
    coord_flip() +
    scale_fill_manual(values = c("down" = "#aaafb0", "up" = "#e56e58")) +
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 11),
          strip.text = element_text(size = 11), 
          strip.background = element_rect(colour="black", fill="#EDEDED"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("plots_all/KO1_t07_methylation_expression.pdf", device = "pdf", width = 6, height = 2.5)
