# libraries
library(tidyverse)
library(gganatogram)
library(hpar)
library(RColorBrewer)
library(ggsci)
library(pheatmap)
library(biomaRt)
library(VennDiagram)
library(EnhancedVolcano)
library(UpSetR)
library(ggbreak) 

# 1. Get input data --------------------------------------------------------------------------------------------------------

# get input files: areas (Mascot) intensities (MaxQuant) and LFQ (MaxQuant)
MS_results_files <- list.files(path = "input_proteins", pattern = ".csv", full.names = T)

# read files
MS_results_read <- map(MS_results_files, read_delim, delim = ";")

# put together Mascot cyt and chrom
MS_Mascot_Area <- full_join(MS_results_read[[2]], MS_results_read[[1]], by = "Uniprot")
# convert NA to 0
MS_Mascot_Area[is.na(MS_Mascot_Area)] <- 0

# list of files
MS_results <- list(MS_Mascot_Area, MS_results_read[[3]], MS_results_read[[4]])
# names
MS_names <- c("MS_Mascot_Area", "MS_MaxQuant_Int", "MS_MaxQuant_LFQ")

# get ensembl ID 
annotation_file <- read_tsv("input_proteins/HUMAN_9606_idmapping.txt", col_names = c("Uniprot", "Database", "Annotation"))
annotation_file_ens <- annotation_file %>%
    filter(Database == "Ensembl") %>%
    dplyr::select(Uniprot, Annotation)
colnames(annotation_file_ens) <- c("Uniprot", "Ensembl")

# join ensembl ID to MS results to have only 1 Ensembl IDs / protein
# put ensembl ID first
# remove Ensembl = NA
get_ens <- function(df) {
    df2 <- cbind(df, Ensembl=annotation_file_ens[match(df$Uniprot, annotation_file_ens$Uniprot),"Ensembl"])
    df2 <- df2 %>%
        dplyr::select(Ensembl, everything()) %>%
        filter(!is.na(Ensembl))
}
MS_results_ens <- map(MS_results, get_ens)

# check length
map(MS_results, nrow) # Area 4495 Int 3602 LFQ 3602
map(MS_results_ens, nrow) # Area 4433 Int 3559 LFQ 3559

# check for duplicated Uniprots
map(MS_results_ens, ~ table(duplicated(.$Uniprot))) # all unique

#check for duplicated ensembl IDs
map(MS_results_ens, ~ table(duplicated(.$Ensembl))) # 1, 2, 2 duplicated genes, respectively 
map(MS_results_ens,  ~ .$Ensembl[duplicated(.$Ensembl)]) # 1: ENSG00000120802 and 2 and 3: ENSG00000186184 and ENSG00000120802 
# these genes give rise to 2 different proteins

# remove duplicated genes
rem_dup_genes <- function(df) {
    df %>% 
        filter(Ensembl != c("ENSG00000186184")) %>%
        filter(Ensembl != c("ENSG00000120802"))
}
MS_results_ens <- map(MS_results_ens, rem_dup_genes)
map(MS_results_ens, nrow) # Area 4430 Int 3555 LFQ 3555

# abundances files
save(MS_results_ens, file = "input_proteins/MS_results_ens.RData")

# 2. MS quality control ---------------------------------------------------------------------------------------------------

# annotate MS results
subloc_annot <- function(df) {
    ens <- unique(df$Ensembl)
    hpa_annot <- getHpa(ens, hpadata = "hpaSubcellularLoc") %>%
        dplyr:::select(Gene, Main.location)
    df_annot <- df %>%
        left_join(hpa_annot, by = c("Ensembl" = "Gene"))
    df_annot
}
MS_results_annot <- map(MS_results_ens, subloc_annot)
map(MS_results_annot, str)

# write the MS localisation annotation
write_annot <- function(df_annot, files) {
    write_csv(df_annot, paste0("HPA_annotation/HPA_annot_", files, ".csv"))
}
walk2(MS_results_annot, MS_names, write_annot)

# function to edit annotation
# separate more than 1 annotation in main location (1) and sublocation (2) and consider only the main location
# mutate organs to overlap with gganatogram
edit_annot <- function(df) {
    df %>% separate(Main.location, c("Main.location", "Sub.location"), ";") %>%
        mutate(organ = case_when(
            Main.location == "Cytosol" ~ "cytosol",
            Main.location == "Intermediate filaments" ~ "intermediate_filaments",
            Main.location == "Actin filaments" ~ "actin_filaments",
            Main.location == "Focal adhesion sites" ~ "focal_adhesion_sites",
            Main.location == "Centrosome" ~ "centrosome",
            Main.location == "Microtubules" ~ "microtubules",
            Main.location == "Microtubules ends" ~ "microtubules_ends",
            Main.location == "Vesicles" ~ "secreted_proteins",
            Main.location == "Lipid droplets" ~ "lipid_droplets",
            Main.location == "Lysosomes" ~ "lysosomes",
            Main.location == "Peroxisomes" ~ "peroxisomes",
            Main.location == "Endosomes" ~ "endosomes",
            Main.location == "Endoplasmic reticulum" ~ "endoplasmic_reticulum",
            Main.location == "Golgi apparatus" ~ "golgi_apparatus",
            Main.location == "Nuclear bodies" ~ "nuclear_bodies",
            Main.location == "Nuclear membrane" ~ "nuclear_membrane",
            Main.location == "Nucleoplasm" ~ "nucleoplasm",
            Main.location == "Nuclear speckles" ~ "nuclear_speckles",
            Main.location == "Nucleoli" ~ "nucleoli",
            Main.location == "Nucleoli fibrillar center" ~ "nucleoli_fibrillar_center",
            Main.location == "Rods & Rings" ~ "rods_and_rings",
            Main.location == "Mitochondria" ~ "mitochondria",
            Main.location == "Plasma membrane" ~ "plasma_membrane",
            TRUE ~ "NA")) 
}

# edit MS files
MS_results_annot2 <- map(MS_results_annot, edit_annot)

# Get NUMBER and ABUNDANCES of proteins

# calculate means of abundances
calc_mean <- function(df) {
    df_mean <- df %>%
        mutate(cyt_MTHFD2 = (df[[3]] + df[[4]] + df[[5]]) / 3) %>%
        mutate(cyt_IgG = (df[[6]] + df[[7]] + df[[8]]) / 3) %>%
        mutate(chrom_MTHFD2 = (df[[9]] + df[[10]] + df[[11]]) / 3) %>%
        mutate(chrom_IgG = (df[[12]] + df[[13]] + df[[14]]) / 3)
    df_mean <- df_mean[ , c("Ensembl", "Uniprot", "Main.location", "organ", "cyt_MTHFD2", "cyt_IgG", "chrom_MTHFD2", "chrom_IgG")]
    df_mean
}
MS_results_annot_mean <- map(MS_results_annot2, calc_mean)

# modify dataframe 
mod_annot <- function(df) {
    df <- df %>%
        pivot_longer(cols = c(5:8), names_to = "Sample", values_to = "Abundance") 
    df
}
MS_results_annot_mod <- map(MS_results_annot_mean, mod_annot)

# remove proteins with abundance 0
rem_ab0 <- function(df) {
    df <- df %>%
        filter(Abundance != 0)
    df
}
MS_results_annot_mod2 <- map(MS_results_annot_mod, rem_ab0)
map(MS_results_annot_mean, nrow) # 4430, 3555, 3555
map(MS_results_annot_mod, nrow) # 17720, 14220, 14220
map(MS_results_annot_mod2, nrow) # 8945, 9137, 6745

# function to get summary of localisations and abundances and remove organ=NA
organs_summary_abund <- function(df) {
    df <- df %>%
        group_by(Sample, Main.location) %>%
        summarise(n = n(), sum_abund = sum(Abundance)) %>%
        mutate(sum_abund_log = log10(sum_abund)) %>%
        filter(Main.location != "NA")
}
MS_results_organs_abund <- map(MS_results_annot_mod2, organs_summary_abund)
map(MS_results_organs_abund, str)

# calculate total and and total abundance in each sample 
get_total_n_abund <- function(df) {
    df <- df %>%
        group_by(Sample) %>%
        summarise(total_sample_n = sum(n), total_abund = sum(sum_abund)) %>%
        mutate(total_abund_log = log10(total_abund))
    df
}
total_n_abund <- map(MS_results_organs_abund, get_total_n_abund)

# add total n and total abundances
add_total_n_abund <- function(df, total) {
    df <- df %>%
        left_join(total, by = "Sample") %>%
        mutate(rel_n = n/total_sample_n*100)  %>%
        mutate(rel_abund = sum_abund/total_abund*100)
}
MS_results_organs_abund_total <- map2(MS_results_organs_abund, total_n_abund, add_total_n_abund)

# order factor
MS_results_organs_abund_total[[1]]$Sample <- factor(MS_results_organs_abund_total[[1]]$Sample, 
                                                    levels = c("cyt_MTHFD2", "cyt_IgG", "chrom_MTHFD2", "chrom_IgG"), 
                                                    labels = c("Cytosol MTHFD2 IP", "Cytosol IgG IP", "Chromatin MTHFD2 IP", "Chromatin IgG IP"))
MS_results_organs_abund_total[[2]]$Sample <- factor(MS_results_organs_abund_total[[2]]$Sample, 
                                                    levels = c("cyt_MTHFD2", "cyt_IgG", "chrom_MTHFD2", "chrom_IgG"),
                                                    labels = c("Cytosol MTHFD2 IP", "Cytosol IgG IP", "Chromatin MTHFD2 IP", "Chromatin IgG IP"))
MS_results_organs_abund_total[[3]]$Sample <- factor(MS_results_organs_abund_total[[3]]$Sample, 
                                                    levels = c("cyt_MTHFD2", "cyt_IgG", "chrom_MTHFD2", "chrom_IgG"), 
                                                    labels = c("Cytosol MTHFD2 IP", "Cytosol IgG IP", "Chromatin MTHFD2 IP", "Chromatin IgG IP"))

# color factor
MS_results_organs_abund_total[[1]]$Main.location <- factor(MS_results_organs_abund_total[[1]]$Main.location , 
                                                           levels = c("Cell Junctions", "Focal adhesion sites", "Plasma membrane", "Cytosol", 
                                                                      "Cytoplasmic bodies", "Actin filaments", "Microtubules", "Microtubule ends", 
                                                                      "Intermediate filaments", "Centrosome", "Centriolar satellite", "Endoplasmic reticulum", 
                                                                      "Golgi apparatus", "Mitochondria", "Lysosomes", "Endosomes", "Peroxisomes", 
                                                                      "Lipid droplets", "Rods & Rings", "Vesicles", "Nuclear membrane", "Nucleoplasm", 
                                                                      "Nuclear bodies", "Nuclear speckles", "Nucleoli", "Nucleoli fibrillar center", 
                                                                      "Nucleoli rim", "Cytokinetic bridge", "Cleavage furrow", "Midbody", "Midbody ring", 
                                                                      "Kinetochore", "Mitotic spindle", "Mitotic chromosome"))
MS_results_organs_abund_total[[2]]$Main.location <- factor(MS_results_organs_abund_total[[2]]$Main.location , 
                                                           levels = c("Cell Junctions", "Focal adhesion sites", "Plasma membrane", "Cytosol", 
                                                                      "Cytoplasmic bodies", "Actin filaments", "Microtubules", "Microtubule ends", 
                                                                      "Intermediate filaments", "Centrosome", "Centriolar satellite", "Endoplasmic reticulum", 
                                                                      "Golgi apparatus", "Mitochondria", "Lysosomes", "Endosomes", "Peroxisomes", "Lipid droplets", 
                                                                      "Rods & Rings", "Vesicles", "Nuclear membrane", "Nucleoplasm", "Nuclear bodies", "Nuclear speckles", 
                                                                      "Nucleoli", "Nucleoli fibrillar center", "Nucleoli rim", "Cytokinetic bridge", "Cleavage furrow", 
                                                                      "Midbody", "Midbody ring", "Kinetochore", "Mitotic spindle", "Mitotic chromosome"))
MS_results_organs_abund_total[[3]]$Main.location <- factor(MS_results_organs_abund_total[[3]]$Main.location , 
                                                           levels = c("Cell Junctions", "Focal adhesion sites", "Plasma membrane", "Cytosol", 
                                                                      "Cytoplasmic bodies", "Actin filaments", "Microtubules", "Microtubule ends", 
                                                                      "Intermediate filaments", "Centrosome", "Centriolar satellite", "Endoplasmic reticulum", 
                                                                      "Golgi apparatus", "Mitochondria", "Lysosomes", "Endosomes", "Peroxisomes", "Lipid droplets", 
                                                                      "Rods & Rings", "Vesicles", "Nuclear membrane", "Nucleoplasm", "Nuclear bodies", "Nuclear speckles", 
                                                                      "Nucleoli", "Nucleoli fibrillar center", "Nucleoli rim", "Cytokinetic bridge", "Cleavage furrow", 
                                                                      "Midbody", "Midbody ring", "Kinetochore", "Mitotic spindle", "Mitotic chromosome"))

# Define the colors you want
ann_colors <- c(`Cell Junctions` = "#003300", `Focal adhesion sites` = "#006633", `Plasma membrane` = "#00CC33", `Cytosol` = "#36F736", 
                `Cytoplasmic bodies` = "#66FF99", `Actin filaments` = "#FFFFCC", `Microtubules` = "#FFFF99", `Microtubule ends` = "#FFFF00", 
                `Intermediate filaments` = "#CCCC66", `Centrosome` = "#FF9933", `Centriolar satellite` = '#E88B00', `Endoplasmic reticulum` = "#CC6633", 
                `Golgi apparatus` = "#FF0000", `Mitochondria` = "#990000", `Lysosomes` = "#660000", `Endosomes` = "#994C00", `Peroxisomes` = "#8B7F72", 
                `Lipid droplets` = "#947A5F", `Rods & Rings` = "#ba9570", `Vesicles` = "#FFCC99", `Nuclear membrane` = "#000066", `Nucleoplasm` = "#3535C7", 
                `Nuclear bodies` = "#3366CC", `Nuclear speckles` = "#3399CC", `Nucleoli` = "#33CCFF", `Nucleoli fibrillar center` = "#99CCFF", 
                `Nucleoli rim` = "#d2e0fc", `Cytokinetic bridge` = "#FFCCFF", `Cleavage furrow` = "#FFCCE5", `Midbody` = "#FF66FF", `Midbody ring` = "#FF3399", 
                `Kinetochore` = '#C959EB', `Mitotic spindle` = '#9900CC', `Mitotic chromosome` = '#660099')

# barplot titles
barplot_titles <- c("MTHFD2 PD-MS Localisation (Mascot, Area)", "MTHFD2 PD-MS Localisation (MaxQuant, Int)", "MTHFD2 PD-MS Localisation (MaxQuant, LFQ)")

# bar plot n
plot_barplot_n <- function(df, title, files) {
    ggplot(df, aes(x = Main.location, y = rel_n, fill = Main.location)) +
        geom_col(alpha = 0.7) +
        facet_wrap(~Sample, ncol = 4) +
        coord_flip() +
        scale_fill_manual(values = ann_colors) +
        theme_bw() +
        theme(legend.position = "none", strip.background = element_rect(colour="black", fill="#EDEDED"),  panel.grid.minor = element_line(size = 0.01),
              panel.grid.major = element_line(size = 0.03)) +
        ggtitle(title) +
        labs(y = "\nN of proteins in each compartment / Total N of proteins in sample (%)", x = "Main Location")
    ggsave(paste0("plots/PD-MS_localisation_n_", files, ".pdf"), device = "pdf", height = 5, width = 8)
} 
pwalk(list(MS_results_organs_abund_total, barplot_titles, MS_names), plot_barplot_n)

# bar plot abundance
plot_barplot_abund <- function(df, title, files) {
    ggplot(df, aes(x = Main.location, y = rel_abund, fill = Main.location)) +
        geom_col(alpha = 0.7) +
        facet_wrap(~Sample, ncol = 4) +
        coord_flip() +
        scale_fill_manual(values = ann_colors) +
        theme_bw() +
        theme(legend.position = "none", strip.background = element_rect(colour="black", fill="#EDEDED"), panel.grid.minor = element_blank(),
              panel.grid.major = element_blank()) +
        ggtitle(title) +
        labs(y = "\nAbundance of proteins in each compartment / Total abundance of proteins in sample (%)", x = "Main Location")
    ggsave(paste0("plots/PD-MS_localisation_abund_", files, ".pdf"), device = "pdf", height = 5, width = 8)
} 
pwalk(list(MS_results_organs_abund_total, barplot_titles, MS_names), plot_barplot_abund)


