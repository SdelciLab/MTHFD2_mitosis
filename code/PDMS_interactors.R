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

# list interactor files
interactor_files <- list.files(path = "interactors", pattern = ".txt$", full.names = T)
interactor_names <- c( "chromatin_interactors_Mascot", "cytosol_interactors_Mascot", "cytosol_interactors_MaxQuant", "chromatin_interactors_MaxQuant")

# read interactor files
interactors_read <- map(interactor_files, read_tsv)

# get ensembl ID 
annotation_file <- read_tsv("input_proteins/HUMAN_9606_idmapping.txt", col_names = c("Uniprot", "Database", "Annotation"))
annotation_file_ens <- annotation_file %>%
    filter( Database == "Ensembl") %>%
    dplyr::select(Uniprot, Annotation)
colnames(annotation_file_ens) <- c("Uniprot", "Ensembl")

# join ensembl ID to MS results to have only 1 ensembl IDs / protein
# put ensembl ID first
# remove Ensembl = NA
get_ens <- function(df) {
    df2 <- cbind(df, Ensembl=annotation_file_ens[match(df$Uniprot, annotation_file_ens$Uniprot),"Ensembl"])
    df2 <- df2 %>%
        dplyr::select(Ensembl, everything()) %>%
        filter(!is.na(Ensembl))
}
interactors_ens <- map(interactors_read, get_ens)

# remove logOddsScore column because it is only present in MaxQuant output
rem_logOdds <- function(df) {
    df <- df %>%
        dplyr::select(-logOddsScore)
}
interactors_ens[3:4] <- map(interactors_ens[3:4], rem_logOdds)
map(interactors_ens, colnames)

# check length
map(interactors_read, nrow) # chrom Mascot 2017 cyt Mascot 3400 cyt MaxQuant 3016 chrom MaxQuant 1910
map(interactors_ens, nrow) # chrom Mascot 1996 cyt Mascot 3356 cyt MaxQuant 2985 chrom MaxQuant 1891

# check for duplicated uniprots
map(interactors_ens, ~ table(duplicated(.$Uniprot))) # all unique

# check for duplicated ensembl IDs
map(interactors_ens, ~ table(duplicated(.$Ensembl))) # 1, 1, 1, 2 duplicated genes, respectively 
map(interactors_ens,  ~ .$Ensembl[duplicated(.$Ensembl)]) # 1,2,3: ENSG00000120802 and 4: ENSG00000186184 and ENSG00000120802 -
# these genes give rise to 2 different proteins

# remove duplicated genes
rem_dup_genes <- function(df) {
    df %>% 
        filter(Ensembl != c("ENSG00000186184")) %>%
        filter(Ensembl != c("ENSG00000120802"))
}
interactors_ens <- map(interactors_ens, rem_dup_genes)
map(interactors_ens, nrow) # chrom Mascot 1993 cyt Mascot 3354 cyt MaxQuant 2982 chrom MaxQuant 1887

# get the abundances from previous files
load("input_proteins/MS_results_ens.RData")

interactors_ens_Area_cyt <- interactors_ens[[2]] %>%
    left_join(MS_results_ens[[1]], by = "Ensembl")
interactors_ens_Area_chrom <- interactors_ens[[1]] %>%
    left_join(MS_results_ens[[1]], by = "Ensembl")

interactors_ens_Int_cyt <- interactors_ens[[3]] %>%
    left_join(MS_results_ens[[2]], by = "Ensembl")
interactors_ens_Int_chrom <- interactors_ens[[4]] %>%
    left_join(MS_results_ens[[2]], by = "Ensembl")

interactors_ens_LFQ_cyt <- interactors_ens[[3]] %>%
    left_join(MS_results_ens[[3]], by = "Ensembl")
interactors_ens_LFQ_chrom <- interactors_ens[[4]] %>%
    left_join(MS_results_ens[[3]], by = "Ensembl")

# list and names
interactors_cyt <- list(interactors_ens_Area_cyt, interactors_ens_Int_cyt, interactors_ens_LFQ_cyt)
interactors_cyt_names <- c("cytosol_interactors_Area", "cytosol_interactors_Int", "cytosol_interactors_LFQ")

interactors_chrom <- list(interactors_ens_Area_chrom, interactors_ens_Int_chrom, interactors_ens_LFQ_chrom)
interactors_chrom_names <- c("chromatin_interactors_Area", "chromatin_interactors_Int", "chromatin_interactors_LFQ")

# 2. Heatmap of interactors --------------------------------------------------------------------------------------------------------

# check colnames
map(interactors_cyt, colnames)
map(interactors_chrom, colnames)

# filter values BFDR <= 0.2 and FC >= 5
top_inter_f <- function(df) {
    df <- df %>%
        filter(df$BFDR <= 0.2 & df$FoldChange >= 5) %>%
        dplyr:::select(-c(Bait, Uniprot.x, Uniprot.y, Gene_name, Spec, SpecSum, AvgSpec, NumReplicates, 
                          ctrlCounts, AvgP, MaxP, TopoAvgP, TopoMaxP, SaintScore, FoldChange, BFDR, boosted_by))
}
top_int_cyt <- map(interactors_cyt, top_inter_f)
top_int_chrom <- map(interactors_chrom, top_inter_f)
map(top_int_cyt, nrow) # 127 119 and 119
map(top_int_chrom, nrow) # 43 19 and 19 

# intersect
length(intersect(top_int_cyt[[1]]$Ensembl,top_int_chrom[[1]]$Ensembl)) # 13

# check duplicated rows
map(top_int_cyt, ~ table(duplicated(.$Ensembl))) # F
map(top_int_chrom, ~ table(duplicated(.$Ensembl))) # F

# convert values to matrix
# create matrix 
create_matrix <- function(df) {
    df %>%
        column_to_rownames("Ensembl")
}
interactors_cyt_matrix <- map(top_int_cyt, create_matrix) 
interactors_chrom_matrix <- map(top_int_chrom, create_matrix) 

# change 0 values to NA, log10 of abundances, and then after the log all NA to 0
# change 0 values to NA, log10 of abundances, and then after the log all NA to 0
mod_matrix_1 <- function(m) {
    m[m == 0] <- NA
    m_log <- log10(m)
    m_log[is.na(m_log)] <- 0
    m_log
}
interactors_cyt_matrix_log <- map(interactors_cyt_matrix, mod_matrix_1)
interactors_chrom_matrix_log <- map(interactors_chrom_matrix, mod_matrix_1)

# heatmap with clustering to obtain the order of the rows (genes)
heatmap_order <- function(m_log) {
    set.seed(1234)
    order <- pheatmap(m_log, show_rownames =  F, cluster_rows = T, cluster_cols = F)
    order_genes <- rownames(m_log)[order$tree_row$order]
    order_genes
}
interactors_cyt_heatmap <- map(interactors_cyt_matrix_log, heatmap_order)
interactors_chrom_heatmap <- map(interactors_chrom_matrix_log, heatmap_order)

# reorder the matrix, and change 0 to NAs
mod_matrix_2 <- function(m_log, order_genes) {
    m_log_reorder <-  m_log[order_genes, ]
    m_log_reorder[m_log_reorder == 0] <- NA
    pheatmap(m_log_reorder,  show_rownames =  F, cluster_rows = F, cluster_cols = F)
    m_log_reorder
}
interactors_cyt_heatmap_order <- map2(interactors_cyt_matrix_log, interactors_cyt_heatmap, mod_matrix_2)
interactors_chrom_heatmap_order <- map2(interactors_chrom_matrix_log, interactors_chrom_heatmap, mod_matrix_2)

# add annotation main location
add_sub_loc <- function(df) {
    ens <- unique(df$Ensembl)
    hpa_annot <- getHpa(ens, hpadata = "hpaSubcellularLoc") %>%
        dplyr:::select(Gene, Main.location)
    ens_loc <- hpa_annot %>% 
        separate(Main.location, c("Main.location", "Sub.location"), ";") %>%
        dplyr:::select(Gene, Main.location)
    rownames(ens_loc) <- c()
    ens_loc <- ens_loc %>%
        column_to_rownames("Gene")
    ens_loc
}
interactors_cyt_subloc <- map(interactors_cyt, add_sub_loc)
interactors_chrom_subloc <- map(interactors_chrom, add_sub_loc)
map(interactors_cyt_subloc, ~ table(.$Main.location))
map(interactors_chrom_subloc, ~ table(.$Main.location))

# annotation colors
ann_colors <- list(Main.location = c(`Cell Junctions` = "#003300", `Focal adhesion sites` = "#006633", `Plasma membrane` = "#00CC33", 
                                     `Cytosol` = "#36F736", `Cytoplasmic bodies` = "#66FF99", `Actin filaments` = "#FFFFCC", 
                                     `Microtubules` = "#FFFF99", `Microtubule ends` = "#FFFF00", `Intermediate filaments` = "#CCCC66", 
                                     `Centrosome` = "#FF9933", `Centriolar satellite` = '#E88B00', `Endoplasmic reticulum` = "#CC6633", 
                                     `Golgi apparatus` = "#FF0000", `Mitochondria` = "#990000", `Lysosomes` = "#660000", 
                                     `Endosomes` = "#994C00", `Peroxisomes` = "#8B7F72", `Lipid droplets` = "#947A5F", 
                                     `Rods & Rings` = "#ba9570", `Vesicles` = "#FFCC99", `Nuclear membrane` = "#000066", 
                                     `Nucleoplasm` = "#3535C7", `Nuclear bodies` = "#3366CC", `Nuclear speckles` = "#3399CC", 
                                     `Nucleoli` = "#33CCFF", `Nucleoli fibrillar center` = "#99CCFF", `Nucleoli rim` = "#d2e0fc", 
                                     `Cytokinetic bridge` = "#FFCCFF", `Cleavage furrow` = "#FFCCE5", `Midbody` = "#FF66FF", 
                                     `Midbody ring` = "#FF3399", `Kinetochore` = '#C959EB', `Mitotic spindle` = '#9900CC', 
                                     `Mitotic chromosome` = '#660099'))

# add annotation heatmap
plot_heatmap <- function(m_log_reorder, annot_df, files) {
    pdf(file= paste0("plots/heatmap_", files, ".pdf"), width = 6, height = 7)
    pheatmap(m_log_reorder, 
             show_rownames =  F, cluster_rows = F, cluster_cols = F, 
             color = colorRampPalette(c("#FFFFFF", "#006666"))(1000),
             annotation_row = annot_df, annotation_colors = ann_colors,
             na_col = "whitesmoke", fontsize = 10, fontsize_col = 6)
    dev.off()
}
pwalk(list(interactors_cyt_heatmap_order, interactors_cyt_subloc, interactors_cyt_names), plot_heatmap)
pwalk(list(interactors_chrom_heatmap_order, interactors_chrom_subloc, interactors_chrom_names), plot_heatmap)

# plot separately heatmap area of cytosolic and chromatin separately
pdf(file= "plots/heatmap_area_cyt_int.pdf", width = 3, height = 5)
pheatmap(interactors_cyt_heatmap_order[[1]][,1:6], 
         show_rownames =  F, cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(c("#FFCCCC", "#660000"))(1000),
         annotation_row = interactors_cyt_subloc[[1]], 
         annotation_colors = ann_colors, na_col = "whitesmoke", 
         border_color = NA, fontsize = 6, fontsize_col = 6)
dev.off()

pdf(file= "plots/heatmap_area_chr_int.pdf", width = 3, height = 2.5)
pheatmap(interactors_chrom_heatmap_order[[1]][,7:12], 
         show_rownames =  F, cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(c("#FFCCCC", "#660000"))(1000),
         annotation_row = interactors_chrom_subloc[[1]], 
         annotation_colors = ann_colors, na_col = "whitesmoke", 
         border_color = NA,fontsize = 6, fontsize_col = 6)
dev.off()

# 3. Volcano plot with nuclear localization --------------------------------------------------------------------------------------------------------

# read interactors from Mascot and SAINTS
int_cyt <- read_tsv("interactors/2021LB004-SAINT-cytosol.txt")
int_chrom <- read_tsv("interactors/20201LB004-SAINT-chromatin.txt")

# get ensembl ID 
annotation_file <- read_tsv("input_proteins/HUMAN_9606_idmapping.txt", col_names = c("Uniprot", "Database", "Annotation"))
annotation_file_ens <- annotation_file %>%
    filter(Database == "Ensembl") %>%
    dplyr::select(Uniprot, Annotation)
colnames(annotation_file_ens) <- c("Uniprot", "Ensembl")

# join ensembl ID to MS results to have only 1 ensembl IDs / protein
# put ensembl ID first
# remove Ensembl = NA
int_cyt_ens <- cbind(int_cyt, Ensembl=annotation_file_ens[match(int_cyt$Uniprot, annotation_file_ens$Uniprot),"Ensembl"])
int_cyt_ens <- int_cyt_ens %>%
    dplyr::select(Ensembl, everything()) %>%
    filter(!is.na(Ensembl))
int_chrom_ens <- cbind(int_chrom, Ensembl=annotation_file_ens[match(int_chrom$Uniprot, annotation_file_ens$Uniprot),"Ensembl"])
int_chrom_ens <- int_chrom_ens %>%
    dplyr::select(Ensembl, everything()) %>%
    filter(!is.na(Ensembl))

# check length
nrow(int_cyt_ens) # 3356
nrow(int_chrom_ens) # 1996

# HPA annot add subcellular localisation
ens_cyt <- unique(int_cyt_ens$Ensembl)
ens_chr <- unique(int_chrom_ens$Ensembl)

hpa_annot_cyt <- getHpa(ens_cyt, hpadata = "hpaSubcellularLoc") %>%
    dplyr::select(-Gene.name)
hpa_annot_chr <- getHpa(ens_chr, hpadata = "hpaSubcellularLoc") %>%
    dplyr::select(-Gene.name)

int_cyt_local <- int_cyt_ens %>% left_join(hpa_annot_cyt, by = c("Ensembl" = "Gene"))
int_chr_local <- int_chrom_ens %>% left_join(hpa_annot_chr, by = c("Ensembl" = "Gene"))

int_cyt_local_nucl <- int_cyt_local %>%
    mutate(Nuclear = ifelse(str_detect(Main.location, "Nucl"), "Mainly nuclear", 
                            ifelse(!(str_detect(Main.location, "Nucl")) & str_detect(Additional.location, "Nucl"), 
                                   "Additionally nuclear", "Not nuclear")),
           Nuclear2 = ifelse(is.na(Nuclear), "Not available", Nuclear))
int_chr_local_nucl <- int_chr_local %>%
    mutate(Nuclear = ifelse(str_detect(Main.location, "Nucl"), "Mainly nuclear", 
                            ifelse(!(str_detect(Main.location, "Nucl")) & str_detect(Additional.location, "Nucl"), 
                                   "Additionally nuclear", "Not nuclear")),
           Nuclear2 = ifelse(is.na(Nuclear), "Not available", Nuclear))

table(int_cyt_local_nucl$Nuclear) # Main 1354 (40%), Additional 288 (9%), Not nuclear 1248 + NA (51%)
table(int_cyt_local_nucl$Nuclear2) 
table(int_chr_local_nucl$Nuclear) # Main 1114  (56%), Additional 115  (6%), Not nuclear 504 + NA (38%)
table(int_chr_local_nucl$Nuclear2)

# get all nuclear proteins
int_nuclear_cyt_HPA <- int_cyt_local_nucl %>%
    filter(Nuclear2 %in% c("Mainly nuclear", "Additionally nuclear")) %>%
    dplyr::select(Ensembl, Uniprot, Gene_name, SaintScore, FoldChange, BFDR)
int_nuclear_chr_HPA <- int_chr_local_nucl %>%
    filter(Nuclear2 %in% c("Mainly nuclear", "Additionally nuclear")) %>%
    dplyr::select(Ensembl, Uniprot, Gene_name, SaintScore, FoldChange, BFDR) 
int_nuclear_cht_chr_HPA <- rbind(int_nuclear_cyt_HPA, int_nuclear_chr_HPA)
nrow(int_nuclear_cyt_HPA) # 1642
nrow(int_nuclear_chr_HPA) # 1229
nrow(int_nuclear_cht_chr_HPA) # 2871
length(unique(int_nuclear_cht_chr_HPA$Ensembl)) # 2148
write_csv2(int_nuclear_cht_chr_HPA, "candidates/all_interactors_nuclear_HPA_cyt_chr.csv")

# get df to plot
set.seed(1234)
get_volc_df <- function(df) {
    df <- df %>%
        dplyr::select(Gene_name, Ensembl, Uniprot, BFDR, FoldChange, Nuclear2) %>%
        mutate(BFDR2 = abs(jitter(BFDR, factor = 1))) %>%
        mutate(log2FC = log2(FoldChange), `-log10(BFDR)` = -log10(BFDR2),
               sign = ifelse((log2FC >= 2.3 & BFDR2 <= 0.2), TRUE, FALSE)) 
}
interactors_cyt_volcano <- get_volc_df(int_cyt_local_nucl)
interactors_chrom_volcano <- get_volc_df(int_chr_local_nucl)

# factor Nuclear
interactors_cyt_volcano$Nuclear2 <- factor(interactors_cyt_volcano$Nuclear2, 
                                           levels = c("Mainly nuclear", "Additionally nuclear", "Not nuclear", "Not available"))
interactors_chrom_volcano$Nuclear2 <- factor(interactors_chrom_volcano$Nuclear2, 
                                             levels = c("Mainly nuclear", "Additionally nuclear", "Not nuclear", "Not available"))

# put together cytosolic and chromatin data
interactors_cyt_volcano <- interactors_cyt_volcano %>%
    mutate(location = rep("Cytosol", nrow(interactors_cyt_volcano)))
interactors_chr_volcano <- interactors_chrom_volcano %>%
    mutate(location = rep("Chromatin", nrow(interactors_chrom_volcano)))
interactors_cyt_chr_volcano <- rbind(interactors_cyt_volcano, interactors_chr_volcano)
interactors_cyt_chr_volcano$location <- factor(interactors_cyt_chr_volcano$location, levels = c("Cytosol", "Chromatin"))

# plot volcano plot of cytosolic and chromatin interactors
ggplot(interactors_cyt_chr_volcano, aes(x = log2FC, y = `-log10(BFDR)`)) +
    geom_point(data = interactors_cyt_chr_volcano %>% subset(sign == FALSE), alpha = 0.05, size = 3, stroke = 0, aes(color = Nuclear2)) +
    geom_point(data = interactors_cyt_chr_volcano %>% subset(sign == TRUE), alpha = 1, size = 3, stroke = 0, aes(color = Nuclear2)) +
    facet_wrap(~location, nrow = 2) +
    labs(x = "log2FC(MTHFD2 IP/IgG IP)", y = "-log10(BFDR)", color = "") +
    theme_bw() +
    xlim(c(-3, 7)) +
    scale_color_manual(values=c("Mainly nuclear" = "#3535C7", "Additionally nuclear" = "#9ADCFF", 
                                "Not nuclear" = "#FFFF33", "Not available" = "lightgrey")) +
    theme(legend.position = "right") +
    theme(
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 10), 
        strip.background = element_rect(colour="black", fill="#EDEDED"))
ggsave("plots/volcano_nuclear_HPA_cyt_chr_int.pdf", width = 4, height = 4.5)

# numbers of only significant
interactors_cyt_chr_volcano_sign <- interactors_cyt_chr_volcano %>%
    filter(sign == TRUE)

table(interactors_cyt_chr_volcano_sign$Nuclear2[interactors_cyt_chr_volcano_sign$location == "Cytosol"])
table(interactors_cyt_chr_volcano_sign$Nuclear2[interactors_cyt_chr_volcano_sign$location == "Chromatin"])

# 4. Volcano plot and heatmap of top MTHFD2 nuclear interactors--------------------------------------------------------------------------------------------------------

# put together interactors cyt and chrom - Area
all_interactors_Area <- rbind(interactors_cyt[[1]], interactors_chrom[[1]])
table(duplicated(all_interactors_Area$Gene_name)) # 1212
nrow(all_interactors_Area) # 5347

# keep only nuclear interactors by HPA
# we are going to manually add some interactors: MTHFD2 and PIH1D1, PLK1 and PPP6C, that I know are nuclear by Uniprot
all_nuclear_int <- read_csv2("candidates/all_interactors_nuclear_HPA_cyt_chr.csv")
length(unique(all_nuclear_int$Gene_name)) # 2148
all_nuclear_int_genes <- all_nuclear_int$Gene_name
nuclear_interactors_Area <- all_interactors_Area %>%
    filter(Gene_name %in% all_nuclear_int_genes | 
               Gene_name %in% c("MTHFD2", "PIH1D1", "PLK1", "PPP6C"))
length(unique(nuclear_interactors_Area$Gene_name)) # 2150
setdiff(unique(all_nuclear_int$Gene_name), unique(nuclear_interactors_Area$Gene_name)) # "TMPO" "POLR1D", not relevant

# for those duplicated, I want to keep the highest FC
nuclear_interactors_Area_unique <- nuclear_interactors_Area %>%
    arrange(desc(FoldChange)) %>%
    group_by(Gene_name) %>%
    filter(row_number() == 1) %>%
    ungroup()
nrow(nuclear_interactors_Area_unique) # 2150

# get df for heatmap with PSM
heatmap_top_candidates_PSM <- nuclear_interactors_Area_unique %>%
    arrange(desc(FoldChange)) %>%
    dplyr::select(Gene_name, Spec, ctrlCounts, FoldChange, BFDR) %>%
    separate(Spec, c("MTHFD2_rep1", "MTHFD2_rep2", "MTHFD2_rep3"), "\\|") %>%
    separate(ctrlCounts, c("IgG_rep1", "IgG_rep2", "IgG_rep3"), "\\|") %>%
    filter(FoldChange >= 5 & BFDR <= 0.2 | Gene_name %in% c("KMT5A", "DNMT3B"))
nrow(heatmap_top_candidates_PSM) # 87

# create matrix
heatmap_top_candidates_matrix <- heatmap_top_candidates_PSM[, 1:7] %>%
    column_to_rownames("Gene_name")

# check duplicate names
table(duplicated(heatmap_top_candidates_PSM$Gene_name)) # F

# transform PSM to numeric
heatmap_top_candidates_matrix$MTHFD2_rep1 <- as.numeric(heatmap_top_candidates_matrix$MTHFD2_rep1)
heatmap_top_candidates_matrix$MTHFD2_rep2 <- as.numeric(heatmap_top_candidates_matrix$MTHFD2_rep2)
heatmap_top_candidates_matrix$MTHFD2_rep3 <- as.numeric(heatmap_top_candidates_matrix$MTHFD2_rep3)
heatmap_top_candidates_matrix$IgG_rep1 <- as.numeric(heatmap_top_candidates_matrix$IgG_rep1)
heatmap_top_candidates_matrix$IgG_rep2 <- as.numeric(heatmap_top_candidates_matrix$IgG_rep2)
heatmap_top_candidates_matrix$IgG_rep3 <- as.numeric(heatmap_top_candidates_matrix$IgG_rep3)

# genes with PSM 0 -> NA
heatmap_top_candidates_matrix[heatmap_top_candidates_matrix == 0] <- NA

# read annotation
top_genes_anot <- read_delim("candidates/annot_top_nuclear_volcano.csv", delim = ";") %>%
    column_to_rownames("Gene_name")
annot_colors <- list(Function = c(`cell cycle & mitosis` = "#559A73", `cell proliferation & apoptosis` = "#58e797", 
                                  `chromatin & chromosome organization` = "#8D4FC8", `epigenetics` = "#4877D5", 
                                  `transcription regulation` = "#8AD6E6", `DNA damage & repair` = "#EBA0EC", 
                                  `RNA processing & transport` = "#EDB467", `protein processing` = "#F5F63B", 
                                  `intracellular protein transport` = "#a77038", `metabolism` = "#ef5656", `other` = "#CCCCCC"))

# heatmap with clustering to obtain the order of the rows (genes)
pdf(file= "plots/heatmap_top_int_functions.pdf", width = 4, height = 7)
pheatmap(heatmap_top_candidates_matrix, show_colnames =  F, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("#FFCCCC", "#660000"))(1000),
         na_col = "white",  fontsize = 7, fontsize_row = 6,
         annotation_row = top_genes_anot, annotation_colors = annot_colors)
dev.off()

# add variables for Volcano plot
top_genes_anot_df <- top_genes_anot %>%
    rownames_to_column("Gene_name")

# create Volcano df
set.seed(123)
volcano_top_df <- nuclear_interactors_Area_unique %>%
    mutate(BFDR2 = abs(jitter(BFDR, factor = 1))) %>%
    mutate(log2FC = log2(FoldChange), `-log10(BFDR)` = -log10(BFDR2)) %>%
    left_join(top_genes_anot_df, by = "Gene_name") %>%
    mutate(label = ifelse(Gene_name %in% c("MTHFD2", "TPR", "MAD1L1", "KIF4A", "PRMT1", "PRKAA1", 
                                           "PRKAB1", "PRKAG1", "PRKAG2", "KMT5A", "DNMT3B"), Gene_name, ""))

# plot volcano - top
ggplot(volcano_top_df, aes(x = log2FC, y = `-log10(BFDR)`)) +
    geom_point(data = volcano_top_df %>% subset(is.na(Function)), 
               aes(x = log2FC, y = `-log10(BFDR)`), alpha = 0.2, size = 4, color = "#CCCCCC", stroke=NA) +
    geom_point(data = volcano_top_df %>% subset(!is.na(Function) & label == ""), 
               aes(x = log2FC, y = `-log10(BFDR)`, fill = Function), alpha = 1, shape = 21, size = 4, stroke = 0) +
    geom_point(data = volcano_top_df %>% subset(!is.na(Function) & label != ""), 
               aes(x = log2FC, y = `-log10(BFDR)`, fill = Function), alpha = 1, shape = 21, size = 4, color = "black") +
    geom_text_repel(data = volcano_top_df %>% subset(!is.na(Function)), 
                    aes(color = Function, label = label), size = 2.5, na.rm = TRUE, max.overlaps = 30, fontface = "bold") +
    theme_bw() +
    ggtitle("MTHFD2 top nuclear interactors") +
    labs(x = "log2FC(MTHFD2 IP/IgG IP)", y = "-log10(BFDR)") +
    xlim(c(-3, 7)) +
    scale_color_manual(values=c("cell cycle & mitosis" = "#559A73", "cell proliferation & apoptosis" = "#58e797", 
                                "chromatin & chromosome organization" = "#8D4FC8", "epigenetics" = "#4877D5", 
                                "transcription regulation" = "#8AD6E6", "DNA damage & repair" = "#EBA0EC", 
                                "RNA processing & transport" = "#EDB467", "protein processing" = "#F5F63B", 
                                "intracellular protein transport" = "#a77038", "metabolism" = "#ef5656", "other" = "#666666")) +
    scale_fill_manual(values=c("cell cycle & mitosis" = "#559A73", "cell proliferation & apoptosis" = "#58e797", 
                               "chromatin & chromosome organization" = "#8D4FC8", "epigenetics" = "#4877D5", 
                               "transcription regulation" = "#8AD6E6", "DNA damage & repair" = "#EBA0EC", 
                               "RNA processing & transport" = "#EDB467", "protein processing" = "#F5F63B", 
                               "intracellular protein transport" = "#a77038", "metabolism" = "#ef5656", "other" = "#666666")) +
    theme_classic() + 
    theme(
        axis.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave(paste0("plots/volcano_MTHFD2_top_nuclear.pdf"), device = "pdf", width = 6, height = 3.5)

# pie chart of functions
table_functions <- table(volcano_top_df$Function)
table_functions_order <- table_functions[order(table_functions)]
pie_colors_order <- c("#8D4FC8", "#ef5656", "#EBA0EC", "#a77038", "#58e797", "#F5F63B", "#EDB467",
                               "#666666", "#4877D5", "#8AD6E6", "#559A73")
                               
pdf(file= "plots/function_volcano_pie.pdf", width = 5, height = 5)
pie(table_functions_order, labels = names(table_functions_order), col= pie_colors_order, border = NA)
dev.off()
