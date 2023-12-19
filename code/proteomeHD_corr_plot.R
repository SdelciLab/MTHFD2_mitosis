# libraries
library(tidyverse)
library(pheatmap)
library(rrvgo)
library(ggrepel)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)

# get file 
proteome_GO <-  read.delim("ProteomeHD_GO_bio_process_P13995_thresh0.95.txt", sep = "\t")
colnames(proteome_GO) <- c("Term", "Count", "Bonf", "ID")

# get GO df to plot
proteome_GO_plot <- proteome_GO %>%
  mutate(`-log10(Bonf)` = -log10(Bonf)) %>%
  arrange(desc(`-log10(Bonf)`)) 
duplicated(proteome_GO_plot)
nrow(proteome_GO_plot)

# get only top terms
proteome_GO_plot_red <- proteome_GO_plot[1:23,]

# add function
proteome_GO_plot_red <- proteome_GO_plot_red %>%
  mutate(mitotic = case_when(
    str_detect(Term, c("mitotic")) ~ "mitosis",
    str_detect(Term, c("cell cycle")) ~ "cell cycle",
    str_detect(Term, c("anaphase")) ~ "mitosis",
    T ~ "other")) %>%
  mutate(data = rep("GO Biological Process", nrow(proteome_GO_plot_red)))

# plot only top terms 
ggplot(proteome_GO_plot_red[1:12,], aes(reorder(Term, `-log10(Bonf)`), `-log10(Bonf)`, color = mitotic)) +
  geom_point(aes(size = Count)) +
  labs(x = "", col = "Function") +
  scale_color_manual(values = c("#559A73", "#9EFFC8", "darkgrey")) +
  scale_size(range = c(3,8)) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "lightgrey", linewidth = 0.3) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_text(size = 11), legend.text = element_text(size = 8)) +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 8)) +
  theme(title = element_text(size = 11), plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_line(size = 0.01), panel.grid.major = element_line(size = 0.03))
ggsave("plot_GOBP_proteomeHD_MTHFD2_0.95.pdf", height = 3.5, width = 7.8, device = "pdf")



