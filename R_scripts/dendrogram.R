library(ggrepel)
library(data.table)
library(tidyverse)
library(pheatmap)
library(stringr)
library(genefilter)
library(Rtsne)
library(ggrepel)
library(readxl)
library(ggdendro)
library(pvclust)

# read in tpm matrix
s = fread("source_data/tpm_matrix.tsv") %>% column_to_rownames("symbol")
sl = log2(s + 0.01) # log transform TPM values

#filter genes based on variance
mV_25 = varFilter(as.matrix(sl), var.func=IQR, var.cutoff=0.75, filterByQuantile=TRUE) # select top 25 % most variable genes.

# add annotations
meta = fread("source_data/meta_data.txt") %>% dplyr::mutate(lab_dig = paste(lab_no, Abbreviation, sep = '_'))

mV_25_2 = as.data.frame(mV_25) %>% rownames_to_column() %>% 
  pivot_longer(!rowname) %>% 
  inner_join(meta, by = c("name" = "lab_no")) %>% 
  dplyr::select(rowname, lab_dig, value) %>% 
  pivot_wider(names_from = lab_dig, values_from = value) %>% 
  column_to_rownames("rowname")

# cluster data
result <- pvclust(mV_25_2, method.dist='euclidean', method.hclust="complete", nboot=100) # takes the filtered mV matrix as input. TAKES like 1h

#extract positions
dendro_data <- dendro_data(result$hclust, type = "rectangle")

# Plot dendrogram using ggplot2
ggplot() +
  geom_segment(data = segment(dendro_data), 
               aes(x = x, y = y, xend = xend, yend = yend), 
               color = "black", linewidth = 0.25) +
  #geom_text(data = label(dendro_data), 
  #          aes(x = x, y = y, label = label), 
  #          hjust = 0.8, size = 1) +  # Adjust text size and positioning
  theme_minimal() +
  labs(title = "Hierarchical Clustering Dendrogram",
       subtitle = "Clustered data with pvclust",
       x="sample number", y = "Height") + 
  coord_flip()

ggsave("dendrogram.pdf", width = 12, height = 41, bg = "white", dpi = 400)


#plot annotations
meta = fread("meta_data.txt")
p = dendro_data$labels %>% 
  dplyr::mutate(lab_no = strsplit(label, "_") %>% map_chr(., 1)) %>% 
  inner_join(meta) %>% 
  mutate(value = 1)  %>%
  mutate(across(c(label), ~str_replace_all(., "-", "_"))) %>% 
  mutate(across(c(label), ~str_replace_all(., " ", "_"))) %>% 
  ggplot(aes(x = x,
             y = 1,
             color = diagnosis,
             shape = shape))+
  geom_point() + 
  scale_shape_manual(values = c(15,16,17,3,16)) +
  scale_color_manual(values=color_vec_65) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 12) # Increase y-axis tick label size
  ) +  geom_text(aes(x = x, y = 1, label = label), 
            hjust = 0.8, size = 1) + coord_flip()  # Adjust text size and positioning

ggsave(plot=p, "annotation_dendrogram.pdf", width = 12, height = 47, device = "pdf")


