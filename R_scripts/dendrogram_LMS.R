# Load required libraries for this script
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Please note that we can't host data obtained through EGAD00001006628, thus the following script descries how the analysis was performed starting from a TPM count matrix of protein coding genes.
# filter based on variance
mV_25 = varFilter(as.matrix(df_pc), var.func=IQR, var.cutoff=0.75, filterByQuantile=TRUE) # select top 15 % most variable genes. Try a few different

#perform clustering
result <- pvclust(mV_25, method.dist='euclidean', method.hclust="complete", nboot=100) # takes the filtered mV matrix as input. TAKES like 1h
# Convert pvclust result into dendrogram object
dendro_data <- dendro_data(result$hclust, type = "rectangle")

# generate dendrogram
dendro_data$labels %>% 
  inner_join(meta, by = c("label" = "UID")) %>% 
  dplyr::select(label, id)
ggplot() +
  geom_segment(data = segment(dendro_data), 
               aes(x = x, y = y, xend = xend, yend = yend), 
               color = "black", linewidth = 0.25) +
  geom_text(data = dendro_data$labels %>% 
              inner_join(meta2, by = c("label" = "UID")), 
            aes(x = x, y = y, label = `mRNA LMS cluster`), 
            hjust = 0.8, size = 3) +  # Adjust text size and positioning
  theme_minimal() +
  labs(title = "Hierarchical Clustering Dendrogram",
       subtitle = "Clustered data with pvclust",
       x="sample number", y = "Height") + 
  coord_flip()

ggsave("Figure_S5_B.pdf", width = 8, height = 12, bg = "white", dpi = 400)

