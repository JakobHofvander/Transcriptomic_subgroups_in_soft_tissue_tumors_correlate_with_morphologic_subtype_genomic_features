# Load required libraries for this script
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(ggpubr)
})

# read in sil score generated with python script
sil_data = fread("source_data/Silhouette_per_diagnosis.tsv")
  
# add meta data
sil_data = fread("source_data/meta_data.txt") %>%
  dplyr::select(lab_no, Driver, Complexity) %>% 
  inner_join(sil_data, by = "lab_no")

# barplot per diagnosis
p = sil_data %>%
  group_by(Abbreviation) %>%
  mutate(number = n()) %>%
  filter(number > 2) %>%
  mutate(median = median(Silhouette_Score)) %>%
  ggplot(aes(reorder(Abbreviation, median) , Silhouette_Score, fill = Driver)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = c("#A7A9AC","#6DC4A6", "#CA3F00", "#0057FF")) +
  coord_flip() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 12) # Increase y-axis tick label size
  )
ggsave(p, file="Boxplot_sil_score.pdf", width = 6, height = 10, dpi = 600)
#boxplot complexity
sil_data %>%
  filter(Complexity != "NA") %>%
  ggplot(aes(Complexity, Silhouette_Score)) +
  geom_boxplot() + theme_bw() +
  stat_compare_means(method = "t.test", label.x = 1.3, size = 4) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_rect(fill = "white"),  # Set background to white
    plot.background = element_rect(fill = "white")   # Set the entire plot background to white
  )
ggsave("Boxplot_Silscore_Complexity.pdf", height = 4.5, width = 5, dpi = 300)
