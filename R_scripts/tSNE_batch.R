# Load required libraries for this script
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

#create color vector for sequencing year
setwd("/Users/jhofvand/Dropbox/Jakob+Figge/Transcriptomics/source_data_and_code/")

color_vector_batch = fread("source_data/color_vector_batch", header = T) %>% dplyr::select('year_col')
color_vector_batch = as.vector(color_vector_batch$year_col)

name = fread("source_data/color_vector_batch", header = T) %>% dplyr::select(sequencing_year)
names(color_vector_batch) <- name$sequencing_year

### shape vector
shape_vector_batch = fread("source_data/color_vector_batch", header = T) %>% dplyr::select(year_shape)
shape_vector_batch = as.vector(shape_vector_batch$year_shape)
names(shape_vector_batch) <- name$sequencing_year

batch = fread("source_data/meta_data.txt")
batch$sequencing_year = as.factor(batch$sequencing_year)
p = batch %>% 
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             Diagnosis=Diagnosis,
             Abbreviation=Abbreviation,
             color = sequencing_year,
             shape = sequencing_year))+
  geom_point() +
  scale_shape_manual(values=shape_vector_batch) +
  scale_color_manual(values=color_vector_batch) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 12) # Increase y-axis tick label size
  )

ggsave(plot=p, "batch_tsne.704_25_p30.pdf", width = 12, height = 5.5, dpi = 600)
###




