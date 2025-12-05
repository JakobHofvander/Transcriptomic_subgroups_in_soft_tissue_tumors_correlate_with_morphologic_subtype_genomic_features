### read in the tSNE co-ordinates
tSNE_fit_25_p30 = fread("source_data/meta_data.txt")

# color vector for benign and intermediate tumors
color_vector_ben = fread("source_data/color_vector", header = T) %>% dplyr::select('9_color')
color_vector_ben = as.vector(color_vector_ben$`9_color`)
name = fread("source_data/color_vector", header = T) %>% dplyr::select(Abbreviation)
names(color_vector_ben) <- name$Abbreviation

# shape vector for benign and intermediate tumors
shape_vector_ben = fread("source_data/color_vector", header = T) %>% dplyr::select('3_shape')
shape_vector_ben = as.vector(shape_vector_ben$`3_shape`)
names(shape_vector_ben) <- name$Abbreviation

# color vector for malignant tumors
color_vector_mal = fread("source_data/color_vector", header = T) %>% dplyr::select('mal_col')
color_vector_mal = as.vector(color_vector_mal$mal_col)
names(color_vector_mal) <- name$Abbreviation

# shape vector for malignant tumors
shape_vector_mal = fread("source_data/color_vector", header = T) %>% dplyr::select('mal_shape')
shape_vector_mal = as.vector(shape_vector_mal$mal_shape)
names(shape_vector_mal) <- name$Abbreviation


## tSNE plot for benign and intermediate tumors
p = tSNE_fit_25_p30 %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             color = diagnosis,
             shape = diagnosis))+
  geom_point() +
  scale_shape_manual(values=shape_vector_ben) +
  scale_color_manual(values=color_vector_ben) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 12) # Increase y-axis tick label size
  )

ggsave(plot=p, "Benign_tsne.704_25_p30.pdf", width = 10, height = 5.5, dpi = 600)

## tSNE plot for malignant tumors
p = tSNE_fit_25_p30 %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             color = diagnosis,
             shape = diagnosis))+
  geom_point() +
  scale_shape_manual(values=shape_vector_mal) +
  scale_color_manual(values=color_vector_mal) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 12) # Increase y-axis tick label size
  )

ggsave(plot=p, "Malignant_tsne.704_25_p30.pdf", width = 10, height = 5, dpi = 600)

