### just read in the tSNE co-ordinates

setwd("/Users/jhofvand/Dropbox/Jakob+Figge/Transcriptomics/revised_figures")

# fix the metadata to include all the information I guess ...

tSNE_fit_25_p30 = fread("~/Dropbox/Jakob+Figge/Transcriptomics/source_data_and_code/source_data/meta_data4.txt") %>% 
  dplyr::select(!Cluster)
tSNE_fit_25_p30 = read_excel("~/Dropbox/Jakob+Figge/Transcriptomics/source_data_and_code/table1_revised.xlsx") %>% 
  dplyr::mutate(lab_no = strsplit(Sample, "_") %>% map_chr(., 1)) %>% 
  dplyr::select(lab_no, Diagnosis, Abbreviation, Cluster) %>% 
  inner_join(tSNE_fit_25_p30, by = c('lab_no'))
## should save this a a meta data file for github ...

w = tSNE_fit_25_p30 %>% dplyr::select(Abbreviationm, diagnosis) %>% 
  distinct()

# sheet 3 osv!

fread("source_data/color_vector_ben", header = T)
w = read_excel("~/Dropbox/Jakob+Figge/Transcriptomics/source_data_and_code/table1_revised.xlsx", sheet = "Sheet3")

w = w %>% dplyr::select(Abbreviation, diagnosis) %>% 
  distinct()

w = fread("source_data/color_vector_ben", header = T) %>% 
  inner_join(w, by = c("diagnosis") )

w %>% dplyr::select(!diagnosis) %>% 
  dplyr::select(Abbreviation, '8_color':mal_shape) %>% clipr::write_clip()

setwd("~/Dropbox/Jakob+Figge/Transcriptomics/source_data_and_code/")

#### start here
# read in tSNE coordinates
tSNE_fit_25_p30 = fread("source_data/meta_data.txt")

# create color vector for benign tumors
color_vector_ben = fread("source_data/color_vector_ben", header = T) %>% dplyr::select('9_color')
color_vector_ben = as.vector(color_vector_ben$`9_color`)

name = fread("source_data/color_vector_ben", header = T) %>% dplyr::select(Abbreviation)
names(color_vector_ben) <- name$Abbreviation

### shape vector
shape_vector_ben = fread("source_data/color_vector_ben", header = T) %>% dplyr::select('3_shape')
shape_vector_ben = as.vector(shape_vector_ben$`3_shape`)
names(shape_vector_ben) <- name$Abbreviation


### color vector malignant ### ### ### ### ### ### ### ### HERE tomorrow
color_vector_mal = fread("source_data/color_vector_ben", header = T) %>% dplyr::select('mal_col')
color_vector_mal = as.vector(color_vector_mal$mal_col)

name = fread("source_data/color_vector_ben", header = T) %>% dplyr::select(Abbreviation)
names(color_vector_mal) <- name$Abbreviation

### shape vector
shape_vector_mal = fread("source_data/color_vector_ben", header = T) %>% dplyr::select('mal_shape')
shape_vector_mal = as.vector(shape_vector_mal$mal_shape)
names(shape_vector_mal) <- name$Abbreviation

## tSNE benign and intermediate 
p = tSNE_fit_25_p30 %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             Diagnosis=Diagnosis, 
             Class=Class,
             color = Abbreviation,
             shape = Abbreviation))+
  geom_point() +
  #scale_shape_manual(values = c(15,16,17,3,16)) +
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

# Make interactive version of the plot
library(htmlwidgets)
library(plotly)
p_interactive <- ggplotly(p, tooltip = c("color", "Diagnosis", "Class"))
# Save as HTML
saveWidget(p_interactive, "Ben_tsne.704_25_p30.html", selfcontained = TRUE)


## tSNE malignant
p = tSNE_fit_25_p30 %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             Diagnosis=Diagnosis, 
             Class=Class,
             color = Abbreviation,
             shape = Abbreviation))+
  geom_point() +
  #scale_shape_manual(values = c(15,16,17,3,16)) +
  scale_shape_manual(values=shape_vector_mal) +
  scale_color_manual(values=color_vector_mal) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 12) # Increase y-axis tick label size
  )

ggsave(plot=p, "Mal_tsne.704_25_p30.pdf", width = 10, height = 5, dpi = 600)

# Save as interactive HTML
p_interactive <- ggplotly(p, tooltip = c("color", "Diagnosis", "Class"))
saveWidget(p_interactive, "Mal_tsne.704_25_p30_no_lable_v3.html", selfcontained = TRUE)

