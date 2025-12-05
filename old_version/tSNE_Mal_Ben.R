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

#### start here

tSNE_fit_25_p30 = fread("source_data/meta_data.txt")

# color vector benign
fread("color_vector_ben", header = T) %>% 
  inner_join(w, by = 'diagnosis') %>% clipr::write_clip()

w %>% distinct(Diagnosisa) %>% view()

color_vector_ben = fread("color_vector_ben", header = T) %>% dplyr::select('9_color')
color_vector_ben = as.vector(color_vector_ben$`9_color`)

name = fread("color_vector_ben", header = T) %>% dplyr::select(Abbreviationm)
names(color_vector_ben) <- name$Abbreviationm

### shape vector
shape_vector_ben = fread("color_vector_ben", header = T) %>% dplyr::select('3_shape')
shape_vector_ben = as.vector(shape_vector_ben$`3_shape`)
names(shape_vector_ben) <- name$Abbreviationm


### color vector malignant ### ### ### ### ### ### ### ### HERE tomorrow
color_vector_mal = fread("color_vector_ben", header = T) %>% dplyr::select('mal_col')
color_vector_mal = as.vector(color_vector_mal$mal_col)

name = fread("color_vector_ben", header = T) %>% dplyr::select(Abbreviationm)
names(color_vector_mal) <- name$Abbreviationm

### shape vector
shape_vector_mal = fread("color_vector_ben", header = T) %>% dplyr::select('mal_shape')
shape_vector_mal = as.vector(shape_vector_mal$mal_shape)
names(shape_vector_mal) <- name$Abbreviationm



#### plot with text
p = tSNE_fit_25_p30 %>%
  #inner_join(shape) %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             shape = mal_shape,
             color = diagnosis))+
  geom_point() +
  scale_color_manual(values=color_vec_65) +
  #scale_shape_manual(values=shape) +
  scale_shape_manual(values = c(15,16,17,3,16)) +
  theme_bw() +
  #geom_point(data=Unclear_25, aes(tSNE1, tSNE2), fill ='black', color = 'black', shape=1) + # 16 = filled black, 1 = black ring no fill
  geom_text_repel(aes(label = tSNE_fit_25_p30$diagnosis),
                  size = 1,
                  #box.padding = unit(0.25, "lines"),
                  #point.padding = unit(0.25, "lines"),
                  min.segment.length = 0,
                  max.overlaps = Inf)
setwd("/Users/jhofvand/Dropbox/Jakob+Figge/Transcriptomics/revised_figures")
ggsave(plot=p, "Mal.704_25_p30_diagnosis.pdf", width = 16, height = 9, dpi=400)

###### No text
p = tSNE_fit_25_p30 %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             color = Abbreviationm,
             shape = Abbreviationm))+
  geom_point() +
  scale_shape_manual(values = c(15,16,4,3,16)) +
  scale_color_manual(values=color_vec_65) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 12) # Increase y-axis tick label size
  )

ggsave(plot=p, "Mal_tsne.704_25_p30_no_lable2.pdf", width = 12, height = 7, dpi = 600)

########## ser bra ut gor samma for benign imorgon!!

#### plot with text
p = tSNE_fit_25_p30 %>%
  #inner_join(shape) %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             shape = ben_shape,
             color = diagnosis))+
  geom_point() +
  scale_color_manual(values=color_vector_ben) +
  #scale_shape_manual(values=shape) +
  scale_shape_manual(values = c(15,16,17,3,16)) +
  theme_bw() +
  #geom_point(data=Unclear_25, aes(tSNE1, tSNE2), fill ='black', color = 'black', shape=1) + # 16 = filled black, 1 = black ring no fill
  geom_text_repel(aes(label = tSNE_fit_25_p30$diagnosis),
                  size = 1,
                  #box.padding = unit(0.25, "lines"),
                  #point.padding = unit(0.25, "lines"),
                  min.segment.length = 0,
                  max.overlaps = Inf)
setwd("/Users/jhofvand/Dropbox/Jakob+Figge/Transcriptomics/revised_figures")
ggsave(plot=p, "Ben.704_25_p30_diagnosis.pdf", width = 16, height = 9, dpi=400)


as.data.frame(color_vector_ben) %>% clipr::write_clip()
###### No text benign 
p = tSNE_fit_25_p30 %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             Diagnosis=Diagnosisa, 
             Class=Class,
             color = Abbreviationm,
             shape = Abbreviationm))+
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

ggsave(plot=p, "Ben_tsne.704_25_p30_no_lable_v2.pdf", width = 10, height = 5.5, dpi = 600)
ggsave(plot=p, "Ben_tsne.704_25_p30_no_lable_v3.pdf", width = 10, height = 5.5, dpi = 600)

# Make it interactive
library(htmlwidgets)
library(plotly)
p_interactive <- ggplotly(p, tooltip = c("color", "Diagnosis", "Class"))
# Save as HTML
saveWidget(p_interactive, "Ben_tsne.704_25_p30_no_lable_v3.html", selfcontained = TRUE)
getwd()




###### No text malignant ###### ###### ###### ###### ###### 
p = tSNE_fit_25_p30 %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             Diagnosis=Diagnosisa, 
             Class=Class,
             color = Abbreviationm,
             shape = Abbreviationm))+
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

ggsave(plot=p, "Mal_tsne.704_25_p30_no_lable_v3.pdf", width = 10, height = 5, dpi = 600)

p_interactive <- ggplotly(p, tooltip = c("color", "Diagnosis", "Class"))
# Save as HTML
saveWidget(p_interactive, "Mal_tsne.704_25_p30_no_lable_v3.html", selfcontained = TRUE)



###################### batch effect ###################### ###################### ###################### 
### color vector run
setwd("/Users/jhofvand/Dropbox/Jakob+Figge/Transcriptomics/source_data_and_code/source_data/")
color_vector_batch = fread("color_vector_batch3", header = T) %>% dplyr::select('year_col')
color_vector_batch = as.vector(color_vector_batch$year_col)

name = fread("color_vector_batch3", header = T) %>% dplyr::select(sequencing_year)
names(color_vector_batch) <- name$sequencing_year

### shape vector
shape_vector_batch = fread("color_vector_batch3", header = T) %>% dplyr::select(year_shape)
shape_vector_batch = as.vector(shape_vector_batch$year_shape)
names(shape_vector_batch) <- name$sequencing_year

tSNE_fit_25_p30 %>% 
  dplyr::mutate(run = strsplit(UID, "_") %>% map_chr(., 2)) %>% clipr::write_clip()

batch = fread("meta_data5.txt")
batch$sequencing_year = as.factor(batch$sequencing_year)
p = #tSNE_fit_25_p30 %>% 
  #dplyr::mutate(run = strsplit(UID, "_") %>% map_chr(., 2)) %>% 
  batch %>% 
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             Diagnosis=Diagnosis,
             Abbreviation=Abbreviation,
             color = sequencing_year,
             shape = sequencing_year))+
  geom_point() +
  #scale_shape_manual(values = c(15,16,17,3,16)) +
  scale_shape_manual(values=shape_vector_batch) +
  scale_color_manual(values=color_vector_batch) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 12) # Increase y-axis tick label size
  )

ggsave(plot=p, "RUN_tsne.704_25_p30_no_lable_year.pdf", width = 12, height = 5.5, dpi = 600)

p_interactive <- ggplotly(p, tooltip = c("color", "Diagnosis", "Abbreviation"))
# Save as HTML
saveWidget(p_interactive, "Batch.704_25_p30_no_lable_v3.html", selfcontained = TRUE)


######




