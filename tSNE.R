library(ggrepel)
library(data.table)
library(tidyverse)
library(pheatmap)
library(stringr)
library(genefilter)
library(Rtsne)
library(ggrepel)
library(readxl)

#### tpm
s = fread("tpm_matrix.tsv") %>% column_to_rownames("symbol")
sl = log2(s + 0.01) # log transform TPM values

#filter genes based on variance
mV_25 = varFilter(as.matrix(sl), var.func=IQR, var.cutoff=0.75, filterByQuantile=TRUE) # select top 25 % most variable genes.
# tSNE
tSNE_fit_25_p30 = t(mV_25) %>% Rtsne(perplexity=30, max_iter = 3000) #

#extract tSNE components
tSNE_fit_25_p30 = tSNE_fit_25_p30$Y %>%
  as.data.frame()
colnames(tSNE_fit_25_p30) = c("tSNE1", "tSNE2")
tSNE_fit_25_p30 = tSNE_fit_25_p30 %>% mutate(group=colnames(mV_25)) # adds sample_name as rowname
#write_csv(tSNE_fit_25_p30, "tSNE_fit_25_p30.csv")

## tSNE plot
# data
tSNE_fit_25_p30 = fread("tSNE_fit_25_p30.csv")

# color
hex = fread("color_vector.txt")
color_vec_65 = hex$hex
names(color_vec_65) <- hex$diagnosis

# shape
shape = fread("shape.csv")

# plot
p = inner_join(tSNE_fit_25_p30, shape, by = c("lab_no")) %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
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
  )

ggsave(plot=p, "tsne.704_25_p30_no_lable.pdf", width = 12, height = 7, dpi = 600)
