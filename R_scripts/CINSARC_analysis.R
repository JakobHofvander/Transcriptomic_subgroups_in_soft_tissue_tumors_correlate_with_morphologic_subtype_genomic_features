# Load required libraries for this script
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

### ssGSEA ####
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(biomaRt)
library(pheatmap)
library(dplyr)
library(tidyverse)
# read in expressions data
df_pc = fread("source_data/tpm_matrix.tsv") %>% column_to_rownames('symbol')
#read in CINSARC gene list
CINSARC = fread("source_data/CINSARC_67_genes.txt")
CINSARC = as.list(CINSARC)
# perform GSVA
GSVA_CINSARC <- gsva(as.matrix(df_pc), CINSARC, method = "gsva",min.sz=10, max.sz=999999, abs.ranking=FALSE, verbose=TRUE)


## save output as a table
GSVA_CINSARC = as.data.frame(t(GSVA_CINSARC))
GSVA_CINSARC = GSVA_CINSARC %>% rownames_to_column()
colnames(GSVA_CINSARC) = c("lab_no", "CINSARC")
write_tsv(GSVA_CINSARC, "GSVA_CINSARC.txt")

## generate tSNE with CINSARC enrichment values
tSNE_fit_25_p30 = fread("source_data/tSNE_fit_25_p30.csv")

p = tSNE_fit_25_p30 %>% 
  inner_join(GSVA_CINSARC) %>% 
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             color = CINSARC))+
  geom_point() + 
  scale_colour_gradient2(low = "darkblue", mid = "grey", high = "darkred") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 12) # Increase y-axis tick label size
  )
ggsave(plot=p, "tSNE_GSVA_CINSARC.pdf", width = 10, height = 7, dpi = 600)


