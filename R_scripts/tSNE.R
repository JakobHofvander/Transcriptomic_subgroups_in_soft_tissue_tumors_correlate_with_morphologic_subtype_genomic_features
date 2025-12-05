# Load required libraries for this script
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(Rtsne)
})

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
s = fread("source_data/tpm_matrix.tsv") %>% column_to_rownames("symbol")
sl = log2(s + 0.01) # log transform TPM values

#filter genes based on variance
mV_25 = varFilter(as.matrix(sl), var.func=IQR, var.cutoff=0.75, filterByQuantile=TRUE) # select top 25 % most variable genes.
# tSNE
tSNE_fit_25_p30 = t(mV_25) %>% Rtsne(perplexity=30, max_iter = 3000) #

#extract tSNE components
tSNE_fit_25_p30 = tSNE_fit_25_p30$Y %>%
  as.data.frame()
colnames(tSNE_fit_25_p30) = c("tSNE1", "tSNE2")
tSNE_fit_25_p30 = tSNE_fit_25_p30 %>% mutate(group=colnames(mV_25))
