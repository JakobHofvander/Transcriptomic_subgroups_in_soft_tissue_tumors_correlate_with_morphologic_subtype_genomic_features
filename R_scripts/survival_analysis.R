# Load required libraries for this script
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(survival)
})

#### 
library(survminer)
library(survival)
library(data.table)
library(tidyverse)
library(ggplot2)
# read in cinsarc enrichment scores
cin = fread("source_data/GSVA_CINSARC.txt") %>% 
  dplyr::mutate(lab_no = strsplit(UID, "_") %>% map_chr(., 1)) %>% 
  dplyr::select(lab_no, CINSARC)

# read in meta_data
ano = fread("source_data/meta_data.txt")

# create intervals
ano = inner_join(ano, cin) %>% 
  mutate(intervals = case_when(
    CINSARC < -0.3 ~ "Low",
    CINSARC < 0.3 ~ "Med",
    CINSARC < 1 ~ "High"
    ))
# boxplot
ano %>% filter(Complexity != 'NA') %>% 
  mutate(new_class = paste(Class, Complexity, sep = '_')) %>% 
  ggplot(aes(CINSARC, new_class)) +
  geom_boxplot() +
  theme_bw()
ggsave("Figure_3B.pdf", dpi = 400, width = 6, height = 5)

#remove benign tumors from further comparisons
ano = ano %>% filter(Class != "Benign")
ano$MFS_months_10 = as.numeric(ano$MFS_months_10)
ano$Met_10 = as.numeric(ano$Met_10)

# KM curves for the cinsarc intervals 
fit <- survfit(Surv(MFS_months_10, Met_10) ~ intervals, data = ano)
ggsurvplot(fit,
           pval = TRUE,
           conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red", "blue", "#808080"),
           xlab = "Time in months",
           ylab = "Met free probability",
           font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(18),
           pval.coord = c(80, 0.92),
           pval.size = 8,
           xlim = c(0, 120))
ggsave("Figure_3C.pdf", dpi = 400, width = 6, height = 5)

##### KM analysis of SC clusters
SC = ano %>% 
  filter(grepl("Complex", subcluster)) 

SC$MFS_months_10 = as.numeric(SC$MFS_months_10)
SC$Met_10 = as.numeric(SC$Met_10)

#
fit <- survfit(Surv(MFS_months_10, Met_10) ~ subcluster, data = SC)
ggsurvplot(fit,
           pval = TRUE,
           #pval.method = TRUE,
           conf.int = FALSE,
           risk.table = FALSE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red", "blue", "#808080"),
           xlab = "Time in months",
           ylab = "Met free probability",
           font.main = c(16, "bold", "darkblue"),
           font.x = c(18),
           font.y = c(18),
           font.tickslab = c(18),
           pval.coord = c(80, 0.92),
           pval.size = 8,
           xlim = c(0, 120))
ggsave("Figure_3E.pdf", dpi = 400, width = 6, height = 5)

## boxplot SC vs CINSARC
SC %>% 
  ggplot(aes(subcluster, CINSARC)) +
  geom_boxplot(outlier.colour="NA") +
  geom_jitter(size = 1, width = 0.2) +
  #stat_summary(fun = mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, color="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.background = element_rect(fill = "white"),  # Set background to white
        plot.background = element_rect(fill = "white"))  + # Set the entire plot background to white)
  ggtitle("") +
  xlab("") +
  ylab("Enrichment score")

ggsave("Figure_3F.pdf", dpi = 400, width = 6, height = 5)
## forest plots
#drop na
ano_na = ano %>% 
  filter(Class != 'Benign') %>% 
  filter(!grepl("NA", Complexity)) %>% 
  filter(!grepl("NA", Met_10))


## cox model for CINSARC
model = coxph(Surv(MFS_months_10, Met_10) ~ CINSARC + Complexity, data =  ano_na)

# cox model for subclusters
SC_model = coxph(Surv(MFS_months_10, Met_10) ~ CINSARC + subcluster, data =  SC)

##### forest plots
ggforest(model)
ggsave("Figure_S7_A.pdf", dpi = 400, width = 6, height = 5)

##### forest plots CINSARC
ggforest(SC_model)
ggsave("Figure_S7_B.pdf", dpi = 400, width = 6, height = 5) # 


