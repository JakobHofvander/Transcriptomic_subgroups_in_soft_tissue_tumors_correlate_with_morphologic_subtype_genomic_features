# Load required libraries
library(tidyverse)
library(readxl)
library(ggpubr)

# Load TPM matrix and convert to data frame
df_pc <- read_tsv("source_data/tpm_matrix.tsv")

# select samples and genes
markers = c("MDM2")

# load metadata
ano <- read_tsv("meta_data.txt") 

plot_df = df_pc %>% filter(symbol %in% markers) %>% 
  pivot_longer(!symbol) %>% 
  mutate(log2TPM = log2(value + 0.01)) %>% 
  inner_join(ano, by = c("name" = "lab_no")) %>% 
  filter(Abbreviation == "ALT" | diagnosis == "Hibernoma" | diagnosis == "Lipoma" | diagnosis == "Lipoma spindle" | diagnosis == "Lipoblastoma")

# MDM2 boxplot
plot_df %>% 
  filter(symbol == "MDM2") %>% 
  ggplot(aes(diagnosis, log2TPM), fill = Cluster) +
  geom_boxplot(outlier.colour="NA") +
  geom_jitter(size = 1, width = 0.2) +
  #stat_summary(fun = mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, color="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.background = element_rect(fill = "white"),  # Set background to white
        plot.background = element_rect(fill = "white"))  + # Set the entire plot background to white)
  ggtitle("MDM2") +
  xlab("") +
  ylab("log2(TPM)") +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ALT")
ggsave("MDM2_BAT.dpi400.pdf", width = 4, height = 3.5, dpi=400)