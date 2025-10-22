#### boxplot for marker genes in SA lei
setwd("~/Dropbox/Jakob+Figge/Transcriptomics/source_data_and_code/")

df_pc = fread("tpm_matrix.tsv")

ano = fread("meta_data2.txt")


markers = c("ACTG2", "SLMAP", "LMOD1", "CFL2", "MYLK", "ARL4C", "MYOCD")

## select SA lei samples
plot_df = df_pc %>% filter(symbol %in% markers) %>% 
  pivot_longer(!symbol) %>% 
  mutate(log2TPM = log2(value + 0.01)) %>% 
  inner_join(ano, by = c("name" = "lab_no")) %>% 
  filter(Cluster == "SA lei A" | Cluster == "SA lei")

#### plot them
for (g in markers)
{plot=plot_df %>% 
  filter(symbol == g) %>% 
  ggplot(aes(Cluster, log2TPM), fill = Cluster) +
  geom_jitter(size = 1.5, width = 0.2) +
  stat_summary(fun = mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, color="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.background = element_rect(fill = "white"),  # Set background to white
        plot.background = element_rect(fill = "white"))  + # Set the entire plot background to white)
  ggtitle(g) +
  xlab("") +
  ylab("log2(TPM)") +
  stat_compare_means(method = "t.test", label.x = 1.3, size = 4)
ggsave(paste(g,"SA_lei.dpi400.pdf",sep="."), plot = plot, bg = "white", width = 4, height = 3.5, dpi=400)
}


