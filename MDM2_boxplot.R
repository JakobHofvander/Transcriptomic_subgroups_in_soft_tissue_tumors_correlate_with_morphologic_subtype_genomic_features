##MDM2 boxplot in BAT tumors
# select samples and genes
markers = c("MDM2")

plot_df = df_pc %>% filter(symbol %in% markers) %>% 
  pivot_longer(!symbol) %>% 
  mutate(log2TPM = log2(value + 0.01)) %>% 
  inner_join(q, by = c("name" = "UID")) %>% 
  filter(NewDiagnosisFigge_short == "ALT" | NewDiagnosisFigge_short == "Hib" | NewDiagnosisFigge_short == "Lip" | NewDiagnosisFigge_short == "Lip spindle" | NewDiagnosisFigge_short == "Lipblast")

# MDM2 boxplot
plot_df %>% 
  filter(symbol == "MDM2") %>% 
  ggplot(aes(NewDiagnosisFigge_short, log2TPM), fill = Cluster) +
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

################## SA lei marler genes ######
# select samples and genes
q = read_excel("~/Dropbox/Jakob+Figge/Transcriptomics/Sil_score_v2.xlsx") %>% 
  dplyr::select(UID, Cluster)

markers = c("ACTG2", "SLMAP", "LMOD1", "CFL2", "MYLK", "ARL4C", "MYOCD")

plot_df = df_pc %>% filter(symbol %in% markers) %>% 
  pivot_longer(!symbol) %>% 
  mutate(log2TPM = log2(value + 0.01)) %>% 
  inner_join(q, by = c("name" = "UID")) %>% 
  filter(Cluster == "SA lei A" | Cluster == "SA lei")

# boxplot for all marker genes
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



