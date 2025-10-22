

### read in sil score
sil = fread("/Users/jhofvand/Dropbox/BMC/RNAseq_Fresh/clustering_all_samples/Annottaions_Figge_2024/tSNEA_704/Silhouette_per_diagnosis.tsv") %>% dplyr::mutate(lab_no = strsplit(UID, "_") %>% map_chr(., 1))
### read in clinical data
c = read_xlsx("~/Dropbox/Jakob+Figge/Transcriptomics/Manuscript files/Clinical data_250502.xlsx") %>% dplyr::select(rowname, Diagnosisa, Class, Driverb) %>%
  drop_na() %>% dplyr::mutate(lab_no = strsplit(rowname, "_") %>% map_chr(., 1))

c = inner_join(sil, c)

colnames(c) = c("UID", "Sil", "Cluster", "Genome", "lab_no", "rowname", "Diagnosisa", "Class", "Driverb" )
??stat_compare_means
p = c %>%
  group_by(Cluster) %>%
  mutate(number = n()) %>%
  filter(number > 2) %>%
  mutate(median = median(Sil)) %>%
  ggplot(aes(reorder(Cluster, median) , Sil, fill = Driverb)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = c("snow4","aquamarine3", "darkred", "azure3", "blue")) +
  coord_flip()

ggsave(p, file="Boxplot_sil_score.pdf", width = 6, height = 10, dpi = 600)

###

sil = fread("/Users/jhofvand/Dropbox/BMC/RNAseq_Fresh/clustering_all_samples/Annottaions_Figge_2024/tSNEA_704/Silhouette_per_diagnosis.tsv") %>% dplyr::mutate(lab_no = strsplit(UID, "_") %>% map_chr(., 1))
colnames(sil) = c("UID", "Sil", "Cluster",  "Genome", "lab_no")
c = read_xlsx("~/Dropbox/Jakob+Figge/Transcriptomics/Manuscript files/Clinical_data_250515.xlsx") %>% dplyr::select(rowname, Simple_Complex) %>%
  drop_na() %>% dplyr::mutate(lab_no = strsplit(rowname, "_") %>% map_chr(., 1))

inner_join(c, sil) %>%
  filter(Simple_Complex != "NA") %>%
  ggplot(aes(Simple_Complex, Sil)) +
  geom_boxplot() + theme_bw() +
  stat_compare_means(method = "t.test", label.x = 1.3, size = 4) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_rect(fill = "white"),  # Set background to white
    plot.background = element_rect(fill = "white")   # Set the entire plot background to white
  )
ggsave("Boxplot_Sil_Simple_Complex.pdf", height = 4.5, width = 5, dpi = 300)
