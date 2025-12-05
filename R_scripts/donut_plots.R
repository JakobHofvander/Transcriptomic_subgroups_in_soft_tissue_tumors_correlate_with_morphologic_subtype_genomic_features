# read in Silhouette score
sil_data = fread("source_data/Silhouette_per_diagnosis.tsv") %>% 
  group_by(Cluster) %>% 
  mutate(median_sil_in_cluster = median(Silhouette_Score)) %>% 
  dplyr::mutate(lab_no = strsplit(UID, "_") %>% map_chr(., 1)) %>% 
  dplyr::select(UID, lab_no, median_sil_in_cluster) %>% 
  distinct()

# add annotation and remove groups with few samples
sil_data = fread("source_data/meta_data.txt") %>% 
  inner_join(sil_data, by = c("lab_no" = "lab_no")) %>% 
  group_by(Cluster) %>% mutate(n = n()) %>% filter(n > 2) %>% 
  mutate(median_sil_in_cluster = median(Silhouette_Score)) %>% 
  dplyr::select(Cluster, Driver, median_sil_in_cluster) %>% distinct()

# create plot_data_high based on median value
plot_data_high = sil_data %>% filter(median_sil_in_cluster > 0)  %>% 
  arrange(median_sil_in_cluster) %>% 
  group_by(Driver) %>% 
  mutate(n = n()) %>% 
  distinct(Driver, n) %>% 
  ungroup() %>% 
  mutate(f = n / sum(n))

# create plot_data_low based on median value
plot_data_low = sil_data %>% filter(median_sil_in_cluster < 0) %>% 
  arrange(median_sil_in_cluster) %>% 
  group_by(Driver) %>% 
  mutate(n = n()) %>% 
  distinct(Driver, n) %>% 
  ungroup() %>% 
  mutate(f = n / sum(n))

## create low group Donut plot
# Compute the cumulative percentages (top of each rectangle)
plot_data_low$ymax <- cumsum(plot_data_low$f)
# Compute the bottom of each rectangle
plot_data_low$ymin <- c(0, head(plot_data_low$ymax, n=-1))
plot_data_low$labelPosition <- (plot_data_low$ymax + plot_data_low$ymin) / 2

ggplot(plot_data_low, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2, fill=Driverb)) +
  geom_rect() +
  geom_text(x=2.5, aes(y=labelPosition, label=f), size=6, position = position_dodge(width = 2)) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("darkgrey", "darkred")) +
  coord_polar(theta="y") +
  xlim(c(1, 6)) +
  theme_void() +
  theme() 
ggsave("low_sil_donut.pdf", bg = "transparent",  width = 6, height = 5, dpi = 300)

## create high group Donut plot
# Compute the cumulative percentages (top of each rectangle)
plot_data_high$ymax <- cumsum(plot_data_high$f)
# Compute the bottom of each rectangle
plot_data_high$ymin <- c(0, head(plot_data_high$ymax, n=-1))

plot_data_high$labelPosition <- (plot_data_high$ymax + plot_data_high$ymin) / 2

ggplot(plot_data_high, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2, fill=Driverb)) +
  geom_rect() +
  geom_text(x=2.5, aes(y=labelPosition, label=f), size=6, position = position_dodge(width = 2)) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("darkgrey", "aquamarine3", "darkred", "blue")) +
  coord_polar(theta="y") +
  xlim(c(1, 6)) +
  theme_void() +
  theme() 
ggsave("high_sil_donut.pdf", bg = "transparent",  width = 6, height = 5, dpi = 300)


