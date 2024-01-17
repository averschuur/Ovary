### FIGURES ###

### Library ----------------------------------------------------------------------
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


### Heatmap ovary ----------------------------------------------------------------------
sample_cor <- cor(betas_topvar_ovary)
rownames(sample_cor) <- anno_ovary$tumorType

anno_ovary <- anno_ovary %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' ~ 'primary_not_in_teratoma',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary'))


annotation = anno_ovary$annotation

row_ha = rowAnnotation(
  tumortype = annotation)

pdf(file="heatmap_ovary_20240105.pdf", width=9, height=8)
Heatmap(sample_cor, clustering_distance_rows = "euclidean", show_column_names = FALSE, show_row_names = FALSE, right_annotation = row_ha)
dev.off()

### UMAP ovary -------------------------------------------------------------------------
anno_ovary %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = annotation), shape=19, size = 3) +
  #geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Tumor Type") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.ticks.length.y =unit(.05, "cm"),
    legend.position = "bottom",
    #legend.margin=margin(0,0,0,0),
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("UMAP_ovary_only_annotation_20240115.pdf", path= "./plots/", dpi=500)


### Heatmap ovary + ileum ----------------------------------------------------------------------
sample_cor <- cor(betas_topvar_ovary_ileal)
rownames(sample_cor) <- anno_ovary_ileal$tumorType

anno_ovary_ileal <- anno_ovary_ileal %>%
  mutate(localization = case_when(tumorType == 'ovary' ~ 'ovary',
                                  tumorType == 'ilealNET' ~ 'ileal')) %>%
  mutate(type = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                          location == 'primary' ~ 'primary',
                          location == 'metastasis_midgut' ~ 'metastasis',
                          location == 'metastasis_rectum' ~ 'metastasis',
                          location == 'metastasis_pancreas' ~ 'metastasis',
                          location == 'metastasis' ~ 'metastasis'))

localization = anno_ovary_ileal$localization
type = anno_ovary_ileal$type

row_ha = rowAnnotation(
  type = type, 
  localization = localization)


pdf(file="heatmap_ovary_ileal_20240105.pdf", width=10, height=8)
Heatmap(sample_cor, show_column_names = FALSE, show_row_names = FALSE, right_annotation = row_ha)
dev.off()



### UMAP ovary + ileum -------------------------------------------------------------------------

anno_ovary_ileal %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = annotation), shape=19, size = 3) +
  #geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Tumor Type") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.ticks.length.y =unit(.05, "cm"),
    legend.position = "bottom",
    #legend.margin=margin(0,0,0,0),
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("UMAP_ovary_ileal_annotation_20240105.png", path= "./plots/", dpi=500)



### Heatmap ovary + ileum +rectum + pancreas ----------------------------------------------------------------------
sample_cor <- cor(betas_topvar_ovary_ileal_pancreatic_rectal)
rownames(sample_cor) <- anno_ovary_ileal_pancreatic_rectal$tumorType

anno_ovary_ileal_pancreatic_rectal <- anno_ovary_ileal_pancreatic_rectal %>%
  mutate(localization = case_when(tumorType == 'ovary' ~ 'ovary',
                                  tumorType == 'ilealNET' ~ 'ileal',
                                  tumorType == 'panNET' ~ 'pancreatic',
                                  tumorType == 'pulmNET' ~ 'pulmonary',
                                  tumorType == 'rectalNET' ~ 'rectal')) %>%
  mutate(type = case_when(location == 'primary_teratoma' ~ 'PONWT',
                          location == 'primary' & tumorType == "ovary" ~ 'PONNT',
                          location == 'primary' ~ 'primary',
                          location == 'metastasis_midgut' ~ 'NOM',
                          location == 'metastasis_rectum' ~ 'NOM',
                          location == 'metastasis_pancreas' ~ 'NOM',
                          location == 'metastasis' ~ 'metastasis'))

localization = anno_ovary_ileal_pancreatic_rectal$localization
type = anno_ovary_ileal_pancreatic_rectal$type

row_ha = rowAnnotation(
  localization = localization,
  type = type)


pdf(file="heatmap_ovary_ileal_pancreatic_rectal_20240115.pdf", width=10, height=8)
Heatmap(sample_cor, show_column_names = FALSE, show_row_names = FALSE, right_annotation = row_ha)
dev.off()


### UMAP ovary + ileum + pancreas + rectum -------------------------------------------------------------------------

anno_ovary_ileal_pancreatic_rectal %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(col = localization, shape = type), size = 3) +
  scale_shape_manual(values=c(15, 19, 17, 18))+
  #geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Tumor Type") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    axis.ticks.length.y =unit(.05, "cm"),
    legend.position = "bottom",
    #legend.margin=margin(0,0,0,0),
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("UMAP_ovary_ileal_pancreatic_rectal_20240115.pdf", path= "./plots/", dpi=500)
