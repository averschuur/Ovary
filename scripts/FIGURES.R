### FIGURES ###

### Library ----------------------------------------------------------------------
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


### Heatmap ovary ----------------------------------------------------------------------
set.seed(123)
mat = betas_topvar

localization = anno$annotation
column_ha = HeatmapAnnotation(
  localization = localization)

pdf(file="heatmap_ovary_ileal_pancreatic_rectal_euclidean_20240123.pdf", width=10, height=8)
Heatmap(mat, name = "mat", 
        show_column_names = FALSE, show_row_names = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        top_annotation = column_ha)
dev.off()

### UMAP  -------------------------------------------------------------------------
anno %>% 
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
ggsave("UMAP_annotation_20240121.pdf", path= "./plots/", dpi=500)



### Chr 18 loss projected on UMAP------------------------------------------------------------------

#load data
CNV <- read.xlsx("./input/CNV analysis.xlsx")

Chr18 <- CNV %>%
  select(arrayId, Chr18)

anno <- merge(anno, Chr18, by = "arrayId")

# plot umap
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = Chr18), shape=19, size = 3) +
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
ggsave("UMAP_chr18 loss_20240121.pdf", path= "./plots/", dpi=500)


### IHC projected on UMAP------------------------------------------------------------------

#load data
anno <- read.xlsx("./input/anno_IHC.xlsx")

# plot Groeiwijze projected over umap
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = Groeiwijze), shape=19, size = 3) +
  #geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "Groeiwijze") +
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
ggsave("UMAP_Groeiwijze_20240123.pdf", path= "./plots/", dpi=500)

# plot CDX2 projected over umap
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = CDX2), shape=19, size = 3) +
  #geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "CDX2") +
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
ggsave("UMAP_CDX2_20240123.pdf", path= "./plots/", dpi=500)

# plot TTF1 projected over umap
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = TTF1), shape=19, size = 3) +
  #geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "TTF1") +
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
ggsave("UMAP_TTF1_20240123.pdf", path= "./plots/", dpi=500)


# plot SATB2 projected over umap
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = SATB2), shape=19, size = 3) +
  #geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "SATB2") +
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
ggsave("UMAP_SATB2_20240123.pdf", path= "./plots/", dpi=500)


# plot ARX projected over umap
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = ARX), shape=19, size = 3) +
  #geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "ARX") +
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
ggsave("UMAP_ARX_20240121.pdf", path= "./plots/", dpi=500)

# plot ISLET 1 projected over umap
anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = ISLET1), shape=19, size = 3) +
  #geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "ISLET1") +
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
ggsave("UMAP_ISLET1_20240121.pdf", path= "./plots/", dpi=500)



