#### ADDITIONAL ANALYSIS

# open libraries
library("openxlsx")
library("writexl")
library(tidyr)
library(naniar)

### CNV -------------------------------------------------------------------------
CNV <- read.xlsx("./input/CNV analysis.xlsx")

Chr18 <- CNV %>%
  select(arrayId, Chr18)

anno <- merge(anno, Chr18, by = "arrayId")
anno <- merge(anno, Chr18, by = "arrayId")

# plot umap
anno %>% 
  ggplot(aes(umap_x, umap_y, col = Chr18)) +
  geom_point(size = 4) +
  geom_text(aes(label = annotation), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

### IHC ----------------------------------------------------------------------
IHC <- read.xlsx("./input/IHC.xlsx")

IHC <- IHC %>%
  drop_na(sampleName)

IHC <- IHC[,c(2,10, 13:19)]
colnames(IHC) <- c("sampleName","Groeiwijze","SSTR","CDX2","PAX8","ISLET1","ARX","TTF1","SATB2")

anno <- left_join(anno, IHC, by = "sampleName")

na_strings <- c("NA", "N A", "N / A", "N/A", "N/ A", "Not Available", "NOt available")
anno <- anno %>%
  replace_with_na_all(condition = ~.x %in% na_strings)

write.xlsx(anno, "./input/anno_IHC.xlsx")

# open data
anno <- read.xlsx("./input/anno_IHC.xlsx")

# plot umap
anno %>% 
  ggplot(aes(umap_x, umap_y, col = Groeiwijze)) +
  geom_point(size = 4) +
  geom_text(aes(label = annotation), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

anno %>% 
  ggplot(aes(umap_x, umap_y, col = CDX2)) +
  geom_point(size = 4) +
  geom_text(aes(label = annotation), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

anno %>% 
  ggplot(aes(umap_x, umap_y, col = ISLET1)) +
  geom_point(size = 4) +
  geom_text(aes(label = annotation), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

anno_IHC %>% 
  ggplot(aes(umap_x, umap_y, col = ARX)) +
  geom_point(size = 4) +
  geom_text(aes(label = annotation), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

anno %>% 
  ggplot(aes(umap_x, umap_y, col = TTF1)) +
  geom_point(size = 4) +
  geom_text(aes(label = annotation), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

anno %>% 
  ggplot(aes(umap_x, umap_y, col = SATB2)) +
  geom_point(size = 4) +
  geom_text(aes(label = annotation), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

anno %>% 
  ggplot(aes(umap_x, umap_y, col = SSTR)) +
  geom_point(size = 4) +
  geom_text(aes(label = annotation), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

