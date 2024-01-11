library("readxl")
library("writexl")
library(tidyr)
library(naniar)

# prepare data ----------------------------------------------------------------------
#IHC <- read_excel("./input/IHC.xlsx")

#IHC <- IHC %>%
#  drop_na(sampleName)


#anno1 <- anno

#IHC1 <- IHC
#IHC1 <- IHC1[,c(2,6,8,10:11, 13:19)]
#IHC1 <- IHC1[1:16,]

#anno_IHC <- left_join(anno1, IHC1)
#colnames(anno_IHC) <- c("arrayId","source","sampleName","arrayType","tumorType","location",
#                        "avg_beta_filtered","umap_x","umap_y","Label","Location","Origin",
#                        "Groeiwijze","Teratoom","SSTR","CDX2","PAX8","ISLET1","ARX","TTF1",
#                        "SATB2") 

#na_strings <- c("NA", "N A", "N / A", "N/A", "N/ A", "Not Available", "NOt available")
#anno_IHC <- anno_IHC %>%
#  replace_with_na_all(condition = ~.x %in% na_strings)

#write_xlsx(anno_IHC, "./input/anno_IHC.xlsx")

## open data -------------------------------------------------------------------------
anno_IHC <- read.xlsx("./input/anno_IHC.xlsx")

### plot UMAP -----------------------------------------------------------------------
anno_IHC %>% 
  ggplot(aes(umap_x, umap_y, col = Groeiwijze)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_ileal_pancreas_rectum_Groeiwijze.png", path= "./plots/")

anno_IHC %>% 
  ggplot(aes(umap_x, umap_y, col = CDX2)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_ileal_pancreas_rectum_CDX2.png", path= "./plots/")

anno_IHC %>% 
  ggplot(aes(umap_x, umap_y, col = ISLET1)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_ileal_pancreas_rectum_ISLET1.png", path= "./plots/")

anno_IHC %>% 
  ggplot(aes(umap_x, umap_y, col = ARX)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_ileal_pancreas_rectum_ARX.png", path= "./plots/")

anno_IHC %>% 
  ggplot(aes(umap_x, umap_y, col = TTF1)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_ileal_pancreas_rectum_TTF1.png", path= "./plots/")

anno_IHC %>% 
  ggplot(aes(umap_x, umap_y, col = SATB2)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_ileal_pancreas_rectum_SATB2.png", path= "./plots/")

anno_IHC %>% 
  ggplot(aes(umap_x, umap_y, col = SSTR)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_ileal_pancreas_rectum_SSTR.png", path= "./plots/")

anno_IHC %>% 
  ggplot(aes(umap_x, umap_y, col = Chr_18_loss)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), size = 4, nudge_x = 0.8) +
  scale_shape_manual(values=c(18, 15, 19, 17, 0, 1))+
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_ileal_pancreas_rectum_Chr_18_loss.png", path= "./plots/")
