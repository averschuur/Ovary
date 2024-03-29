# last edited 04/01/2023 by AV Verschuur

### Load required packages and sources  ----------------------------------------------

library(tidyverse)
library(minfi)
library(paletteer)
library(Rtsne)

source("./scripts/0_helpers.R")


#### load data and annotation 
raw_epic<- readRDS(file = "./input/raw_epic.rds")
anno <- readRDS("./input/sample_annotation.rds")

# investigate QC measures 1: bisulfite conversion efficiency 
conv_raw_epic <- getControlBeta(raw_epic, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

rm(raw_epic)


## merge conversion scores with annotation
anno <- anno %>% 
  left_join(conv_raw_epic)


# add classification of conversion tot anno
anno <- anno %>%
  mutate(good_score1 = case_when(conversion >= 1 ~ 'good',
                                 conversion >= 0.8 ~ 'medium',
                                 conversion < 0.8 ~ 'bad'))

# stats
anno %>%
  group_by(good_score1) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)*100)

count <- anno %>%
  group_by(source, good_score1) %>%
  count(good_score1)

count %>%
  ggplot(aes(fill = good_score1, n, source)) +
  geom_bar(position="fill", stat="identity")
ggsave("Figure s1_Control Scores by source.pdf", path= "./output/", dpi=500)


# save anno ex bad and medium samples
anno_ebs <- anno %>%
  subset(good_score1 == "good")

saveRDS(anno_ebs, file ="./annotation/sample_annotation_conversion_scoresex_bad_samples.csv")

# Prepare for data visualisation
anno %>%
  mutate(reply = anno$arrayId %in% colnames(betas)) %>%
  group_by(reply, sampleName) %>%
  summarize(reply_sum = sum(reply))

anno <- anno %>% 
  filter(!sampleName == "UMCU_ACC2")

anno1 <- anno %>% 
  filter(tumorType %in% c("PanNET", "ACC", "SPN", "PDAC", "normal", "acc normal", "PanNEC", "PB")) %>% 
  filter(location %in% c("primary", "pancreas")) %>% 
  filter(source != "UMCU")



### create UMAP:

anno %>% 
  ggplot(aes(umap_x, umap_y)) + 
  geom_point(aes(color = tumorType), shape=19, size = 4) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe",
                                     limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                                     breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                                     labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
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
    legend.margin=margin(0,0,0,0),
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid = element_blank()) +
  guides(col = guide_legend(nrow = 1))
ggsave("Figure 1F_UMAP-all_colorSafe.pdf", path= "./output/", dpi=500)


# Figure 1S: UMAP conversion
anno1 %>% 
  ggplot(aes(umap_x, umap_y, col = conversion)) + 
  geom_point(shape=19, size = 5) +
  paletteer::scale_color_paletteer_c("ggthemes::Green-Blue-White Diverging") +
  labs(#title="UMAP Clustering",
    x = "UMAP 1", 
    y = "UMAP 2",
    col = "ControlScores") +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 12),
    axis.title.x = element_text(face = "bold", hjust = 0.5, vjust=0, size = 10),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 10),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    legend.margin=margin(0,0,0,0),
    legend.title = element_text(face = "bold", size = 10, vjust = 0.5),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, 'cm'),
    panel.grid = element_blank()) +
  guides(colour=guide_colourbar(barwidth=15,label.position="bottom", title.vjust=0.9, label.vjust=2))
ggsave("Figure S1_UMAP-conversion.pdf", path= "./output/", dpi=500)

anno1 %>% 
  ggplot(aes(tumorType, conversion, fill=tumorType, col=tumorType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=tumorType), position=position_jitter(0.2),shape = 21,size =2) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe", breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")) +
  labs(#title="Absolute Tumor Purity By Tumor Type",
    x = "", 
    y = "Control Scores") +
  geom_hline(yintercept=0.8, color = "grey", linetype = "dashed") +
  scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
                   labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure s1_tumorType_Control Scores.pdf", path= "./output/", dpi=500)

anno1 %>% 
  ggplot(aes(source, conversion, fill=source, col=source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill=source), position=position_jitter(0.2),shape = 21,size =2) +
  paletteer::scale_fill_paletteer_d("rcartocolor::Safe", 
                                    #breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
  ) +
  paletteer::scale_color_paletteer_d("rcartocolor::Safe", 
                                     #breaks = c("ACC", "normal", "PanNEC", "PanNET", "PB", "PDAC", "SPN")
  ) +
  labs(#title="Absolute Tumor Purity By Tumor Type",
    x = "", 
    y = "Control Scores") +
  geom_hline(yintercept=0.8, color = "grey", linetype = "dashed") +
  #scale_x_discrete(limits = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
  #breaks = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "normal"),
  #labels = c("PDAC", "ACC", "PanNET","PanNEC", "SPN", "PB", "NORMAL")) +
  theme_bw(base_size = 18) +
  theme(#plot.title = element_text(face = "bold", hjust = 0, size = 14),
    axis.title.y = element_text(face = "bold", hjust = 0.5, vjust=2, size = 12),
    axis.text.x = element_text(face = "bold", size = 11, angle = 90, vjust=0.5, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid = element_blank())
ggsave("Figure s1_Control Scores by study.pdf", path= "./output/", dpi=500)