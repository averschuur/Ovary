# libraries

# Load required packages and sources
library(minfi)
library(RFpurify)

library(tidyverse)

library(Rtsne)
library(umap)

library(caret)
library(pheatmap)

source("./scripts/0_helpers.R")



## load sample annotation and filter -------------------------------------------

anno <- readRDS("./input/sample_annotation.rds")


# load unfiltered beta values
betas <- readRDS("./input/betas_everything.rds")
betas <- betas[, anno$arrayId]

## get purity estimates --------------------------------------------------------
absolute <- RFpurify::predict_purity_betas(betas = betas, method = "ABSOLUTE")
estimate <- RFpurify::predict_purity_betas(betas = betas, method = "ESTIMATE")


## obtain bilsufite conversion scores-----------------------------------------------------
raw_epic<- readRDS(file = "./input/raw_epic.rds")

# investigate QC measures 1: bisulfite conversion efficiency 
conversion <- getControlBeta(raw_epic, controls = "BISULFITE CONVERSION II") %>% 
  as_tibble() %>% 
  group_by(arrayId) %>% 
  summarise(conversion = min(value[channel == "Red"]) / max(value[channel == "Green"]))

rm(raw_epic)

# merge conversion scores with annotation
anno <- anno %>% 
  left_join(conversion)

## collect data into tibble ----------------------------------------------------
anno <- anno %>% 
  mutate(absolute = absolute,
         estimate = estimate,
         avg_beta_unfiltered = apply(betas, 2, mean, na.rm = TRUE))


### plot basic statistics for the data set ---------------------------------------

# overview of all variables
anno %>% 
  select(-c(arrayId, source, sampleName, location)) %>% 
  GGally::ggpairs()


# number of cases
anno %>% 
  group_by(tumorType,location) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(fill = location, y =n, x = tumorType)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = "Number of Cases")
ggsave("Figure_count tumorType per location.png", path= "./plots/", dpi=500)

# number of cases by source
anno %>% 
  ggplot(aes(tumorType)) +
  geom_bar(aes(fill = source)) +
  labs(title="Tumor Type By Source",x = "", y = "No. of cases") +
  theme_bw(base_size = 18)
ggsave("Figure_count tumorType by source.png", path= "./plots/", dpi=500)

# tumor purity
anno %>% 
  ggplot(aes(absolute, estimate)) +
  geom_point(size = 3) +
  geom_smooth() +
  theme_bw(base_size = 18) +
  labs(x = "Purity (ABSOLUTE)", y = "Purity (ESTIMATE)")

anno %>% 
  ggplot(aes(tumorType, estimate, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(0.5), alpha=0.7) +
  #geom_jitter(aes(fill=tumorType, shape=tumorType)) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Figure_estimate tumorpurity.png", path= "./plots/", dpi=500)

anno %>% 
  ggplot(aes(tumorType, absolute, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(0.5), alpha=0.9) +
  #geom_jitter(aes(fill=tumorType, shape=tumorType)) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Figure_absolute tumorpurity.png", path= "./plots/", dpi=500)

# average methylation
anno %>% 
  ggplot(aes(tumorType, avg_beta_unfiltered, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.7) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Figure_average methylation.png", path= "./plots/", dpi=500)

# conversion scores
anno %>% 
  ggplot(aes(tumorType, conversion, col = tumorType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.7) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Figure_conversion.png", path= "./plots/", dpi=500)

# save anno
saveRDS(object = anno, file = "./input/sample_annotation.rds")

# clean up
rm(absolute, estimate, betas)



### Ovary only --------------------------------------------------------------------------------------------------------------

# load filtered beta values
betas <- readRDS("./input/betas_filtered.rds")

anno_ovary <- anno %>%
  filter(tumorType == "ovary")
anno_ovary <- anno_ovary[c(1:9, 11:18),]

betas_ovary <- betas[, anno_ovary$arrayId]
probe_var_ovary <- apply(betas_ovary, 1, var)
probes_topvar_ovary <- rownames(betas_ovary)[order(probe_var_ovary, decreasing = TRUE)]
#saveRDS(object = probes_topvar_ovary, file = "./output/top_variable_probes_ovary_20240105.rds")

# pick betas for 5,000 top variable probes
betas_topvar_ovary <- betas_ovary[probes_topvar_ovary[1:5000], ]


# heatmap of sample-wise correlations
sample_cor <- cor(betas_topvar_ovary)

rownames(sample_cor) <- anno_ovary$tumorType
pheatmap::pheatmap(sample_cor, labels_row = anno_ovary$sampleName, show_colnames = FALSE)

# run UMAP
set.seed(45098)
umap_settings <- umap.defaults
umap_settings$n_neighbors = 15
umap_settings$min_dist = 0.1

umap<- umap(d = t(betas_topvar_ovary), config = umap_settings, ret_model = TRUE)

anno_ovary <- anno_ovary %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2],
         Label = sampleName)

#saveRDS(object = umap, file = "./output/umap_model_ovary_20240105.rds")
#saveRDS(object = anno_ovary, file = "./output/sample_annotation_umap_purity_ovary_20240105.rds")


# plot UMAP
anno_ovary %>% 
  ggplot(aes(umap_x, umap_y, col = location)) +
  geom_point(size = 4) +
  geom_text(aes(label = Label), size = 4, nudge_y = 0.1) +
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_only_annotation_20240105.png", path= "./plots/")



## Ovary + ileal NET  ----------------------------------------------------------------

# load filtered beta values
betas <- readRDS("./input/betas_filtered.rds")


anno_ovary_ileal <- anno %>%
  filter(tumorType == "ovary" | tumorType == "ilealNET")
anno_ovary_ileal <- anno_ovary_ileal[c(1:9, 11:21,23:29),]
betas_ovary_ileal <- betas[, anno_ovary_ileal$arrayId]

# determine most variable probes across dataset and subset beta values
probe_var_ovary_ileal <- apply(betas_ovary, 1, var)
probes_topvar_ovary_ileal <- rownames(betas_ovary_ileal)[order(probe_var_ovary_ileal, decreasing = TRUE)]
#saveRDS(object = probes_topvar_ovary_ileal, file = "./output/top_variable_probes_ovary_ileal_20240105.rds")

# pick betas for 5,000 top variable probes
betas_topvar_ovary_ileal <- betas_ovary_ileal[probes_topvar_ovary_ileal[1:5000], ]


# heatmap of sample-wise correlations
sample_cor <- cor(betas_topvar_ovary_ileal)

rownames(sample_cor) <- anno_ovary_ileal$tumorType
pheatmap::pheatmap(sample_cor, labels_row = anno_ovary_ileal$sampleName, show_colnames = FALSE)



# run UMAP
set.seed(45098)
umap_settings <- umap.defaults
umap_settings$n_neighbors = 5
umap_settings$min_dist = 0.2

umap<- umap(d = t(betas_topvar_ovary_ileal), config = umap_settings, ret_model = TRUE)

anno_ovary_ileal <- anno_ovary_ileal %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])

#saveRDS(object = umap, file = "./output/umap_model_ovary_ileal_20240105.rds")
#saveRDS(object = anno_ovary_ileal, file = "./output/sample_annotation_umap_purity_ovary_ileal_20240105.rds")

anno_ovary_ileal <- anno_ovary_ileal %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' ~ 'primary',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary',
                                location == 'metastasis' ~ 'metastasis'))

# plot UMAP
anno_ovary_ileal %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType, shape = annotation)) +
  geom_point(size = 4) +
  #geom_text(aes(label = sampleName), size = 4) +
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_only.png", path= "./plots/")




## all --------------------------------------------------------------------------

anno_ovary_ileal_pancreatic_pulm_rectal <- anno %>%
  filter(tumorType == "ovary" | tumorType == "ilealNET" | tumorType == "rectalNET" | tumorType == "panNET" | tumorType == "pulmNET")
anno_ovary_ileal_pancreatic_pulm_rectal <- anno_ovary_ileal_pancreatic_pulm_rectal[c(1:9, 11:63,66:71),]
betas_ovary_ileal_pancreatic_pulm_rectal <- betas[, anno_ovary_ileal_pancreatic_pulm_rectal$arrayId]

# determine most variable probes across dataset and subset beta values
probe_var_ovary_ileal_pancreatic_pulm_rectal <- apply(betas_ovary_ileal_pancreatic_pulm_rectal, 1, var)
probes_topvar_ovary_ileal_pancreatic_pulm_rectal <- rownames(betas_ovary_ileal_pancreatic_pulm_rectal)[order(probe_var_ovary_ileal_pancreatic_pulm_rectal, decreasing = TRUE)]
saveRDS(object = probes_topvar_ovary_ileal_pancreatic_pulm_rectal, file = "./output/top_variable_probes_ovary_ileal_pancreatic_pulm_rectal_200240105.rds")

# pick betas for 5,000 top variable probes
betas_topvar_ovary_ileal_pancreatic_pulm_rectal <- betas_ovary_ileal_pancreatic_pulm_rectal[probes_topvar_ovary_ileal[1:5000], ]


# heatmap of sample-wise correlations
sample_cor <- cor(betas_topvar_ovary_ileal_pancreatic_pulm_rectal)

rownames(sample_cor) <- anno_ovary_ileal_pancreatic_pulm_rectal$tumorType
pheatmap::pheatmap(sample_cor, labels_row = anno_ovary_ileal_pancreatic_pulm_rectal$sampleName, show_colnames = FALSE)


# run UMAP
set.seed(45098)
umap_settings <- umap.defaults
umap_settings$n_neighbors = 10
umap_settings$min_dist = 0.2

umap<- umap(d = t(betas_topvar_ovary_ileal_pancreatic_pulm_rectal), config = umap_settings, ret_model = TRUE)

anno_ovary_ileal_pancreatic_pulm_rectal <- anno_ovary_ileal_pancreatic_pulm_rectal %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])

#saveRDS(object = umap, file = "./output/umap_model_ovary_ileal_pancreatic_pulm_rectal.rds")
#saveRDS(object = anno, file = "./output/sample_annotation_umap_purity_ovary_ileal_pancreatic_pulm_rectal.rds")

anno_ovary_ileal_pancreatic_pulm_rectal <- anno_ovary_ileal_pancreatic_pulm_rectal %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' ~ 'primary',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary',
                                location == 'metastasis' ~ 'metastasis'))

# plot UMAP
anno_ovary_ileal_pancreatic_pulm_rectal %>% 
  ggplot(aes(umap_x, umap_y, col = tumorType, shape = annotation)) +
  geom_point(size = 4) +
  geom_text(aes(label = sampleName), size = 4) +
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave("UMAP_ovary_ileal_pancreatic_pulm_rectal.png", path= "./plots/")





