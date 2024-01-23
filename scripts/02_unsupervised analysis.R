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

anno <- readRDS("./input/sample.rds")

anno <- anno %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'PONWT',
                                location == 'primary' & tumorType == "ovary" ~ 'PONNT',
                                location == 'metastasis_midgut' ~ 'NOM',
                                location == 'metastasis_rectum' ~ 'NOM',
                                location == 'metastasis_pancreas' ~ 'NOM',
                                location == 'primary' & tumorType == 'ilealNET' ~ 'ileal',
                                location == 'metastasis'& tumorType == 'ilealNET' ~ 'ileal_metastasis',
                                location == 'primary' & tumorType == 'rectalNET' ~ 'rectal',
                                location == 'metastasis'& tumorType == 'rectalNET' ~ 'rectal_metastasis',
                                location == 'primary' & tumorType == 'panNET' ~ 'pancreatic'))


## load unfiltered beta values -------------------------------------------
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


## plot basic statistics for the data set ---------------------------------------

# load anno and select study cases
anno <- readRDS("./input/sample_annotation_20240121.rds")

anno <- anno %>%
  filter(location == "primary" | location == "primary_teratoma" | location == "metastasis_midgut" | location == "metastasis_rectum" | location == "metastasis_pancreas")
anno <- anno[c(1:9, 11:27, 29:38, 40:55, 58:60),] # remove NEC and duplicate cases

# number of cases
anno %>% 
  group_by(annotation) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(fill = annotation, y = n, x = annotation)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = "Number of Cases")

# number of cases by source
anno %>% 
  ggplot(aes(annotation)) +
  geom_bar(aes(fill = source)) +
  labs(title="Tumor Type By Source",x = "", y = "No. of cases") +
  theme_bw(base_size = 18)

# tumor purity
anno %>% 
  ggplot(aes(absolute, estimate)) +
  geom_point(size = 3) +
  geom_smooth() +
  theme_bw(base_size = 18) +
  labs(x = "Purity (ABSOLUTE)", y = "Purity (ESTIMATE)")

anno %>% 
  ggplot(aes(annotation, estimate, col = annotation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(0.5), alpha=0.7) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

anno %>% 
  ggplot(aes(annotation, absolute, col = annotation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(0.5), alpha=0.9) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# average methylation
anno %>% 
  ggplot(aes(annotation, avg_beta_unfiltered, col = annotation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.7) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# conversion scores
anno %>% 
  ggplot(aes(annotation, conversion, col = annotation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.7) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# save anno
#saveRDS(object = anno, file = "./input/sample_annotation_20240121.rds")

# clean up
rm(absolute, estimate, betas)



## Unsupervised analysis --------------------------------------------------------------------------

# load filtered beta values
betas <- readRDS("./input/betas_filtered.rds")

# load anno and select study cases
anno <- readRDS("./input/sample_annotation_20240121.rds")

anno <- anno %>%
  filter(location == "primary" | location == "primary_teratoma" | location == "metastasis_midgut" | location == "metastasis_rectum" | location == "metastasis_pancreas")
anno <- anno[c(1:9, 11:27, 29:38, 40:55, 58:60),] # remove NEC and duplicate cases

# obtain betas from selected cases
betas <- betas[, anno$arrayId]

# determine most variable probes across dataset and subset beta values
probe_var <- apply(betas, 1, var)
probes_topvar <- rownames(betas)[order(probe_var, decreasing = TRUE)]
#saveRDS(object = probes_topvar, file = "./output/top_variable_probes_20240121.rds")
probes_topvar <- readRDS("./output/top_variable_probes_20240121.rds")

# pick betas for 5,000 top variable probes
betas_topvar <- betas[probes_topvar[1:5000], ]


# heatmap of sample-wise correlations
sample_cor <- cor(betas_topvar, method = "pearson")

sample_cor <- cor(betas_topvar, method = "spearman")


rownames(sample_cor) <- anno$annotation
pheatmap::pheatmap(sample_cor, labels_row = anno$sampleName, show_colnames = FALSE)


# run UMAP
set.seed(45098)
umap_settings <- umap.defaults
umap_settings$n_neighbors = 15
umap_settings$min_dist = 0.1

umap<- umap(d = t(betas_topvar), config = umap_settings, ret_model = TRUE)

anno <- anno %>% 
  mutate(umap_x = umap$layout[, 1], 
         umap_y = umap$layout[, 2])

#saveRDS(object = umap, file = "./output/umap_model_20240121.rds")
#saveRDS(object = anno, file = "./output/sample_annotation_umap_purity_20240121.rds")

umap <- readRDS("./output/umap_model_20240121.rds")
anno <- readRDS("./output/sample_annotation_umap_purity_ovary_20240121.rds")


# plot UMAP
anno %>% 
  ggplot(aes(umap_x, umap_y, col = annotation)) +
  geom_point(size = 4) +
  #geom_text(aes(label = sampleName), size = 4) +
  theme_classic(base_size = 24) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

