library("sesame")

data <- sesameDataGet('HM450.76.TCGA.matched')
cf <- DMR(betas_ovary, anno_ovary$sampleName, ~type)


# libraries

library(minfi)
library(tidyverse)
library(umap)
library(clusteval)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

source("./scripts/0_helpers.R")

### OVARY --------------------------------------------------------------------------
# import sample annotation
sample_anno <- readRDS(file = "./output/sample_annotation_umap_purity_ovary_20240105.rds")
sample_anno <- sample_anno[1:16,]
  
# import betas 
betas <- readRDS("./input/betas_filtered.rds")
betas <- betas[, sample_anno$arrayId]

# import most variable probes
topvar_probes <- readRDS("./output/top_variable_probes_ovary_20240105.rds")
topvar_probes <- topvar_probes[1:5000]
betas_topvar <- betas[topvar_probes, ]
#betas_topvar <- betas[sample(1:nrow(betas), size = 5000), ]


# import EPIC annotation
anno_array <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno_array <- anno_array[rownames(betas), ] %>% 
  as_tibble(rownames = "probe_id")

# clean up the array annotation a little bit
anno_array <- anno_array %>% 
  mutate(topvar = ifelse(probe_id %in% topvar_probes, "topvar", "other"), 
         Relation_to_Island = str_replace(Relation_to_Island, 
                                          pattern = "[N,S]{1}_", replacement = ""))


# investigate distribution of most variable probes -----------------------------

# extract probe locations
topvar_location <- anno_array %>% 
  select(topvar, Relation_to_Island) %>% 
  group_by(topvar, Relation_to_Island) %>% 
  summarise(n = n()) %>%
  group_by(topvar) %>% 
  mutate(prop = n / sum(n) *100) %>% 
  ungroup
#saveRDS(object = topvar_location, file = "./output/topvar_locations_20240105.rds")

# plot
topvar_location %>% 
  ggplot(aes(Relation_to_Island, prop, fill = topvar)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 24) +
  labs(x = NULL, y = "Proportion of probes (%)")

# Chi-Square test (p-val = 4.281e-08)
topvar_location %>% 
  pivot_wider(id_cols = Relation_to_Island, names_from = topvar, values_from = n) %>% 
  select(other, topvar) %>% 
  chisq.test()

# what about gene bodies ? -----------------------------------------------------

# extract probes which are unequivocally associated with one gene
gb <- anno_array %>% 
  filter(UCSC_RefGene_Accession != "") %>% 
  filter(str_detect(UCSC_RefGene_Accession, pattern = ";", negate = TRUE))

# count and normalize w.r.t. to gene position
gb_counts <- gb %>% 
  mutate(cpg_island = ifelse(Relation_to_Island == "Island", "CpG Island", "Non-Island")) %>% 
  select(UCSC_RefGene_Group, cpg_island, topvar) %>% 
  dplyr::rename("location" = UCSC_RefGene_Group) %>% 
  group_by(location, cpg_island, topvar) %>% 
  summarise(n = n()) %>% 
  group_by(topvar) %>% 
  mutate(prop = n/sum(n) * 100)

gb_counts %>% 
  ggplot(aes(location, prop, fill = topvar)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 24) +
  facet_wrap(nrow = 2, facets = vars(cpg_island), scales = "free_y") +
  labs(x = NULL, y = "Proportion of probes (%)")


# look at differential methylation ---------------------------------------------

prom_anno <- gb %>% 
  filter(Relation_to_Island == "Island", 
         str_detect(UCSC_RefGene_Group, pattern = "TSS")) %>% 
  group_by(UCSC_RefGene_Accession) %>% 
  mutate(n = n()) %>% 
  filter(n > 3)

betas_prom <- betas[prom_anno$probe_id, sample_anno$arrayId]
betas_prom_avg <- apply(betas_prom, 2, FUN = function(x) 
  tapply(x, INDEX = prom_anno$UCSC_RefGene_Name, FUN = mean, na.rm = TRUE))

prom_umap <- umap(t(betas_prom_avg))

sample_anno <- sample_anno %>% 
  mutate(umap_prom_x = prom_umap$layout[, 1], 
         umap_prom_y = prom_umap$layout[, 2])

sample_anno %>% 
  ggplot(aes(umap_prom_x, umap_prom_y, col = annotation)) +
  geom_point(size = 3) +
  theme_minimal(base_size = 18) +
  labs(x = "Umap 1", y = "Umap 2")

diff_meth_groups <- c("primary_teratoma", "primary_not_in_teratoma")

anno_comp <- sample_anno %>% 
  filter(annotation %in% diff_meth_groups)

diff_meth <- genefilter::rowttests(betas_prom_avg[, anno_comp$arrayId], 
                                   fac = as.factor(anno_comp$annotation))
diff_meth$mean_pnet = apply(betas_prom_avg[, sample_anno$annotation == diff_meth_groups[1]], 1, mean)
diff_meth$mean_norm = apply(betas_prom_avg[, sample_anno$annotation == diff_meth_groups[2]], 1, mean)

diff_meth <- diff_meth %>% 
  as_tibble(rownames = "gene_symbol") %>% 
  mutate(pval_log = -log10(p.value))

sum(diff_meth$p.value < 0.001)

diff_meth %>% 
  mutate(filtered_label = ifelse(pval_log > 10, gene_symbol, "")) %>% 
  ggplot(aes(dm, pval_log)) +
  geom_point() +
  geom_text(aes(label = filtered_label), nudge_y = 0.5) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

diff_meth %>% 
  ggplot(aes(dm)) +
  geom_histogram(bins = 200) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

write_csv(x = diff_meth, file = "./output/diff_meth_teratoma_vs_non_teratoma.csv")

### metastasis to ovary vs primary not in teratoma
diff_meth_groups <- c("metastasis_to_ovary", "primary_not_in_teratoma")

anno_comp <- sample_anno %>% 
  filter(annotation %in% diff_meth_groups)

diff_meth <- genefilter::rowttests(betas_prom_avg[, anno_comp$arrayId], 
                                   fac = as.factor(anno_comp$annotation))
diff_meth$mean_pnet = apply(betas_prom_avg[, sample_anno$annotation == diff_meth_groups[1]], 1, mean)
diff_meth$mean_norm = apply(betas_prom_avg[, sample_anno$annotation == diff_meth_groups[2]], 1, mean)

diff_meth <- diff_meth %>% 
  as_tibble(rownames = "gene_symbol") %>% 
  mutate(pval_log = -log10(p.value))

sum(diff_meth$p.value < 0.001)

diff_meth %>% 
  mutate(filtered_label = ifelse(pval_log > 10, gene_symbol, "")) %>% 
  ggplot(aes(dm, pval_log)) +
  geom_point() +
  geom_text(aes(label = filtered_label), nudge_y = 0.5) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

diff_meth %>% 
  ggplot(aes(dm)) +
  geom_histogram(bins = 200) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

write_csv(x = diff_meth, file = "./output/diff_meth_metastasis_vs_non_teratoma.csv")

### OVARY_NOT_IN_TERATOMA + ILEAL -------------------------------------------------------------------

# import sample annotation
sample_anno <- readRDS(file = "./output/sample_annotation_umap_purity_ovary_ileal_20240105.rds")
sample_anno <- sample_anno %>%
  filter(type == 'primary')

# import betas 
betas <- readRDS("./input/betas_filtered.rds")
betas <- betas[, sample_anno$arrayId]

# import most variable probes
topvar_probes <- readRDS("./output/top_variable_probes_ovary_ileal_200240105.rds")
topvar_probes <- topvar_probes[1:5000]
betas_topvar <- betas[topvar_probes, ]
#betas_topvar <- betas[sample(1:nrow(betas), size = 5000), ]


# import EPIC annotation
anno_array <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno_array <- anno_array[rownames(betas), ] %>% 
  as_tibble(rownames = "probe_id")

# clean up the array annotation a little bit
anno_array <- anno_array %>% 
  mutate(topvar = ifelse(probe_id %in% topvar_probes, "topvar", "other"), 
         Relation_to_Island = str_replace(Relation_to_Island, 
                                          pattern = "[N,S]{1}_", replacement = ""))


# investigate distribution of most variable probes -----------------------------

# extract probe locations
topvar_location <- anno_array %>% 
  select(topvar, Relation_to_Island) %>% 
  group_by(topvar, Relation_to_Island) %>% 
  summarise(n = n()) %>%
  group_by(topvar) %>% 
  mutate(prop = n / sum(n) *100) %>% 
  ungroup
#saveRDS(object = topvar_location, file = "./output/topvar_locations_20240105.rds")

# plot
topvar_location %>% 
  ggplot(aes(Relation_to_Island, prop, fill = topvar)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 24) +
  labs(x = NULL, y = "Proportion of probes (%)")

# Chi-Square test (p-val = 4.281e-08)
topvar_location %>% 
  pivot_wider(id_cols = Relation_to_Island, names_from = topvar, values_from = n) %>% 
  select(other, topvar) %>% 
  chisq.test()

# what about gene bodies ? -----------------------------------------------------

# extract probes which are unequivocally associated with one gene
gb <- anno_array %>% 
  filter(UCSC_RefGene_Accession != "") %>% 
  filter(str_detect(UCSC_RefGene_Accession, pattern = ";", negate = TRUE))

# count and normalize w.r.t. to gene position
gb_counts <- gb %>% 
  mutate(cpg_island = ifelse(Relation_to_Island == "Island", "CpG Island", "Non-Island")) %>% 
  select(UCSC_RefGene_Group, cpg_island, topvar) %>% 
  dplyr::rename("location" = UCSC_RefGene_Group) %>% 
  group_by(location, cpg_island, topvar) %>% 
  summarise(n = n()) %>% 
  group_by(topvar) %>% 
  mutate(prop = n/sum(n) * 100)


# look at differential methylation ---------------------------------------------

prom_anno <- gb %>% 
  filter(Relation_to_Island == "Island", 
         str_detect(UCSC_RefGene_Group, pattern = "TSS")) %>% 
  group_by(UCSC_RefGene_Accession) %>% 
  mutate(n = n()) %>% 
  filter(n > 3)

betas_prom <- betas[prom_anno$probe_id, sample_anno$arrayId]
betas_prom_avg <- apply(betas_prom, 2, FUN = function(x) 
  tapply(x, INDEX = prom_anno$UCSC_RefGene_Name, FUN = mean, na.rm = TRUE))

prom_umap <- umap(t(betas_prom_avg))

sample_anno <- sample_anno %>% 
  mutate(umap_prom_x = prom_umap$layout[, 1], 
         umap_prom_y = prom_umap$layout[, 2])


### primary ovary vs primary ileal
diff_meth_groups <- c("ovary", "ileal")

anno_comp <- sample_anno %>% 
  filter(localization %in% diff_meth_groups)

diff_meth <- genefilter::rowttests(betas_prom_avg[, anno_comp$arrayId], 
                                   fac = as.factor(anno_comp$localization))
diff_meth$mean_pnet = apply(betas_prom_avg[, sample_anno$localization == diff_meth_groups[1]], 1, mean)
diff_meth$mean_norm = apply(betas_prom_avg[, sample_anno$localization == diff_meth_groups[2]], 1, mean)

diff_meth <- diff_meth %>% 
  as_tibble(rownames = "gene_symbol") %>% 
  mutate(pval_log = -log10(p.value))

sum(diff_meth$p.value < 0.001)
# [1] 7

diff_meth %>% 
  mutate(filtered_label = ifelse(pval_log > 10, gene_symbol, "")) %>% 
  ggplot(aes(dm, pval_log)) +
  geom_point() +
  geom_text(aes(label = filtered_label), nudge_y = 0.5) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

diff_meth %>% 
  ggplot(aes(dm)) +
  geom_histogram(bins = 200) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

write_csv(x = diff_meth, file = "./output/diff_primary_ovary_vs_primary_ilel.csv")


### OVARY_IN_TERATOMA + RECTAL -------------------------------------------------------------------

#anno_ovary_rectal <- anno %>%
#  filter(tumorType == "ovary" & location == "primary_teratoma" | tumorType == "rectalNET")
#anno_ovary_rectal <- anno_ovary_rectal[c(1:11, 13:16),]

# determine most variable probes across dataset and subset beta values
#betas <- readRDS("./input/betas_filtered.rds")
#betas_ovary_rectal <- betas[, anno_ovary_rectal$arrayId]
#probe_var_ovary_rectal <- apply(betas_ovary_rectal, 1, var)
#probes_topvar_ovary_rectal <- rownames(betas_ovary_rectal)[order(probe_var_ovary_rectal, decreasing = TRUE)]
#saveRDS(object = probes_topvar_ovary_rectal, file = "./output/top_variable_probes_ovary_rectal_20240105.rds")


# import sample annotation
sample_anno <- anno_ovary_rectal

# import betas 
betas <- readRDS("./input/betas_filtered.rds")
betas <- betas[, sample_anno$arrayId]

# import most variable probes
topvar_probes <- readRDS("./output/top_variable_probes_ovary_rectal_20240105.rds")
topvar_probes <- topvar_probes[1:5000]
betas_topvar <- betas[topvar_probes, ]
#betas_topvar <- betas[sample(1:nrow(betas), size = 5000), ]


# import EPIC annotation
anno_array <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno_array <- anno_array[rownames(betas), ] %>% 
  as_tibble(rownames = "probe_id")

# clean up the array annotation a little bit
anno_array <- anno_array %>% 
  mutate(topvar = ifelse(probe_id %in% topvar_probes, "topvar", "other"), 
         Relation_to_Island = str_replace(Relation_to_Island, 
                                          pattern = "[N,S]{1}_", replacement = ""))


# investigate distribution of most variable probes -----------------------------

# extract probe locations
topvar_location <- anno_array %>% 
  select(topvar, Relation_to_Island) %>% 
  group_by(topvar, Relation_to_Island) %>% 
  summarise(n = n()) %>%
  group_by(topvar) %>% 
  mutate(prop = n / sum(n) *100) %>% 
  ungroup
#saveRDS(object = topvar_location, file = "./output/topvar_locations_20240105.rds")

# plot
topvar_location %>% 
  ggplot(aes(Relation_to_Island, prop, fill = topvar)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 24) +
  labs(x = NULL, y = "Proportion of probes (%)")

# Chi-Square test (p-val = 4.281e-08)
topvar_location %>% 
  pivot_wider(id_cols = Relation_to_Island, names_from = topvar, values_from = n) %>% 
  select(other, topvar) %>% 
  chisq.test()

# what about gene bodies ? -----------------------------------------------------

# extract probes which are unequivocally associated with one gene
gb <- anno_array %>% 
  filter(UCSC_RefGene_Accession != "") %>% 
  filter(str_detect(UCSC_RefGene_Accession, pattern = ";", negate = TRUE))

# count and normalize w.r.t. to gene position
gb_counts <- gb %>% 
  mutate(cpg_island = ifelse(Relation_to_Island == "Island", "CpG Island", "Non-Island")) %>% 
  select(UCSC_RefGene_Group, cpg_island, topvar) %>% 
  dplyr::rename("location" = UCSC_RefGene_Group) %>% 
  group_by(location, cpg_island, topvar) %>% 
  summarise(n = n()) %>% 
  group_by(topvar) %>% 
  mutate(prop = n/sum(n) * 100)


# look at differential methylation ---------------------------------------------

prom_anno <- gb %>% 
  filter(Relation_to_Island == "Island", 
         str_detect(UCSC_RefGene_Group, pattern = "TSS")) %>% 
  group_by(UCSC_RefGene_Accession) %>% 
  mutate(n = n()) %>% 
  filter(n > 3)

betas_prom <- betas[prom_anno$probe_id, sample_anno$arrayId]
betas_prom_avg <- apply(betas_prom, 2, FUN = function(x) 
  tapply(x, INDEX = prom_anno$UCSC_RefGene_Name, FUN = mean, na.rm = TRUE))

prom_umap <- umap(t(betas_prom_avg))

sample_anno <- sample_anno %>% 
  mutate(umap_prom_x = prom_umap$layout[, 1], 
         umap_prom_y = prom_umap$layout[, 2])


### ovary in teratoma vs rectal net
diff_meth_groups <- c("ovary", "rectalNET")

anno_comp <- sample_anno %>% 
  filter(tumorType %in% diff_meth_groups)

diff_meth <- genefilter::rowttests(betas_prom_avg[, anno_comp$arrayId], 
                                   fac = as.factor(anno_comp$tumorType))
diff_meth$mean_pnet = apply(betas_prom_avg[, sample_anno$tumorType == diff_meth_groups[1]], 1, mean)
diff_meth$mean_norm = apply(betas_prom_avg[, sample_anno$tumorType == diff_meth_groups[2]], 1, mean)

diff_meth <- diff_meth %>% 
  as_tibble(rownames = "gene_symbol") %>% 
  mutate(pval_log = -log10(p.value))

sum(diff_meth$p.value < 0.001)
# [1] 3

diff_meth %>% 
  mutate(filtered_label = ifelse(pval_log > 10, gene_symbol, "")) %>% 
  ggplot(aes(dm, pval_log)) +
  geom_point() +
  geom_text(aes(label = filtered_label), nudge_y = 0.5) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

diff_meth %>% 
  ggplot(aes(dm)) +
  geom_histogram(bins = 200) +
  theme_bw(base_size = 18) +
  labs(x = "Methylation change (diff. beta)", y = "Significance (-log10 p-value)")

write_csv(x = diff_meth, file = "./output/diff_primary_ovary_in_teratoma_vs_primary_rectal.csv")
