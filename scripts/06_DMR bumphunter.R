#https://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#bumphunter-to-find-differentially-methylated-regions-dmrs

# libraries
library(tidyverse)
library(minfi)
library(bumphunter)

options(scipen = 999)

# load anno data
anno <- readRDS("./input/sample_annotation.rds")

anno <- anno %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' & tumorType == "ovary" ~ 'primary_not_in_teratoma',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary',
                                location == 'primary' & tumorType == 'ilealNET' ~ 'primary_ileal',
                                location == 'metastasis'& tumorType == 'ilealNET' ~ 'ileal_metastasis',
                                location == 'primary' & tumorType == 'rectalNET' ~ 'primary_rectal',
                                location == 'metastasis'& tumorType == 'rectalNET' ~ 'rectal_metastasis',
                                location == 'primary' & tumorType == 'pulmNET' ~ 'primary_pulmonall',
                                location == 'primary' & tumorType == 'pulmNEC' ~ 'primary_pulmonall'))

# load unfiltered beta values
betas <- readRDS("./input/betas_everything.rds")
betas_back_up <- betas

### DMP finder ---------------------------------------------------------------------------------------------

# 1. primary teratoma vs primary not in teratoma
anno_ovary <- anno %>%
  filter(annotation == "primary_not_in_teratoma" | annotation == "primary_teratoma")
anno_ovary <- anno_ovary[c(1:9, 11),] # remove NEC and two duplicate samples
betas <- betas[, anno_ovary$arrayId]

# DMP finder
teratoma  <- anno_ovary$annotation
dmp <- dmpFinder(betas, pheno = teratoma  , type = "categorical")
head(dmp)

dmp0.05 <- dmp %>%
  filter(pval < 0.05)
#50508

# 2. primary not in teratoma vs ileal metastasis
anno_ovary <- anno %>%
  filter(annotation == "primary_not_in_teratoma" | annotation == "metastasis_to_ovary")
anno_ovary <- anno_ovary[c(1:5, 7:12),] # remove NEC and two duplicate samples
betas <- betas_back_up
betas <- betas[, anno_ovary$arrayId]

# DMP finder
teratoma  <- anno_ovary$annotation
dmp <- dmpFinder(betas, pheno = teratoma  , type = "categorical")
head(dmp)

dmp0.05 <- dmp %>%
  filter(pval < 0.05)
#31372


# 3. primary not in teratoma vs ileal primary
anno_ovary <- anno %>%
  filter(annotation == "primary_not_in_teratoma" | annotation == "primary_ileal")
anno_ovary <- anno_ovary[c(1:5, 7:16),] # remove NEC and two duplicate samples
betas <- betas_back_up
betas <- betas[, anno_ovary$arrayId]

# DMP finder
teratoma  <- anno_ovary$annotation
dmp <- dmpFinder(betas, pheno = teratoma  , type = "categorical")
head(dmp)

dmp0.05 <- dmp %>%
  filter(pval < 0.05)
#87773

# 4. metastasis_to_ovary vs ileal primary
anno_ovary <- anno %>%
  filter(annotation == "metastasis_to_ovary" | annotation == "primary_ileal")
betas <- betas_back_up
betas <- betas[, anno_ovary$arrayId]

# DMP finder
teratoma  <- anno_ovary$annotation
dmp <- dmpFinder(betas, pheno = teratoma  , type = "categorical")
head(dmp)

dmp0.05 <- dmp %>%
  filter(pval < 0.05)
#76702

# 5. primary_teratoma vs rectal primary
anno_ovary <- anno %>%
  filter(annotation == "primary_teratoma" | annotation == "primary_rectal")
anno_ovary <- anno_ovary[c(1:9, 11:14),] # remove NEC and two duplicate samples
betas <- betas_back_up
betas <- betas[, anno_ovary$arrayId]

# DMP finder
teratoma  <- anno_ovary$annotation
dmp <- dmpFinder(betas, pheno = teratoma  , type = "categorical")
head(dmp)

dmp0.05 <- dmp %>%
  filter(pval < 0.05)
#49162



### DMR finder----------------------------------------------------------------------------------------
# load raw data
raw_epic <- readRDS("./input/raw_epic_ovary.rds")

detP <- detectionP(raw_epic)

remove <- apply(detP, 1, function (x) any(x > 0.01))
ovary <- preprocessFunnorm(raw_epic)
ovary_back_up <- ovary




# 1. primary teratoma vs primary not in teratoma
anno_ovary <- anno %>%
  filter(annotation == "primary_not_in_teratoma" | annotation == "primary_teratoma")
anno_ovary <- anno_ovary[c(1:9, 11),] # remove NEC and two duplicate samples
ovary <- ovary[, anno_ovary$arrayId]

# define phenotype of interest
pheno <- anno_ovary$annotation
designMatrix <- model.matrix(~ pheno)

# set cutoff
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=0, type="Beta")

# include permutations
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=1000, type="Beta")

names(dmrs)
head(dmrs$table, n=3)
dim(dmrs$table)
#[1] 27 14 > 27 DMRs, waarvan 1 pval <0,05
dmr_df <- as.data.frame(dmrs$table)


# 2. primary not in teratoma vs ileal metastasis
anno_ovary <- anno %>%
  filter(annotation == "primary_not_in_teratoma" | annotation == "metastasis_to_ovary")
anno_ovary <- anno_ovary[c(1:5, 7:12),] # remove NEC and two duplicate samples
ovary <- ovary_back_up
ovary <- ovary[, anno_ovary$arrayId]

# define phenotype of interest
pheno <- anno_ovary$annotation
designMatrix <- model.matrix(~ pheno)

# set cutoff
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=0, type="Beta")

# include permutations
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=1000, type="Beta")

names(dmrs)
head(dmrs$table, n=3)
dim(dmrs$table)
# 48, waarvan 5 p val <0.05
dmr_df <- as.data.frame(dmrs$table)

# 3. primary not in teratoma vs ileal primary
anno_ovary <- anno %>%
  filter(annotation == "primary_not_in_teratoma" | annotation == "primary_ileal")
anno_ovary <- anno_ovary[c(1:5, 7:16),] # remove NEC and two duplicate samples
ovary <- ovary_back_up
ovary <- ovary[, anno_ovary$arrayId]

# define phenotype of interest
pheno <- anno_ovary$annotation
designMatrix <- model.matrix(~ pheno)

# set cutoff
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=0, type="Beta")

# include permutations
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=1000, type="Beta")

names(dmrs)
head(dmrs$table, n=3)
dim(dmrs$table)
# 33 waarvan 4 p val <0.05
dmr_df <- as.data.frame(dmrs$table)


# 4. metastasis_to_ovary vs ileal primary
anno_ovary <- anno %>%
  filter(annotation == "metastasis_to_ovary" | annotation == "primary_ileal")
ovary <- ovary_back_up
ovary <- ovary[, anno_ovary$arrayId]

# define phenotype of interest
pheno <- anno_ovary$annotation
designMatrix <- model.matrix(~ pheno)

# set cutoff
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=0, type="Beta")

# include permutations
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=1000, type="Beta")

names(dmrs)
head(dmrs$table, n=3)
dim(dmrs$table)
# 31 waarvan 6 p val <0.05
dmr_df <- as.data.frame(dmrs$table)


# 5. primary_teratoma vs rectal primary
anno_ovary <- anno %>%
  filter(annotation == "primary_teratoma" | annotation == "primary_rectal")
anno_ovary <- anno_ovary[c(1:9, 11:14),] # remove NEC and two duplicate samples
ovary <- ovary_back_up
ovary <- ovary[, anno_ovary$arrayId]

# define phenotype of interest
pheno <- anno_ovary$annotation
designMatrix <- model.matrix(~ pheno)

# set cutoff
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=0, type="Beta")

# include permutations
dmrs <- bumphunter(ovary, design = designMatrix, 
                   cutoff = 0.5, B=1000, type="Beta")

names(dmrs)
head(dmrs$table, n=3)
dim(dmrs$table)
# 56 waarvan 6 p val <0.05
dmr_df <- as.data.frame(dmrs$table)
