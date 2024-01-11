### load required packages -----------------------------------------------------

library(tidyverse)
library(minfi)
library(doParallel)



# load annotation --------------------------------------------------------------

anno_files <- list.files(path = "./annotation/",
                         pattern = ".csv",
                         full.names = TRUE)


anno <- lapply(as.list(anno_files), read_csv)
anno <- Reduce(f = bind_rows, x = anno)


### prepare methylation data -----------------------------------------------------


idats <- list.files(path = "T:/pathologie/PRL/Groep-Brosens/METHYLATIE DATA/UMCU/Idat Files/",
                    recursive = TRUE,
                    full.names = TRUE, 
                    pattern = "_Grn.idat")

# infer platform from filesize
file_size <- file.info(idats)$size

idats_epic <- idats[file_size > 10000000]
idats_450k <- idats[file_size < 10000000]

raw_epic <- read.metharray(basenames = idats_epic, force = TRUE)
raw_epic_ovary <- raw_epic[, t(anno$arrayId)]
#raw_450k <- read.metharray(basenames = idats_450k)
saveRDS(object = raw_epic, file = "./input/raw_epic.rds")
saveRDS(object = raw_epic_ovary, file = "./input/raw_epic_ovary.rds")

# preprocess data
preprocessed_epic <- raw_epic %>%
  preprocessNoob(dyeMethod= "single")


# combine 450k and EPIC data, save full (combined) data to disk
betas_epic <- getBeta(preprocessed_epic)


# clean up
rm(raw_epic)
rm(idats, idats_epic, idats_450k, file_size)

# identify probes (1) on X/Y chr, (2) SNPs/multimappers or (3) are crossreactive

# sex chromosomes
platform_anno <- getAnnotation(preprocessed_epic)
filtered_probes <- list()
filtered_probes$xy <- which(platform_anno$chr %in% c("chrX","chrY"), arr.ind = TRUE)
filtered_probes$xy <- rownames(preprocessed_epic)[filtered_probes$xy]

# SNPs and multimappers
betas_epic_filtered <- preprocessed_epic %>% 
  mapToGenome %>% 
  dropLociWithSnps %>% 
  getBeta()
filtered_probes$snp <- intersect(rownames(betas_epic), rownames(betas_epic_filtered))
filtered_probes$snp <- which(rownames(betas_epic) %in% filtered_probes$snp, arr.ind = TRUE)
filtered_probes$snp <- rownames(betas_epic)[-filtered_probes$snp]

# cross reactive
filtered_probes$xr <- read.csv("https://github.com/sirselim/illumina450k_filtering/blob/master/48639-non-specific-probes-Illumina450k.csv",  
                               )
filtered_probes$xr <- filtered_probes$xr %>% as.character()

# clean up and save to file
saveRDS(object = filtered_probes, file = "./input/filtered_probes_list.rds")
rm(platform_anno)

# filter data 
probes_remove <- Reduce(f = union, x = filtered_probes)
probes_remove <- intersect(probes_remove, rownames(betas_epic_filtered))
probes_remove <- match(probes_remove, rownames(betas_epic_filtered))

betas_epic_filtered <- betas_epic[-probes_remove, ]


## double check annotation vs. beta values, save files to disk -----------------

dim(betas_epic)
dim(betas_epic_filtered)
colnames(betas_epic)

all(colnames(betas_epic) == colnames(betas_epic_filtered))
all(colnames(betas_epic) %in% anno$arrayId)
all(anno$arrayId %in% colnames(betas_epic))

## save data -------------------------------------------------------------------

# anno
saveRDS(object = anno, file = "./input/sample_annotation.rds")

# all betas
betas_epic <- betas_epic[, anno$arrayId]
saveRDS(object = betas_epic, file = "./input/betas_everything.rds")

# filtered betas
betas_epic_filtered <- betas_epic_filtered[, anno$arrayId]
saveRDS(object = betas_epic_filtered, 
        file = "./input/betas_filtered.rds")

# clean up
rm(betas_epic, betas_epic_filtered, filtered_probes, probes_remove)
rm(preprocessed_epic)
