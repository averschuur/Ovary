# https://bioconductor.org/packages/release/bioc/vignettes/DMRcate/inst/doc/DMRcate.pdf

# libraries
library(DMRcate)
library(tidyverse)
library(minfi)


raw_epic <- readRDS("./input/raw_epic.rds")

detP <- detectionP(raw_epic)

remove <- apply(detP, 1, function (x) any(x > 0.01))
ovary <- preprocessFunnorm(raw_epic)

ovary_back_up <- ovary


ovary <- ovary[!rownames(ovary) %in% names(which(remove)),]


# select ovary dataset
anno <- readRDS("./input/sample_annotation.rds")
ovary <- ovary_back_up
ovary <- ovary[, t(anno$arrayId)]

ovaryms <- getM(ovary)
nrow(ovaryms)

ovaryms.noSNPs <- rmSNPandCH(ovaryms, dist=2, mafcut=0.05)
nrow(ovaryms.noSNPs)

dim(ovary)
dim(ovaryms)
dim(ovaryms.noSNPs)


ovary <- ovary[rownames(ovaryms.noSNPs),]
colnames(ovaryms.noSNPs) <- colnames(ovary)
assays(ovary)[["M"]] <- ovaryms.noSNPs
assays(ovary)[["Beta"]] <- ilogit2(ovaryms.noSNPs)

### select ovary primary tumors only ------------------------------------------------
anno_ovary <- anno %>%
  filter(tumorType == "ovary")
anno_ovary <- anno_ovary[c(1:9, 11:18),]

anno_ovary <- anno_ovary %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' ~ 'primary_not_in_teratoma',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary'))

#anno_ovary <- anno_ovary %>%
#  filter(annotation == "primary_not_in_teratoma" | annotation == "metastasis_to_ovary")

ovary <- ovary[, t(anno_ovary$arrayId)]
ovaryms <- ovaryms[, t(anno_ovary$arrayId)]
ovaryms.noSNPs <- ovaryms.noSNPs[, t(anno_ovary$arrayId)]

type <- factor(anno_ovary$annotation)
design <- model.matrix(~type)
myannotation <- cpg.annotate("array", ovary, arraytype = "EPIC",
                             analysis.type="ANOVA", design=design, fdr = 0.05, coef=2)

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)

results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results.ranges

groups <- c("primary_not_in_teratoma" ="magenta", "primary_teratoma" ="forestgreen", "metastasis_to_ovary" = "blue")
cols <- groups[as.character(type)]
cols

DMR.plot(ranges=results.ranges, dmr=99886, CpGs=getBeta(ovary), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")


