# https://bioconductor.org/packages/release/bioc/vignettes/DMRcate/inst/doc/DMRcate.pdf

# libraries
library(DMRcate)
library(tidyverse)
library(minfi)

# load raw data
raw_epic <- readRDS("./input/raw_epic_ovary.rds")

detP <- detectionP(raw_epic)

remove <- apply(detP, 1, function (x) any(x > 0.01))
ovary <- preprocessFunnorm(raw_epic)
ovary_back_up <- ovary

ovary <- ovary[!rownames(ovary) %in% names(which(remove)),]


### select ovary primary and secondary tumors - ANOVA analysis ------------------------------------------------
anno <- readRDS("./input/sample_annotation.rds")

anno_ovary <- anno %>%
  filter(tumorType == "ovary")
anno_ovary <- anno_ovary[c(1:9, 11:18),] # remove NEC and two duplicate samples

anno_ovary <- anno_ovary %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' ~ 'primary_not_in_teratoma',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary'))

#anno_ovary <- anno_ovary %>%
#  filter(annotation == "primary_not_in_teratoma" | annotation == "metastasis_to_ovary")

ovary <- ovary[, t(anno_ovary$arrayId)]
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

# annotate matrix of M-values
type <- factor(anno_ovary$annotation)
design <- model.matrix(~type)
myannotation <- cpg.annotate("array", ovary, arraytype = "EPIC",
                             analysis.type="ANOVA", design=design, fdr = 0.05, coef=2)

# obtain DMRs
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)

# convert DMRs in GRanges object
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results.ranges

# plot results
groups <- c("primary_not_in_teratoma" ="magenta", "primary_teratoma" ="forestgreen", "metastasis_to_ovary" = "blue")
cols <- groups[as.character(type)]
cols

DMR.plot(ranges=results.ranges, dmr=1, CpGs=getBeta(ovary), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")




### Ovary primary tumors only (primary_teratoma vs primary_not_in_teratoma) - differential analysis ------------------------------------------------------------------------------------

anno <- readRDS("./input/sample_annotation.rds")
ovary <- ovary_back_up

anno_ovary <- anno %>%
  filter(tumorType == "ovary")
anno_ovary <- anno_ovary[c(1:9, 11:18),] # remove NEC and two duplicate samples

anno_ovary <- anno_ovary %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' ~ 'primary_not_in_teratoma',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary'))

anno_ovary <- anno_ovary %>%
  filter(annotation == "primary_not_in_teratoma" | annotation == "primary_teratoma")

ovary <- ovary[, t(anno_ovary$arrayId)]
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

# annotate matrix of M-values
type <- factor(anno_ovary$annotation)
design <- model.matrix(~type)
myannotation <- cpg.annotate("array", ovary, arraytype = "EPIC",
                             analysis.type="differential", design=design, fdr = 0.1, coef=2)
#Your contrast returned no individually significant probes. 
#Try increasing the fdr. 
# Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.

# obtain DMRs
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
#   The FDR you specified in cpg.annotate() returned no significant CpGs, hence there are no DMRs.
# Try specifying a value of 'pcutoff' in dmrcate() and/or increasing 'fdr' in cpg.annotate().

# convert DMRs in GRanges object
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results.ranges

# plot results
groups <- c("primary_not_in_teratoma" ="magenta", "primary_teratoma" ="forestgreen")
cols <- groups[as.character(type)]
cols

DMR.plot(ranges=results.ranges, dmr=1, CpGs=getBeta(ovary), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")



### Primary_not_in_teratoma vs metastasis_to_ovary - differential analysis ------------------------------------------------------------------------------------
anno <- readRDS("./input/sample_annotation.rds")
ovary <- ovary_back_up

anno_ovary <- anno %>%
  filter(tumorType == "ovary")
anno_ovary <- anno_ovary[c(1:9, 11:18),] # remove NEC and two duplicate samples

anno_ovary <- anno_ovary %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' ~ 'primary_not_in_teratoma',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary'))

anno_ovary <- anno_ovary %>%
  filter(annotation == "primary_not_in_teratoma" | annotation == "metastasis_to_ovary")

ovary <- ovary[, t(anno_ovary$arrayId)]
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

# annotate matrix of M-values
type <- factor(anno_ovary$annotation)
design <- model.matrix(~type)
myannotation <- cpg.annotate("array", ovary, arraytype = "EPIC",
                             analysis.type="differential", design=design, fdr = 0.05, coef=2)
# Your contrast returned 71 individually significant probes; a small but real effect. 
# Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors.


# obtain DMRs
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)

# convert DMRs in GRanges object
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results.ranges
#GRanges object with 3 ranges and 8 metadata columns:
#  seqnames            ranges strand   |   no.cpgs min_smoothed_fdr    Stouffer     HMFDR      Fisher   maxdiff
#<Rle>         <IRanges>  <Rle>        | <integer>        <numeric>   <numeric> <numeric>   <numeric> <numeric>
#[1]     chr7 94285520-94287242      * |        52     8.04638e-261 5.86776e-10 0.1126885 1.59944e-07  0.493884
#[2]    chr19 57351791-57352542      * |        13     1.29533e-125 7.69392e-05 0.0864563 4.77013e-04  0.506253
#[3]    chr11   2721409-2721632      * |         6     5.68071e-113 5.71960e-04 0.0678882 3.39952e-03  0.523635
#meandiff overlapping.genes
#<numeric>       <character>
#[1]  0.271228       PEG10, SGCE
#[2]  0.304917 MIMT1, ZIM2, PEG3
#[3]  0.355594             KCNQ1
#-
#  seqinfo: 3 sequences from an unspecified genome; no seqlengths


# plot results
groups <- c("primary_not_in_teratoma" ="magenta", "metastasis_to_ovary" = "blue")
cols <- groups[as.character(type)]
cols

DMR.plot(ranges=results.ranges, dmr=1, CpGs=getBeta(ovary), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")



### Primary_not_in_teratoma vs primary_ileal - differential analysis ------------------------------------------------------------------------------------
anno <- readRDS("./input/sample_annotation.rds")
ovary <- ovary_back_up

anno_ovary <- anno %>%
  filter(tumorType == "ovary" | tumorType == "ilealNET")
anno_ovary <- anno_ovary[c(1:9, 11:29),] # remove NEC and two duplicate samples

anno_ovary <- anno_ovary %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' & tumorType == "ovary" ~ 'primary_not_in_teratoma',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary',
                                location == 'primary' & tumorType == "ilealNET" ~ 'primary_ileal',
                                location == 'metastasis'& tumorType == "ilealNET" ~ 'ileal_metastasis'))

anno_ovary <- anno_ovary %>%
  filter(annotation == "primary_not_in_teratoma" | annotation == "primary_ileal")

ovary <- ovary[, t(anno_ovary$arrayId)]
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

# annotate matrix of M-values
type <- factor(anno_ovary$annotation)
design <- model.matrix(~type)
myannotation <- cpg.annotate("array", ovary, arraytype = "EPIC",
                             analysis.type="differential", design=design, fdr = 0.05, coef=2)
#Your contrast returned 825 individually significant probes. We recommend the default setting of pcutoff in dmrcate().


# obtain DMRs
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)

# convert DMRs in GRanges object
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results.ranges
#GRanges object with 67 ranges and 8 metadata columns:
#seqnames              ranges strand |   no.cpgs min_smoothed_fdr    Stouffer      HMFDR      Fisher   maxdiff
#<Rle>           <IRanges>  <Rle> | <integer>        <numeric>   <numeric>  <numeric>   <numeric> <numeric>
#[1]     chr7   94285004-94287242      * |        61     1.22779e-258 5.24441e-22  0.0330729 1.52855e-20  0.414708
#[2]    chr19   57349204-57353128      * |        34      4.84627e-90 3.79481e-15  0.0292736 6.80550e-12  0.467709
#[3]    chr14 101291440-101294623      * |        29      1.68140e-99 1.23438e-09  0.0165842 1.92533e-10 -0.514302
#[4]     chr7 130130288-130133110      * |        48     9.57710e-105 2.21912e-12  0.0921838 1.19214e-08  0.352924
#[5]    chr11     2018635-2021103      * |        38      8.64945e-94 8.50384e-11  0.0763011 1.22477e-08 -0.416437
#...      ...                 ...    ... .       ...              ...         ...        ...         ...       ...
#[63]     chr6     3849272-3850106      * |        21      8.99477e-27  0.00657053 0.15307736   0.0463996  0.262690
#[64]     chr4     1166528-1167423      * |        15      3.75070e-22  0.02435028 0.15749888   0.0478490 -0.495560
#[65]     chr6   32118204-32118457      * |        13      7.30269e-17  0.01125367 0.16395837   0.0528915 -0.311694
#[66]     chr6   30710755-30711863      * |        31      6.83032e-37  0.03128309 0.18408314   0.0884601 -0.333961
#[67]    chr22   20104664-20105375      * |        12      2.30158e-21  0.76636464 0.00151386   0.2831024  0.268566
#meandiff      overlapping.genes
#<numeric>            <character>
#  [1]   0.220852            PEG10, SGCE
#[2]   0.241347      MIMT1, ZIM2, PEG3
#[3]  -0.253550                   MEG3
#[4]   0.204660                   MEST
#[5]  -0.228153                    H19
#...        ...                    ...
#[63]  0.1365978   RP11-420L9.4, FAM50B
#[64] -0.1815520                  SPON2
#[65] -0.1632077                  PRRT1
#[66] -0.1207491 XXbac-BPG252P9.10, I..
#[67]  0.0217722         RANBP1, TRMT2A
#-------
#seqinfo: 20 sequences from an unspecified genome; no seqlengths

# plot results
groups <- c("primary_not_in_teratoma" ="magenta", "primary_ileal" = "blue")
cols <- groups[as.character(type)]
cols

DMR.plot(ranges=results.ranges, dmr=1, CpGs=getBeta(ovary), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")



### Primary_teratoma vs primary_rectal - differential analysis ------------------------------------------------------------------------------------
anno <- readRDS("./input/sample_annotation.rds")
ovary <- ovary_back_up

anno_ovary <- anno %>%
  filter(tumorType == "ovary" | tumorType == "rectalNET")
anno_ovary <- anno_ovary[c(1:9, 11:24, 27:30),] # remove NEC and two duplicate samples

anno_ovary <- anno_ovary %>%
  mutate(annotation = case_when(location == 'primary_teratoma' ~ 'primary_teratoma',
                                location == 'primary' & tumorType == "ovary" ~ 'primary_not_in_teratoma',
                                location == 'metastasis_midgut' ~ 'metastasis_to_ovary',
                                location == 'metastasis_rectum' ~ 'metastasis_to_ovary',
                                location == 'metastasis_pancreas' ~ 'metastasis_to_ovary',
                                location == 'primary' & tumorType == "rectalNET" ~ 'primary_rectal',
                                location == 'metastasis'& tumorType == "rectalNET" ~ 'rectal_metastasis'))

anno_ovary <- anno_ovary %>%
  filter(annotation == "primary_teratoma" | annotation == "primary_rectal")

ovary <- ovary[, t(anno_ovary$arrayId)]
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

# annotate matrix of M-values
type <- factor(anno_ovary$annotation)
design <- model.matrix(~type)
myannotation <- cpg.annotate("array", ovary, arraytype = "EPIC",
                             analysis.type="differential", design=design, fdr = 0.05, coef=2)

# obtain DMRs
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)

# convert DMRs in GRanges object
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results.ranges
#GRanges object with 2 ranges and 8 metadata columns:
#seqnames              ranges strand |   no.cpgs min_smoothed_fdr    Stouffer      HMFDR      Fisher   maxdiff
#<Rle>           <IRanges>  <Rle> | <integer>        <numeric>   <numeric>  <numeric>   <numeric> <numeric>
#[1]     chr7 156798258-156798931      * |         4     2.03581e-100 3.75621e-05 0.00056886 4.94708e-05 -0.687841
#[2]     chr7   94285873-94286669      * |        32      3.72944e-87 4.24450e-03 0.27720900 1.39937e-01  0.310704
#meandiff overlapping.genes
#<numeric>       <character>
#[1] -0.469348              MNX1
#[2]  0.227394             PEG10
#-------
#seqinfo: 1 sequence from an unspecified genome; no seqlengths


# plot results
groups <- c("primary_teratoma" ="forestgreen", "primary_rectal" = "blue")
cols <- groups[as.character(type)]
cols

DMR.plot(ranges=results.ranges, dmr=1, CpGs=getBeta(ovary), what="Beta",
         arraytype = "EPIC", phen.col=cols, genome="hg19")
