### CNV analysis ####

library(openxlsx)
library("writexl")
library(tidyr)
library(naniar)

## open data -------------------------------------------------------------------------
CNV <- read.xlsx("./input/CNV analysis.xlsx")

Chr18 <- CNV %>%
  select(arrayId, Chr18)

anno_ovary <- merge(anno_ovary, Chr18, by = "arrayId")
anno_ovary_ileal_pancreatic_rectal <- merge(anno_ovary_ileal_pancreatic_rectal, Chr18, by = "arrayId")

