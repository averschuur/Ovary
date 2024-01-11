
## script to generate CNV plots based on EPIC data ##
parameters <- list()
## folder with the idat files of the samples that need to be analyzed ##
parameters$dataFolder <- "L:/pathologie/PRL/Groep-Brosens/2. Anna Vera/15. ovaryNET/Ovary/data"
## folder with a set of healthy control samples to use for normalisation ##
## control sampels were downloaded here: https://www.ncbi.nlm.nih.gov/gds/?term=GSE110530[ACCN]%20AND%20gsm[ETYP]
## this is a set of 12 male blood samples, this also means that female samples that are processed will
## seem to have a gain of chrX and a loss of chrY since they will always be compared to male samples 
parameters$controlFolder <- "L:/pathologie/PRL/Groep-Brosens/2. Anna Vera/15. ovaryNET/Ovary/data/ControlMetSet rds"
#parameters$controlFolder <- "L:/pathologie/Moleculair/Diagnostiek/Methylatie array/UMC_arrayTools/data/controlUMC/"

## required packages and installation ##
packages <- c("shiny","BiocManager")
bioPackages <- c("conumee","minfi","minfiData","minfiDataEPIC","IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

if (length(setdiff(packages, rownames(installed.packages()))) > 0){
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(bioPackages, rownames(installed.packages()))) > 0){
  BiocManager::install(setdiff(bioPackages, rownames(installed.packages())))  
}

library("minfi") 
library("conumee") 
library("minfiData") 
library('minfiDataEPIC')

## function to load folders ##
loadEPICFolders <- function(){
  inputChoicesEPIC <- list.dirs(parameters$dataFolder,recursive = F,full.names = F)
}

## function that returns a list of unique samples based on the idat files present in a folder ##
loadEPICsamples <- function(folder){
  files <- list.files(folder,recursive = F,full.names=F,pattern = ".idat")
  files <- sub("_Grn.idat","",files)
  files <- sub("_Red.idat","",files)
  return(unique(files))
}

## function that loads samples en performs preprocessing ##
loadAndProcessSamples <- function(dataFolder,sampleName){
  patient<-read.metharray(file.path(dataFolder,paste0(sampleName,"_Grn.idat")))
  gmSet = mapToGenome(preprocessRaw(patient))
  patient <- preprocessIllumina(patient)
  return(patient)
}

## function that loads and preprocesses the control samples, and also defines detailed regions that are highlighted in the CNV plot ##
addCustomDetails <- T
loadAndProcessControls <- function(controlFolder = parameters$controlFolder,addCustomDetails = F){
  control <- read.metharray.exp(controlFolder)
  gmSet <- mapToGenome(preprocessRaw(control))
  control <- preprocessIllumina(control)
  data("exclude_regions")
  data("detail_regions")
  ## if require custom detail regions can be added here, these should be ranges ##
if(addCustomDetails){
  detailMat <- cbind(as.matrix(detail_regions@ranges), as.matrix(detail_regions@elementMetadata))
  #    detailMat <- detailMat[-6,]
  #    EGFRrange <- cbind(c(55209000,55215500),rep(6500,2))
  #    EGFRrange <- rbind(EGFRrange,c(55224000,55275032-55224000))
  #    EGFRrange <- rbind(c(55084726,55200000-55084726),EGFRrange)
  #    rownames(EGFRrange) <- paste0("EGFR-",c("I","II-IV","V-VII","VII-XXVIII"))
  #    PLAG1range <- cbind(c(57060000,57070000,57080000,57090000,57100000,57110000,57120000,57130000,57140000),rep(10000,9))
  #    rownames(PLAG1range) <- paste("PLAG1-",c(1:9),sep="")
  #    MYCNrange <- t(as.matrix(c(16082000,6000)))
  #    rownames(MYCNrange) <- "MYCN-1"
  CDK4_range <- as.matrix(c(58140000, 10000))
  rownames(CDK4_range) <- "CDK4"
  PAX_FOXO_range <- rbind(c(222199887, 222298998 - 222199887),c(18630846,18748866-18630846),c(40555667,40666641-4055566))
  rownames(PAX_FOXO_range) <- c("PAX3","PAX7","FOXO1")
  
  #    EGFRrange <- rbind(EGFRrange,PLAG1range,MYCNrange, PAX_FOXO_range)
  EGFRrange <- rbind(PAX_FOXO_range,CDK4_range)
  rr <- IRanges(start = c(unlist(detailMat[,1]),EGFRrange[,1]),width = c(unlist(detailMat[,2]),EGFRrange[,2]),names = c(unlist(detailMat[,3]),rownames(EGFRrange)))
  #    seqNames <- c(as.character(detail_regions@seqnames[-6]),paste0("chr",c(rep("7",4),rep("8",9),rep("2", 2),"1","13")))          
  seqNames <- c(as.character(detail_regions@seqnames[-6]),paste0("chr",c("1","13","12")))          
  detail_regions <- GRanges(seqnames = seqNames,ranges = rr,name=c(unlist(detailMat[,3]),rownames(EGFRrange)))
}
  anno <- CNV.create_anno(array_type = "EPIC", exclude_regions = exclude_regions,detail_regions = detail_regions, chrXY = TRUE)
  refData <- list("control" = control, "anno" = anno)
  return(refData)
}

## function that segments the data and generates a CNV plot and an IGV browsable file ##
segmentData <- function(patient,refData,dataFolder){
  minfi.data <- CNV.load(patient) 
  minfi.controls <- CNV.load(refData$control)
  patient <- mapToGenome(patient)
  refData$anno@probes <- subsetByOverlaps(refData$anno@probes, granges(patient))
  x <- CNV.fit(minfi.data, minfi.controls, refData$anno) 
  x <- CNV.bin(x) 
  x <- CNV.detail(x) 
  x <- CNV.segment(x,alpha=0.005)
  igvData <- as.data.frame(cbind(as.character(refData$anno@bins@seqnames),refData$anno@bins@ranges@start,refData$anno@bins@ranges@width+refData$anno@bins@ranges@start,names(x@bin$ratio),x@bin$ratio))
  headerLine <- paste0("#track type=bedGraph name=",sampleNames(patient)," description=center_label visibility=full autoScale=off graphType=points viewLimits=-1:1 windowingFunction=none smoothingWindow=off")
  file <- paste0(dataFolder,"\\",sampleNames(patient),".igv")
  cat(headerLine, '\n',  file = file)
  write.table(igvData, file, append = T, col.names = F,row.names = F,quote = F,sep="\t")
  igvDetail <- as.data.frame(cbind(as.character(refData$anno@detail@seqnames),refData$anno@detail@ranges@start,refData$anno@detail@ranges@width+refData$anno@detail@ranges@start,names(x@detail$ratio),x@detail$ratio))
  file <- paste0(dataFolder,"\\",sampleNames(patient),".detail.igv")
  cat(headerLine, '\n',  file = file)
  write.table(igvDetail, file, append = T, col.names = F,row.names = F,quote = F,sep="\t")
  
  return(x)
}


