#Libraries to Load:
library(data.table)
library(lme4)
library(stringr)
library(boot)
library(MASS)
library(Matrix)
library(parallel) #core parallelization.
library(methods) #S4 class OOP methods.


#Get functions and classes:
source("SCENTfunctions.R")


#Initialize directories: (Example)
input_atac <- "./RData/Data/pbmc_multimodal.atac.rds"
input_mrna <- "./RData/Data/pbmc_multimodal.rna.rds"
input_meta <- "./RData/Data/pbmc_multimodal.meta.rds"
input_gene_peak <- "./RData/Data/qced_Tnk.G2P.txt"
output <- "./Testing/Output/test_output.txt"

options(stringsAsFactors = F)

#Read-in Necessary Files:
atac <- readRDS(input_atac)
mrna <- readRDS(input_mrna)
meta <- readRDS(input_meta)
gene_peak <- read.table(input_gene_peak)


#####1st way######
library(data.table)
library(Matrix)
library(stringr)
library(Hmisc)
library(R.utils)
# bedtools should also be in the path to run this

####Using the SCENT Object:
SCENT_obj <- CreateSCENTObj(rna = mrna, atac = atac, meta.data = meta,
                            covariates = c("log(nCount_RNA)","percent.mito"), 
                            celltypes = "newCT")
str(SCENT_obj)


genebed_loc <- "./RData/Data/GeneBody_500kb_margin_chr.bed" # path to the file which we will place in the github

SCENT_obj <- CreatePeakToGeneList(SCENT_obj, genebed = genebed_loc,
                                  nbatch = 1000,tmpfile="./temporary_atac_peak.bed",
                                  intersectedfile="./temporary_atac_peak_intersected.bed.gz")
str(SCENT_obj, max.level = 2)

saveRDS(SCENT_obj, file = "./Testing/Output/SCENT_obj.rds") #Takes a couple minutes.


####Main Q.):
#How to deal with differences in types for peak.info.
SCENT_obj2 <- SCENT_obj
SCENT_obj2@peak.info <- SCENT_obj2@peak.info[["1001544.0"]]



#####2nd way#######
####Using the SCENT Object:
SCENT_obj <- CreateSCENTObj(rna = mrna, atac = atac, meta.data = meta,
                            peak.info = gene_peak,
                            covariates = c("log(nCount_RNA)","percent.mito"), 
                            celltypes = "newCT")

#Of the set of peak gene pairs: example of pick a specific set of pairs to test: 
#Example: (first 10 gene-peak pairs)
SCENT_obj@peak.info <- SCENT_obj@peak.info[1:10,]

#Run SCENT algorithm of Tnk cell type and use 6 cores for parallelization:
SCENT_obj <- SCENT_algorithm(SCENT_obj, "Tnk", 6)

write.table(SCENT_obj@SCENT.result, file = output, quote=F, row=F, sep="\t")