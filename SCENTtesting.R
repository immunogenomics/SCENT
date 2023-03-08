library(data.table)
library(lme4)
library(stringr)
library(boot)
library(MASS)
library(Matrix)
library(parallel) #core parallelization.
library(methods) #S4 class OOP methods.


#Get functions and classes:
source("SCENTfunctions_update.R")

input_atac <- "./RData/Data/pbmc_multimodal.atac.rds"
input_mrna <- "./RData/Data/pbmc_multimodal.rna.rds"
input_meta <- "./RData/Data/pbmc_multimodal.meta.rds"
gene_peak <- "./RData/Data/qced_Tnk.G2P.txt"
celltype <- "Tnk"
output <- "./Testing/Output/test_output.txt"
num_cores <- 6
covariates <- "log(nUMI) + sample + batch"

options(stringsAsFactors = F)

atac <- readRDS(input_atac)
mrna <- readRDS(input_mrna)
meta <- readRDS(input_meta)
chunkinfo <- read.table(gene_peak)
colnames(chunkinfo) <- c("gene","peak")



#####TESTING of the SCENT Object:
###Instantiate the Object:
SCENT_obj <- CreateSCENTObj(rna = mrna, atac = atac, meta.data = meta,
                            peak.info = chunkinfo,
                            covariates = c("log(nCount_RNA)","percent.mito"), 
                            celltypes = "newCT")

####Testing if errors can be caught:
#Swapping: dataframes
SCENT_obj_error1 <- CreateSCENTObj(rna = atac, atac = mrna, meta.data = meta,
                                   covariates = covariates, peak.info = chunkinfo)

#Subtract a column from one dataframe:
SCENT_obj_error2 <- CreateSCENTObj(rna = mrna, atac = atac[,-1], meta.data = meta,
                                   covariates = covariates, peak.info = chunkinfo)

#Swap chunkinfo columns:
chunkinfo[,c(1,2)]  <- chunkinfo[,c(2,1)]
SCENT_obj_error3 <- CreateSCENTObj(rna = mrna, atac = atac, meta.data = meta,
                                   covariates = covariates, peak.info = chunkinfo)