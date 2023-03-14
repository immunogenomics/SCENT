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

####### INPUTS 
#Obtain arguments: (from Cluster)
node = commandArgs(trailingOnly = T)[1] # JOB ARRAY number: node usage
cores = commandArgs(trailingOnly = T)[2] # Number of Cores
SCENTobj_rds = commandArgs(trailingOnly = T)[3] # RDS object file type
celltype = commandArgs(trailingOnly = T)[4] # CellType
output_dir  = commandArgs(trailingOnly = T)[5] # Output of each text file to a specific folder


###Example of inputs:
# node <- 1
# cores <- 6
# celltype <- "Tnk"
# SCENTobj_rds <- "./Testing/Output/SCENT_obj.rds"
# output_dir <- "./Testing/Output/"


#### Load:
SCENT_obj <- readRDS(SCENTobj_rds)

#### Get the corresponding dataframe from the list:
SCENT_obj@peak.info <- SCENT_obj@peak.info[[node]]

#### Run SCENT algorithm of Tnk cell type and use 6 cores for parallelization:
SCENT_obj <- SCENT_algorithm(SCENT_obj, celltype, cores)

#### Output SCENT results for each gene-peak pair block.
filename <- paste0(output_dir,"SCENTresult_",node,".txt")

write.table(SCENT_obj@SCENT.result, file = filename, row.names = F, col.names = T)