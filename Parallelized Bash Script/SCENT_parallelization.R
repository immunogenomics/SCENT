#Libraries to Load:
library(SCENT)

####### INPUTS
#Obtain arguments: (from Cluster)
node = as.integer(commandArgs(trailingOnly = T)[1]) # integer. JOB ARRAY number: node usage
cores = as.integer(commandArgs(trailingOnly = T)[2]) # integer. Number of Cores
SCENTobj_rds = commandArgs(trailingOnly = T)[3] # character. RDS object file type
celltype = commandArgs(trailingOnly = T)[4] # character. CellType
regr = commandArgs(trailingOnly = T)[5] # character. Regression Type
bin = as.logical(commandArgs(trailingOnly = T)[6]) # logical. Binarize ATAC counts
output_dir  = commandArgs(trailingOnly = T)[7] # character. Output of each text file to a specific folder

###Example of inputs from the bash script: parallelizedSCENT.sh
# node <- 1
# cores <- 6
# celltype <- "Tnk"
# SCENTobj_rds <- "./Testing/Output/SCENT_obj.rds"
# output_dir <- "./Testing/Output/"


#### Load:
SCENT_obj <- readRDS(SCENTobj_rds)

#### Get the corresponding dataframe from the list:
SCENT_obj@peak.info <- SCENT_obj@peak.info.list[[node]]

#### Run SCENT algorithm of Tnk cell type and use 6 cores for parallelization:
SCENT_obj <- SCENT_algorithm(SCENT_obj, celltype, cores, regr, bin)

#### Output SCENT results for each gene-peak pair block.
filename <- paste0(output_dir,"/SCENTresult_",node,".txt")

write.table(SCENT_obj@SCENT.result, file = filename, row.names = F, col.names = T, quote = F)
