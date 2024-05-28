#!/bin/bash

#BSUB -J SCENT[1-100] #Number of job arrays dependent on number of peak-gene pair batches
#BSUB -q big     #node for multi-parallelized threading and nodes
#BSUB -M 18000   #18 GB
#BSUB -n 6       #6 cores
#BSUB -o Output_%J_%I.out  #output file %J is job %I is job array index
#BSUB -e Error_%J_%I.err   #error file %J is job %I is job array index


module load R
Rscript SCENT_parallelization.R $LSB_JOBINDEX ${num_cores} ${file_SCENT_obj} ${celltype} ${regr} ${bin} ${output_dir}

