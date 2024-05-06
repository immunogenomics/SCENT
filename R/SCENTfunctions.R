## define functions
#' Interpolate a p-value from quantiles that should be "null scaled"
#'
#' @param q bootstrap quantiles, centered so that under the null, theta = 0
#' @return two-sided p-value
#' @export
interp_pval = function(q) {
  R = length(q)
  tstar = sort(q)
  zero = findInterval(0, tstar)
  if(zero == 0 || zero == R) return(2/R) # at/beyond extreme values
  pval = 2*min(zero/R, (R-zero)/R)
  pval
}


#' Derive a p-value from a vector of bootstrap samples using the "basic" calculation
#'
#' @param obs observed value of parameter (using actual data)
#' @param boot vector of bootstraps
#'
#' @return p-value
#' @export
basic_p = function(obs, boot, null = 0){
  interp_pval(2*obs - boot - null)
}


#' Perform poisson regression: exprs ~ peak + covariates
#'
#' @param data contains expr values and associated peak and covariates for a gene.
#' @param idx rows of the data to use: argument for boot function (bootstrapping)
#' @param formula user defined formula based on initialization in CreateSCENTObj Constructor
#'
#' @return vector: (coefficient of the peak effect on gene, variance of peak effect on gene)
#' @export
assoc_poisson = function(data, idx = seq_len(nrow(data)), formula){
  gg = glm(formula, family = 'poisson', data = data[idx,,drop = FALSE])
  c(coef(gg)['atac'], diag(vcov(gg))['atac'])
}


#' Perform Mixed Linear regression: exprs ~ peak + (1|RandomEffect) + FixedEffect_covariates
#'
#' @param data contains expr values and associated peak and covariates for a gene.
#' @param idx rows of the data to use: argument for boot function (bootstrapping)
#' @param formula user defined formula based on initialization in CreateSCENTObj Constructor
#'
#' @return vector: (coefficient of the peak effect on gene, variance of peak effect on gene)
#' @export
assoc_lmm = function(data, idx = seq_len(nrow(data)), formula){
  gg = lme4::lmer(formula, data = data[idx,,drop = FALSE])
  c( lme4::fixef(gg)['atac'], diag( as.matrix( lme4::vcov.merMod( gg, fixed.only =T ) ) )['atac'] )
}


#' Perform negative binomial regression: exprs ~ peak + covariates
#'
#' @param data contains expr values and associated peak and covariates for a gene.
#' @param idx rows of the data to use: argument for boot function (bootstrapping)
#' @param formula user defined formula based on initialization in CreateSCENTObj Constructor
#'
#' @return vector: (coefficient of the peak effect on gene, variance of peak effect on gene)
#' @export
assoc_negbin = function(data, idx = seq_len(nrow(data)), formula){
  gg = glm.nb(formula, data = data[idx,,drop = FALSE])
  c(coef(gg)['atac'], diag(vcov(gg))['atac'])
}

#' Validity and Type Checking for CreateSCENTObject Constructor
#'
#' @param object SCENT object constructed from class CreateSCENTObject
#'
#' @return None OR Errors dependent on if the object follows the guidelines for SCENT
#' RNA: matrix of (genes x cells)
#' ATAC: matrix of (peaks x cells)
#' @export
check_dimensions <- function(object){
  errors <- character()

  #Check dimensionality of cells:
  num_cells_rna <- lengths(object@rna@Dimnames)[2]
  num_cells_atac <- lengths(object@atac@Dimnames)[2]

  num_genes <- lengths(object@rna@Dimnames)[1]
  num_peaks <- lengths(object@atac@Dimnames)[1]

  #Check if the number of cells match between rna and atac matrix.
  if(num_cells_rna != num_cells_atac){
    msg <- paste("Error: The num of cells in scRNA matrix is: ", num_cells_rna,
                 " and the num of cells in scATAC matrix is: ", num_cells_atac,
                 ". These should EQUAL EACH OTHER, please check to make sure",
                 " both matrices for scRNA and scATAC are read in as",
                 " (genes x cells) and (peaks x cells), respectively. ")
    errors <- c(errors, msg)
  }


  #Most likely the number of peaks is greater than the number of genes if not WARN.
  if(num_peaks < num_genes){
    warning(paste("Warning: in general there are more peaks found through ATAC",
                   " than genes. Currently you have number of peaks =", num_peaks,
                   " and number of genes =",num_genes))
  }

  #If peak.info is present check the following:
  if(!(length(object@peak.info) == 0)){
    #Check if genes correspond between rna matrix and peak.info dataframe:
    if(!all(object@peak.info[[1]] %in% object@rna@Dimnames[[1]])){
      msg <- paste("The gene names in the peak.info dataframe is NOT a subset of the gene names in",
                   " the scRNA matrix")
      errors <- c(errors, msg)
    }


    #Check if peaks correspond between atac matrix and peak.info dataframe:
    if(!all(object@peak.info[[2]] %in% object@atac@Dimnames[[1]])){
      msg <- paste("The peak ids in the peak.info dataframe is NOT a subset of the peak names in",
                   " the scATAC matrix")
      errors <- c(errors, msg)
    }
  }


  ###Additional things to check:
  #Check if meta.data table with covariates has the correct cell column names
  #Check if covariates are in the columns of meta.data
  if (length(errors) == 0) TRUE else errors
}



#' SCENT Class Constructor
#'
#' @slot rna dgCMatrix. scRNAseq matrix read as a sparse matrix
#' @slot atac dgCMatrix. scATACseq matrix read as a sparse matrix
#' @slot meta.data data.frame. Metadata table with covariates and a cell ID column ("cell")
#' @slot peak.info data.frame. Dataframe that contains gene-peak pairs for SCENT to search through
#' @slot peak.info.list list. List of dataframes that contain gene-peak pairs to parallelize through
#' @slot covariates character. Assign covariates that are needed for the analysis. Must be names that are in the columns of meta.data
#' @slot random.effects character. Assign covariates that are to be treated as random effects if using randome effects modeling. if not set to NULL
#' @slot celltypes character. Assign celltype column from meta.data
#' @slot SCENT.result data.frame. Initialized as empty. Becomes a table of resultant significant gene peak pairs
#'
#' @return SCENT object to use for further analysis
#' @export
CreateSCENTObj <- setClass(
  Class = "SCENT",
  slots = c(
    rna = 'dgCMatrix',
    atac = 'dgCMatrix',
    meta.data = 'data.frame',
    peak.info = 'data.frame',  ###Must be gene (1st column) then peak (2nd column)
    peak.info.list = 'list',
    covariates = 'character',
    random.effects = 'character',
    celltypes = 'character',
    SCENT.result = 'data.frame'
  ),
  validity = check_dimensions
)

#' SCENT Algorithm: Poisson Regression with Empirical P-values through Bootstrapping.
#'
#' @param object SCENT object
#' @param celltype character. User specified cell type defined in celltypes column of meta.data
#' @param ncores numeric. Number of cores to use for Parallelization
#' @param regr character. Regression type: "poisson" or "negbin" for Poisson regression and Negative Binomial regression, respectively
#' @param bin logical. TRUE to binarize ATAC counts. FALSE to NOT binarize ATAC counts
#'
#' @return SCENT object with updated field SCENT.results
#' @export
SCENT_algorithm <- function(object, celltype, ncores, regr = "poisson", bin = TRUE){
  res <- data.frame()
  for (n in 1:nrow(object@peak.info)){ ####c(1:nrow(chunkinfo))
    message( paste0( "Working on ", n, " of ", nrow(object@peak.info), " Potential Interactions."))
    gene <- as.character( object@peak.info[n,1] ) #GENE is FIRST COLUMN OF PEAK.INFO
    this_peak <- as.character( object@peak.info[n,2] ) #PEAK is SECOND COLUMN OF PEAK.INFO
    atac_target <- data.frame(cell = colnames(object@atac), atac = object@atac[this_peak,])

    #binarize peaks:
    if(bin){
      atac_target[atac_target$atac>0,]$atac<-1
    }

    mrna_target <- object@rna[gene,]
    df <- data.frame(cell=names(mrna_target),exprs=as.numeric(mrna_target))
    df<-merge(df,atac_target,by="cell")
    df<-merge(df,object@meta.data,by="cell")

    df2 <- df[df[[object@celltypes]] == celltype,]

    nonzero_m  <- length( df2$exprs[ df2$exprs > 0] ) / length( df2$exprs )
    nonzero_a  <- length( df2$atac[ df2$atac > 0] ) / length( df2$atac )
    if(nonzero_m > 0.05 & nonzero_a > 0.05){
      #Run Regression Once Before Bootstrapping:
      res_var <- "exprs"
      pred_var <- c("atac", object@covariates)
      if( length( object@random.effects ) == 0 ){
        formula <- as.formula(paste(res_var, paste(pred_var, collapse = "+"), sep = "~"))  
      }else{
        # Build Random Effects Model
        pred_fixed <- pred_var[ !(pred_var %in% object@random.effects) ]
        pred_rand <- paste0( '(1|', pred_var[ pred_var %in% object@random.effects ], ')')
        formula <- as.formula(paste(res_var, paste(c(pred_rand, pred_fixed), collapse = "+"), sep = "~"))
      }
      
      #Estimated Coefficients Obtained without Bootstrapping:
      if(regr == "poisson"){
        base = glm(formula, family = 'poisson', data = df2)
        coefs<-summary(base)$coefficients["atac",]
        assoc <- assoc_poisson
      } else if (regr == "negbin"){
        base = glm.nb(formula, data = df2)
        coefs<-summary(base)$coefficients["atac",]
        assoc <- assoc_negbin
      } else if (regr == "lmm"){
        base = lme4::lmer(formula, data = df2)
        coefs <- summary(base)$coefficients["atac",]
        coefs['Pr(>|z|)'] <- 2 * stats::pt(abs(as.numeric(coefs['t value'])), stats::df.residual(base), lower.tail = FALSE)
        assoc <- assoc_lmm
      }
      
      ###Iterative Bootstrapping Procedure: Estimate the Beta coefficients and associate a 2-sided p-value.
      bs = boot::boot(df2,assoc, R = 100, formula = formula, stype = 'i', parallel = "multicore", ncpus = ncores)
      p0 = basic_p(bs$t0[1], bs$t[,1])
      if(p0<0.1){
        bs = boot::boot(df2,assoc, R = 500, formula = formula,  stype = 'i', parallel = "multicore", ncpus = ncores)
        p0 = basic_p(bs$t0[1], bs$t[,1])
      }
      if(p0<0.05){
        bs = boot::boot(df2,assoc, R = 2500, formula = formula,  stype = 'i', parallel = "multicore", ncpus = ncores)
        p0 = basic_p(bs$t0[1], bs$t[,1])
      }
      if(p0<0.01){
        bs = boot::boot(df2,assoc, R = 25000, formula = formula,  stype = 'i', parallel = "multicore", ncpus = ncores)
        p0 = basic_p(bs$t0[1], bs$t[,1])
      }
      #if(p0<0.001){
      #  bs = boot::boot(df2,assoc, R = 50000, formula = formula, stype = 'i', parallel = "multicore", ncpus = ncores)
      #  p0 = basic_p(bs$t0[1], bs$t[,1])
      #}
      out <- data.frame(gene=gene,peak=this_peak,beta=coefs[1],se=coefs[2],z=coefs[3],p=coefs[4],boot_basic_p=p0)
      message( paste0( n, " Result: ", gene, " ", this_peak, " ", coefs[1], " ", coefs[2], " ", coefs[3], " ", coefs[4], " ", p0))
      res<-rbind(res,out)
    }
  }

  #Update the SCENT.result field of the constructor in R:
  object@SCENT.result <- res
  return(object)
}



#' Creating Cis Gene-Peak Pair Lists to Parallelize Through
#'
#' @param object SCENT object
#' @param genebed character. File directory for bed file that contains 500 kb windows for each gene
#' @param nbatch numeric. Number of batches to produce: Length of the list
#' @param tmpfile character. Location of temporary file.
#' @param intersectedfile character. Location of intersected file.
#'
#' @return SCENT object with updated field of peak.info.list
#' @export
CreatePeakToGeneList <- function(object,genebed="/path/to/GeneBody_500kb_margin.bed",nbatch,tmpfile="./temporary_atac_peak.bed",intersectedfile="./temporary_atac_peak_intersected.bed.gz"){
  peaknames <- rownames(object@atac) # peak by cell matrix
  peaknames_r <- gsub(":","-",peaknames) # in case separator included ":"
  peaknames_r <- gsub("_","-",peaknames_r) # in case separator included "_"
  peak_bed <- data.frame(chr = str_split_fixed(peaknames_r,"-",3)[,1], start = str_split_fixed(peaknames_r,"-",3)[,2], end = str_split_fixed(peaknames_r,"-",3)[,3], peak=peaknames)
  write.table(peak_bed,tmpfile,quote=F,row=F,col=F,sep="\t")
  system(paste("bedtools intersect -a",genebed,"-b ",tmpfile, " -wa -wb -loj | gzip -c >", intersectedfile))
  system(paste("rm ", tmpfile))
  d <- fread(intersectedfile,sep="\t")
  d<-data.frame(d)
  d <- d[d$V5 != ".",]

  #Obtain gene to peak pairs.
  cis.g2p <- d[c("V4","V8")]
  colnames(cis.g2p) <- c("gene","peak")

  cis.g2p$index <- 1:nrow(cis.g2p)
  cis.g2p$batch_index <- cut2(cis.g2p$index, g = nbatch, levels.mean = TRUE)
  cis.g2p_list <- split(cis.g2p, f = cis.g2p$batch_index)
  cis.g2p_list <- lapply(cis.g2p_list, function(x) x[(names(x) %in% c("peak", "gene"))])
  names(cis.g2p_list) <- 1:length(cis.g2p_list)
  # Update the SCENT.peak.info field of the constructor in R:
  object@peak.info.list <- cis.g2p_list
  return(object)
}




