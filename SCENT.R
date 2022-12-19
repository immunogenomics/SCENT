library(data.table)
library(lme4)
library(stringr)
library(boot)
library(MASS)
library(Matrix)

## define functions

#' Interpolate a p-value from quantiles that should be "null scaled"
#'
#' @param q bootstrap quantiles, centered so that under the null, theta = 0
#' @return two-sided p-value
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
basic_p = function(obs, boot, null = 0){
  interp_pval(2*obs - boot - null)
}


assoc_poisson = function(data, idx = seq_len(nrow(data))){
  gg = glm(exprs ~ atac + percent_mito + log(nUMI) + sample + batch, family = 'poisson', data = data[idx,,drop = FALSE])
  c(coef(gg)['atac'], diag(vcov(gg))['atac'])
}


## load data

input_atac <- commandArgs(trailingOnly = T)[1]
input_mrna <- commandArgs(trailingOnly = T)[2]
input_meta <- commandArgs(trailingOnly = T)[3]
gene_peak <- commandArgs(trailingOnly = T)[4]
celltype <- commandArgs(trailingOnly = T)[5]
output <- commandArgs(trailingOnly = T)[6]

options(stringsAsFactors = F)

atac <- readRDS(input_atac)
mrna <- readRDS(input_mrna)
meta <- readRDS(input_meta)

chunkinfo <- readRDS(gene_peak)
colnames(chunkinfo)<-c("gene","peak")

res<-data.frame()
for(n in c(1:nrow(chunkinfo))){
    gene=chunkinfo$gene[n]
    this_peak=chunkinfo$peak[n]
    atac_target<-data.frame(cell=colnames(atac),atac=as.numeric(atac[this_peak,]))
    # binarize peaks
    atac_target[atac_target$atac>0,]$atac<-1
    mrna_target<-mrna[gene,]
    df <- data.frame(cell=names(mrna_target),exprs=as.numeric(mrna_target))
    df<-merge(df,atac_target,by="cell")
    df<-merge(df,meta,by="cell")
    df2 <- df[df$celltype==celltype,]
    nonzero_m  <- length( df2$exprs[ df2$exprs > 0] ) / length( df2$exprs )
    nonzero_a  <- length( df2$atac[ df2$atac > 0] ) / length( df2$atac )
    if(nonzero_m > 0.05 & nonzero_a > 0.05){
      # poisson
      base = glm(exprs ~ atac + percent_mito + log(nUMI) + sample + batch, family = 'poisson', data = df2)
      coefs<-summary(base)$coefficients["atac",]
      bs = boot::boot(df2,assoc_poisson, R = 100, stype = 'i')
      p0 = basic_p(bs$t0[1], bs$t[,1])
      if(p0<0.1){
        bs = boot::boot(df2,assoc_poisson, R = 500, stype = 'i')
        p0 = basic_p(bs$t0[1], bs$t[,1])
      }
      if(p0<0.05){
        bs = boot::boot(df2,assoc_poisson, R = 2500, stype = 'i')
        p0 = basic_p(bs$t0[1], bs$t[,1])
      }
      if(p0<0.01){
        bs = boot::boot(df2,assoc_poisson, R = 25000, stype = 'i')
        p0 = basic_p(bs$t0[1], bs$t[,1])
      }
      if(p0<0.001){
        bs = boot::boot(df2,assoc_poisson, R = 50000, stype = 'i')
        p0 = basic_p(bs$t0[1], bs$t[,1])
      }
      out <- data.frame(gene=gene,peak=this_peak,beta=coefs[1],se=coefs[2],z=coefs[3],p=coefs[4],boot_basic_p=p0)
      res<-rbind(res,out)
    }}

write.table(res, output, quote=F, row=F, sep="\t")
