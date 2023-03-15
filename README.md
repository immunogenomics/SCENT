

# SCENT

Single-Cell ENhancer Target gene mapping using multimodal data with ATAC + RNA

(*beta version*)

The manuscript will soon appear at medRxiv! (Sakaue et al. "**Tissue-specific enhancer-gene maps from multimodal single-cell data identify causal disease alleles**")



### Overview

SCENT uses single-cell multimodal data (e.g., 10X Multiome RNA/ATAC) and links ATAC-seq peaks (putative enhancers) to their target genes by modeling association between chromatin accessibility and gene expression across individual single cells.

<div align="center">
<img src="https://github.com/immunogenomics/SCENT/blob/e0f14cd59a7a148d94383e2a825f3546e2045d41/fig/cover_image.png" width=90%>
</div>


We use Poisson regression to associate gene expression (raw) count and (binarized) peak accessibility, and estimate errors in coefficients by bootstrapping framework to control for type I error.




### Installation of SCENT Package

You can install the development version of SCENT from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("immunogenomics/SCENT")
```


### Requirements

The SCENT package will automatically install CRAN R packages. The packages below will go into your `R`.

- `methods`
- `data.table`
- `lme4`
- `stringr`
- `boot`
- `MASS`
- `Matrix`
- `parallel`

The SCENT package also requires command-line tool, bedtools, for developing a list of: gene-peak pair dataframes to parallelize through.
- `https://github.com/arq5x/bedtools2`


### Example usage

Vignettes are posted in this github repo to show 2 potential uses of the SCENT package.

### 1.) Using SCENT interactively for testing small sets of gene-peak associations

`SCENT_interactive.Rmd` vignette contains an example of using the SCENT package to generate results on small sets of gene-peak associations. 

In summary, the main functionality is the SCENT object construction:

```r
SCENT_obj <- CreateSCENTObj(rna = mrna, atac = atac, meta.data = meta,
                            peak.info = gene_peak,
                            covariates = c("log(nCount_RNA)","percent.mito"), 
                            celltypes = "newCT")
```

Followed by SCENT algorithm:
```r
SCENT_obj <- SCENT_algorithm(object = SCENT_obj, celltype = "Tcells", ncores = 6)
```
Where the user specifies the celltype for association analysis and the number of cores for parallelization of bootstrapping.

The output of the SCENt algorithm will be contained in the field:
```r
SCENT_obj@SCENT.result
```
which can be saved as a textfile for further downstream analysis.


Further information on Inputs and Outputs of SCENT are detailed below:

#### Inputs To SCENT Object:

| #    | Argument name (format)       | Descriptions                                                 |
| ---- | ---------------------------- | ------------------------------------------------------------ |
| 1    | rna (.rds)            | A gene-by-cell matrix from multimodal RNA-seq data. This is a raw count matrix without any normalization. The column names should be the gene names used in the `file_gene_peak_tested` file. This can be sparse matrix format. |
| 2    | atac (.rds)           | A peak-by-cell matrix from multimodal ATAC-seq data. The row names should be the peak names used in the `file_gene_peak_tested` file. The column names are the cell names which should be the same names used in `rna_matrix` and the `cell`column of `metafile`. The matrix may not be binarized while it will be binarized within the function. This can be sparse matrix format. |
| 3    | meta.data (.txt)              | A metadata for cells (rows are cells, and cell names should be in the column named as "cell"; see below example). Additionally, this text should include covariates to use in the model. Examples include: % mitochondrial reads, nUMI, sample, and batch as covariates. |
| 4    | peak.info (.txt) | A textfile indicating which gene-peak pairs you want to test in this chunk (see below example). We highly recommend splitting gene-peak pairs into many chunks to increase computational efficiency. |
| 5    | covariates (character)        | A vector of character fields that denote the covariates listed in the meta.data. For example, a set of covariates can be: %mitochondrial reads, nUMI, sample, and batch. Additionally the user can specify transformations to the covariates such as log transformation on nUMI counts for direct usage in the SCENT algorithm invoking poisson glm. |
| 6    | cell_type (character)        | A name of the cell type you want to test in this association analysis. This should be corresponding to the `celltype` column of the `metafile`. |

Alternatives: The peak.info field can be left blank and created using the CreatePeakToGeneList function in the SCENT package. This function requires the user to specify a bed file that specifies ~500 kb windows of multiple gene loci to identify cis gene-peak pairs to test.



#### Example Formats: 
The example format of  `file_gene_peak_tested` file in text format.

```bash
$ head ${file_gene_peak_tested} 
A1BG	chr19-57849279-57850722
A1BG	chr19-57888160-57889279
A1BG	chr19-57915851-57917093
A1BG	chr19-57934422-57935603	
```

We usually only select peaks of which the center falls within 500 kb from the target gene (*cis* analysis). Also, while we have a function to QC peaks and genes so that they are present in at least 5% of all cells within `SCENT.R`, it is more efficient to only include these QCed peaks and genes in  `file_gene_peak_tested`  to reduce the number of tests.


The example format of  `metafile` file in rds format.

```r
meta <- readRDS(metafile)
head(meta)

                                 cell nUMI percent_mito   sample   batch
AAACAGCCAAGGAATC-1 AAACAGCCAAGGAATC-1 8380   0.01503428 sample_1 batch_a
AAACAGCCAATCCCTT-1 AAACAGCCAATCCCTT-1 3771   0.02207505 sample_1 batch_a
AAACAGCCAATGCGCT-1 AAACAGCCAATGCGCT-1 6876   0.01435579 sample_1 batch_a
AAACAGCCACACTAAT-1 AAACAGCCACACTAAT-1 1733   0.03881841 sample_1 batch_a
AAACAGCCACCAACCG-1 AAACAGCCACCAACCG-1 5415   0.01600768 sample_1 batch_a
AAACAGCCAGGATAAC-1 AAACAGCCAGGATAAC-1 2759   0.02485340 sample_1 batch_a
                   celltype
AAACAGCCAAGGAATC-1    Tcell
AAACAGCCAATCCCTT-1    Tcell
AAACAGCCAATGCGCT-1    Tcell
AAACAGCCACACTAAT-1    Tcell
AAACAGCCACCAACCG-1    Tcell
AAACAGCCAGGATAAC-1    Tcell
```


### Output of SCENT (SCENT.result field)

```bash
$ head ${file_output}
gene	peak	beta	se	z	p	boot_basic_p
A1BG	chr19-57849279-57850722	0.587060911718621	0.227961010352348	2.57526894977009	0.0100162168431262	0.0192
A1BG	chr19-57888160-57889279	-0.0842330294127105	0.232845263030106	-0.3617553920425660.717534829528597	0.688
A1BG	chr19-57915851-57917093	-0.00971211792633636	0.225020479431863	-0.0431610400566990.965573161660521	1
A1BG	chr19-57934422-57935603	0.0136752444069743	0.249810124611214	0.05474255468331160.956343566437322	0.968
```

Each column indicates ...

| Column       | Descriptions                                                 |
| ------------ | ------------------------------------------------------------ |
| gene         | The gene(-peak) pair in each test statistics                 |
| peak         | The (gene-)peak pair in each test statistics                 |
| beta         | The regression coefficient from primary Poisson regression   |
| se           | The standard error  from primary Poisson regression          |
| z            | The Z score from primary Poisson regression                  |
| p            | The raw p value from primary Poisson regression              |
| boot_basic_p | The bootstrap p value calculated from bootstrapping analyses |



### 2.) Using SCENT with parallelized jobs.


`SCENT_parallelization.R` is the code necessary for running parallelized SCENT jobs.
This code needs a `SCENT Object rds` file that contains a list of gene-peak pairs. 
To generate this object please follow the SCENT_parallelize.Rmd vignette file.

The corresponding bash script `parallelizedSCENT.sh` contains a parallelization scheme that is 
dependent on the amount of gene-peak pair batches that is user defined (for context please refer to the
SCENT_parallelize.Rmd vignette). The main part of the bash script contains the line:

```bash
Rscript SCENT_parallelization.R $LSB_JOBINDEX ${num_cores} ${file_SCENT_obj} ${celltype} ${output_dir}
```

Arguments in the bash file are user specified as follows:

|#      | Argument Name | Descriptions |
| ----  | ------------- | ------------ |
|1    | LSB_JOBINDEX   | jobarray index specified by BSUB -J |
|2    | num_cores      | number of cores (ex. 6) to parallelize to the SCENT algorithm |
|3    | file_SCENT_obj | SCENT object that contains atac_matrix, rna_matrix, metafile, peak_gene_list, etc. To run the SCENT algorithm |
|4    | celltype       |  User specified celltype (ex. "Tcells") to run the SCENT algorithm |
|5    | output_dir     | User specified directory to output the SCENT results to aggregate once parallization is completed.
|


### Contact

Saori Sakaue ssakaue@broadinstitute.org
