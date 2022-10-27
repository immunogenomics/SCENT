

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




### Installation

In order to get started with `SCENT`, you can clone this repository.

```{bash}
$ git clone https://github.com/immunogenomics/SCENT
$ cd ./SCENT
```



### Requirements

You should have installed these packages into your `R`.

- `data.table`
- `lme4`
- `stringr`
- `boot`
- `MASS`
- `Matrix`



### Example usage

`SCENT.R` is the basic code necessary for running SCENT.

```bash
$ Rscript SCENT.R ${atac_matrix} ${rna_matrix} ${metafile} ${file_gene_peak_tested} ${cell_type} ${file_output}
```

SCENT takes several inputs that have to be formatted as follows.

### Input

| #    | Argument name (format)       | Descriptions                                                 |
| ---- | ---------------------------- | ------------------------------------------------------------ |
| 1    | atac_matrix (.rds)           | A peak-by-cell matrix from multimodal ATAC-seq data. The row names should be the peak names used in the `file_gene_peak_tested` file. The column names are the cell names which should be the same names used in `rna_matrix` and the `cell`column of `metafile`. The matrix may not be binarized while it will be binarized within the function. This can be sparse matrix format. |
| 2    | rna_matrix (.rds)            | A gene-by-cell matrix from multimodal RNA-seq data. This is a raw count matrix without any normalization. The column names should be the gene names used in the `file_gene_peak_tested` file. This can be sparse matrix format. |
| 3    | metafile (.rds)              | A metadata for cells (rows are cells, and cell names should be in the column named as "cell"; see below example). Currently the model includes % mitocondrial reads, nUMI, sample, and batch as covariates (see line 35 and 74). |
| 4    | file_gene_peak_tested (.txt) | A textfile indicating which gene-peak pairs you want to test in this chunk (see below example). We highly recommend splitting gene-peak pairs into many chunks to increase computational efficiency. |
| 5    | cell_type (character)        | A name of the cell type you want to test in this association analysis. This should be corresponding to the `celltype` column of the `metafile`. |
| 6    | file_output (character)      | A name of the output file.                                   |



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

`cell` , `nUMI` , and `celltype` columns must exist.

You can modify the line 35 and 74 of `SCENT.R` to edit covariate names and/or include other covariates into the model.



### Output

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



### Contact

Saori Sakaue ssakaue@broadinstitute.org
