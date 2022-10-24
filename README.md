# SCENT
Single-Cell ENhancer Target gene mapping using multimodal data with ATAC + RNA

(beta version)



### Overview

SCENT links ATAC-seq peaks (putative enhancers) with target gene by modeling association between chromatin accessibility and gene expression across individual single cells.

<div align="center">
<img src="https://github.com/immunogenomics/SCENT/blob/e0f14cd59a7a148d94383e2a825f3546e2045d41/fig/cover_image.png" width=90%>
</div>





### Installation

In order to get started with `SCENT`, you can clone this repository.

```{bash}
$ git clone https://github.com/immunogenomics/SCENT
$ cd ./SCENT
```



### Example usage

`SCENT.R` is the basic code necessary for running SCENT.

```bash
$ Rscript SCENT.R ${atac_matrix} ${rna_matrix} ${metafile} ${file_gene_peak_tested} ${cell_type} ${file_output}
```

SCENT takes several inputs that have to be formatted as follows.

### Input

| #    | Argument name (format)       | Descriptions                                                 |
| ---- | ---------------------------- | ------------------------------------------------------------ |
| 1    | atac_matrix (.rds)           | A peak-by-cell matrix from multimodal ATAC-seq data.         |
| 2    | rna_matrix (.rds)            | A gene-by-cell matrix from multimodal RNA-seq data.          |
| 3    | metafile (.rds)              | A metadata for cells (rows are cells, and cell names should be in the column named as "cell"). Currently the model includes % mitocondrial reads, nUMI, sample, and batch as covariates. |
| 4    | file_gene_peak_tested (.txt) | A textfile indicating which gene-peak pairs you want to test in this chunk (see below example). We highly recommend splitting gene-peak pairs into many chunks to increase computational efficiency. |
| 5    | cell_type (character)        | A name of the cell type you want to test in this             |
| 6    | file_output (character)      | A name of the output file.                                   |

```bash
$ head ${file_gene_peak_tested} 
A1BG	chr19-57849279-57850722
A1BG	chr19-57888160-57889279
A1BG	chr19-57915851-57917093
A1BG	chr19-57934422-57935603	
```



### Output

```bash
$ head ${file_output}

gene	peak	beta	se	z	p	boot_basic_p
A1BG	chr19-57849279-57850722	0.587060911718621	0.227961010352348	2.57526894977009	0.0100162168431262	0.0192
A1BG	chr19-57888160-57889279	-0.0842330294127105	0.232845263030106	-0.3617553920425660.717534829528597	0.688
A1BG	chr19-57915851-57917093	-0.00971211792633636	0.225020479431863	-0.0431610400566990.965573161660521	1
A1BG	chr19-57934422-57935603	0.0136752444069743	0.249810124611214	0.05474255468331160.956343566437322	0.968
```



| Column       | Descriptions                                                 |
| ------------ | ------------------------------------------------------------ |
| gene         | The gene(-peak) pair in each test statistics                 |
| peak         | The (gene-)peak pair in each test statistics                 |
| beta         | The regression coefficient from primary Poisson regression   |
| se           | The standard error  from primary Poisson regression          |
| z            | The Z score from primary Poisson regression                  |
| p            | The raw p value from primary Poisson regression              |
| boot_basic_p | The bootstrap p value calculated from bootstrapping analyses |

