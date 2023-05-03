# scASfind - mining alternative splicing patterns in single-cell atlas

## Introduction
The emergence of large single cell RNA-seq(scRNA-seq) datasets facilitates biological discoveries at single-cell resolution. The information encoded in single cell atlases should be readily accessible to users from all disciplines. The Hemberg lab has developed `scfind`, a search engine for gene expression patterns in cell atlases [scfind](https://github.com/hemberg-lab/scfind).  

`scASfind` is built on top of `scfind` and it adapts `scfind` to **alternative splicing** analysis. Through querying with splicing features or cell types, `scASfind` conducts rapid searching and returns feature-enriched cell types or cell type-specific splicing signatures, respectively. `scASfind` also discovers cell type-specific mutually exclusive exon pairs and clustered spliced-in or spliced-out of a block of exons by enumerating over all possible combinations of splicing events. 

Quantification of single cell alternative splicing events is performed by [Whippet](https://github.com/timbitz/Whippet.jl) integrated in the workflow [MicroExonator](https://github.com/hemberg-lab/MicroExonator). scASfind builds custom indices for any single-cell dataset using the output files of `Whippet`.

## Installation
To install or run `scASfind`, use the following codes in an R session:
```
install.packages("devtools")
devtools::install_github("hemberg-lab/scASfind")
library(scASfind)
```
## User guide
Please refer to the `scASfind` package [Vignette](https://github.com/hemberg-lab/scASfind/blob/master/vignettes/example_workflow.ipynb) for a detailed user guide.
The example scASfind index can be downloaded [here](https://drive.google.com/file/d/1dsaNOJKJN6PQQME1s0CMRPheiR94qkjv/view?usp=sharing) (too big).

## Contact
Please contact Yuyao Song (ys585@cam.ac.uk) for any enquiries or bug reports.
