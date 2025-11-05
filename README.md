
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scClustAnnot

<!-- badges: start -->

<!-- badges: end -->

An R package integrate into [Seurat
V5](https://github.com/satijalab/seurat), meant to aid clustering and
annotation step in scRNA-seq analysis

## Description

`scClustAnnot` is an R package to perform clustering and annotation step
in scRNA-seq work flow. Currently, clustering and annotation in
scRNA-seq analysis requires many repetitive manual works, including
resolution tuning, identifying differentially expressed genes (DEGs),
and interpreting DEGs biological functions at given context, which is
time consuming and subject to researcher’s bias. There exists some
automatic annotation tool, however, mostly done by map to previously
annotated data. In that case, the quality of previously annotated data
is a great concern. This package aims to automate above manual task to
an extent, so that researchers’ bias can be reduced and analysis can be
performed more effectively. The `scClustAnnot` package was developed
using `R version 4.5.0 (2025-04-11 ucrt)`, platform
`Microsoft Windows 11 25H2 (64 bit)`

## Installation

You can install the development version of scClustAnnot like so:

``` r
> install.packages("devtools")
> library("devtools")
> devtools::install_github("DDMMMAA/scClustAnnot", build_vignettes = TRUE)
> library("scClustAnnot")
```

## Overview

``` r
ls("package:scClustAnnot")
data(package = "scClustAnnot") 
browseVignettes("scClustAnnot")
```

## Contributions

This package is an individual project developed and maintained solely by
Jiaqi Ma, contributions include project planing, literature review,
documentation, implementation, etc. The package is heavily based on the
`Seurat` framework. In `ClusterUnderRes()` function `Clustree` were
utilized to generate basic cluster tree. In `SingleRAnnotate()`
function, reference scRNA-seq annotation was accessed via `Celldex` API,
`SingleR` were utilized for mapping query scRNA-seq data to reference
data. In `FindClusterDeg()` funcion, `Seurat::FindAllMarkers()` play the
role of identify cluster specific markers. Development of this package
was assisted by ChatGPT 5. All generated code and documentation were
manually reviewed, refined, and tested by the author before inclusion in
the final package.

## References

…

## Acknowledgements

This package was developed as part of an assessment for 2025 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `scClustAnnot` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the GitHub issues.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scClustAnnot)
#> Warning: replacing previous import 'SingleR::NovershternHematopoieticData' by
#> 'celldex::NovershternHematopoieticData' when loading 'scClustAnnot'
#> Warning: replacing previous import 'SingleR::DatabaseImmuneCellExpressionData'
#> by 'celldex::DatabaseImmuneCellExpressionData' when loading 'scClustAnnot'
#> Warning: replacing previous import 'SingleR::MouseRNAseqData' by
#> 'celldex::MouseRNAseqData' when loading 'scClustAnnot'
#> Warning: replacing previous import 'SingleR::MonacoImmuneData' by
#> 'celldex::MonacoImmuneData' when loading 'scClustAnnot'
#> Warning: replacing previous import 'SingleR::HumanPrimaryCellAtlasData' by
#> 'celldex::HumanPrimaryCellAtlasData' when loading 'scClustAnnot'
#> Warning: replacing previous import 'SingleR::BlueprintEncodeData' by
#> 'celldex::BlueprintEncodeData' when loading 'scClustAnnot'
#> Warning: replacing previous import 'SingleR::ImmGenData' by
#> 'celldex::ImmGenData' when loading 'scClustAnnot'
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
