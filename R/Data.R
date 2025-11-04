#' 10X Genomics Peripheral Blood Mononuclear Cells' (PBMC) scRNA-seq data
#'
#' Consist of 2638 cell's metadata and RNA expression data
#' post quality control and normalization from Seurat's vignettes
#' (https://satijalab.org/seurat/articles/pbmc3k_tutorial)
#'
#' @source 10X Genomics, Inc
#'
#' @format A Seurat object with the following slots filled:
#' \describe{
#'  \item{assays}{RNA expression data}
#'  \item{meta.data}{Cell level metadata}
#'  \item{active.ident}{current cell identify}
#'  \item{version}{Seurat verison used to create this obejct}
#'  \item{commands}{Command history}
#' }
#' @examples
#' \dontrun{
#'  pbmc
#' }
"pbmc"

