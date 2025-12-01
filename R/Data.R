#' 10X Genomics Peripheral Blood Mononuclear Cells' (PBMC) scRNA-seq data
#'
#' Consist of 2638 cell's metadata and RNA expression data
#' post quality control from Seurat's vignettes
#' (https://satijalab.org/seurat/articles/pbmc3k_tutorial)
#'
#' @source 10X Genomics, Inc
#'
#' @format A Seurat object with the following slots filled:
#' \describe{
#'   \item{assays}{RNA expression data (contains 'counts' and 'data' layers)}
#'   \item{meta.data}{Basic cell-level metadata (e.g., nFeature_RNA, nCount_RNA)}
#'   \item{active.ident}{Default identity (typically set to project name, as no clusters are defined)}
#'   \item{version}{Seurat version used to create this object}
#'   \item{commands}{Command history (subsetting and pre-processing steps)}
#' }
#' @examples
#' \dontrun{
#'  pbmc
#' }
"pbmc"

#' Random Subset of 10X Genomics Peripheral Blood Mononuclear Cells (PBMC) scRNA-seq data
#'
#' Consist of 500 cell's metadata and RNA expression data
#' post quality control from Seurat's vignettes
#' (https://satijalab.org/seurat/articles/pbmc3k_tutorial)
#'
#' @source 10X Genomics, Inc
#'
#' @format A Seurat object with the following slots filled:
#' \describe{
#'   \item{assays}{RNA expression data (contains 'counts' and 'data' layers)}
#'   \item{meta.data}{Basic cell-level metadata (e.g., nFeature_RNA, nCount_RNA)}
#'   \item{active.ident}{Default identity (typically set to project name, as no clusters are defined)}
#'   \item{version}{Seurat version used to create this object}
#'   \item{commands}{Command history (subsetting and pre-processing steps)}
#' }
#' @usage sub_pbmc
#' @examples
#' \dontrun{
#'  sub_pbmc
#' }
"sub_pbmc"
