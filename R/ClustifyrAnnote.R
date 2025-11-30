#' Annotate cell using Clustifyr package
#'
#' Utilize Clustifyr, annotate cell automatically at cluster level by map to
#' selected previously annotated scRNA-seq data. Annotation result will be stored
#' in metadata slot as Clustifyr.labels
#'
#' @param obj a Seurat object after pre-processing, PCA, and Neighbors constructed
#' @param ref_mat A reference expression matrix (genes as rows, cell types as columns).
#' @param cluster_col The name of the metadata column containing the cluster identities to annotate. Defaults to "seurat_clusters".
#'
#' @return Returns an Seurat object with Clustofyr annotation result added into
#' metadata slot as Clustifyr.labels
#'
#' @examples
#'
#' # Example 1 using subset of pbmc data
#' # Complete pre-processing prior clustering
#' sub_pbmc <- Seurat::ScaleData(sub_pbmc)
#' sub_pbmc <- Seurat::RunPCA(sub_pbmc)
#' sub_pbmc <- Seurat::FindNeighbors(sub_pbmc, dims = 1:15, reduction = "pca")
#' sub_pbmc <- Seurat::FindClusters(sub_pbmc, resolution = 0.5)
#' sub_pbmc <- Seurat::RunUMAP(sub_pbmc, dims = 1:15)
#' # perform annotation using Human hematopoietic cell microarray
#' # (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24759/matrix/GSE24759_series_matrix.txt.gz)
#' sub_pbmc <- ClustifyrAnnote(sub_pbmc, clustifyrdatahub::ref_hema_microarray())
#'
#' \dontrun{
#' # Example 2 using pbmc data
#' # Complete pre-processing prior clustering
#' pbmc <- Seurat::ScaleData(pbmc)
#' pbmc <- Seurat::RunPCA(pbmc)
#' pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:15, reduction = "pca")
#' pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
#' pbmc <- Seurat::RunUMAP(pbmc, dims = 1:15)
#' # perform annotation using Human hematopoietic cell microarray
#' # (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24759/matrix/GSE24759_series_matrix.txt.gz)
#' pbmc <- ClustifyrAnnote(pbmc, clustifyrdatahub::ref_hema_microarray())
#' }
#'
#' @references
#' Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology, 42(2), 293–304. \href{https://doi.org/10.1038/s41587-023-01767-y}{Link}
#'
#' Fu, R., Gillen, A. E., Sheridan, R. M., Tian, C., Daya, M., Hao, Y., Hesselberth, J. R., & Riemondy, K. A. (2020). clustifyr: An R package for automated single-cell RNA sequencing cluster classification. F1000Research. \href{https://doi.org/10.12688/f1000research.22969.2}{Link}
#'
#' @export
#' @import Seurat
#' @import clustifyr
#' @import clustifyrdatahub

ClustifyrAnnote <- function(
    obj,
    ref_mat,
    cluster_col = "seurat_clusters"
) {
  # Parameter checking
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.", call. = FALSE)
  }

  # Check if cluster_col exists in metadata
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop(paste0("Column '", cluster_col, "' not found in object metadata. Run FindClusters first."), call. = FALSE)
  }

  # Check ref_mat
  if (missing(ref_mat) || !is.matrix(ref_mat) && !is.data.frame(ref_mat)) {
    stop("`ref_mat` must be a matrix or data frame of reference expression data.", call. = FALSE)
  }

  # Function logic

  # 1. Run clustifyr
  # We force obj_out = FALSE so we get the correlation matrix needed for cor_to_call
  res_mat <- clustifyr::clustify(
    input = obj,
    cluster_col = cluster_col,
    ref_mat = ref_mat,
    obj_out = FALSE
  )

  # 2. Call the cell types (Best fit per cluster)
  clusterIdent <- clustifyr::cor_to_call(res_mat)

  # Ensure the results are sorted by cluster ID for clean mapping
  clusterIdent <- clusterIdent[order(clusterIdent$cluster), ]

  # 3. Create a lookup vector (The dictionary)
  # Structure: names = Cluster IDs, values = Cell Type Labels
  lookup_table <- stats::setNames(
    as.character(clusterIdent$type),    # Values: Cell Types
    as.character(clusterIdent$cluster)  # Keys: Cluster IDs
  )

  # 4. Apply the mapping
  # Assign to "Clustifyr.labels" slot
  obj[["Clustifyr.labels"]] <- factor(unname(lookup_table[as.character(obj[[cluster_col]][,1])]))

  return(obj)
}

#[END]
