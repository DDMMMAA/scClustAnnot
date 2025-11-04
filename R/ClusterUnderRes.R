#' Cluster Seurat object under given range of resolution
#'
#' Utilize Clustree, plot clustering tree showing the relationship
#' between clustering under various resolution and return clustered Seurat object
#' with clustering result stored in "metadata" slot with naming convention
#' "RNA_snn_res.<resolution>".
#'
#' @param obj a Seurat object after pre-processing, PCA, and Neighbors constructed
#' @param start The start of range of resolution (>= 0)
#' @param end The end of range of resolution, inclusive (>= start)
#' @param margin The margin of range (positive rational number)
#' @param showPlot whether display the clustering tree plot (T/F)
#'
#' @return Returns an Seurat object with clustering result stored in "metadata"
#' slot with naming convention "RNA_snn_res.<resolution>".
#'
#' @examples
#' # example using pbmc data
#' # Complete pre-processing prior clustering
#' pbmc <- Seurat::ScaleData(pbmc)
#' pbmc <- Seurat::RunPCA(pbmc)
#' pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:15, reduction = "pca")
#' # perform clustering using default parameter without plotting
#' pbmc <- ClusterUnderRes(pbmc, showPlot = FALSE)
#'
#' @references
#' Zappia, L., & Oshlack, A. (2018). Clustering trees: A visualization for evaluating clusterings at multiple resolutions. GigaScience, 7(7), giy083. \href{https://doi.org/10.1093/gigascience/giy083}{Link}
#'
#' Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology, 42(2), 293–304. \href{https://doi.org/10.1038/s41587-023-01767-y}{Link}
#'
#' @export
#' @import Seurat
#' @import clustree
#' @import ggplot2
#' @import ggraph


ClusterUnderRes <- function(obj, start = 0, end = 1.5, margin = 0.1, showPlot = TRUE) {
  # Precondition checking
  if (!inherits(obj, "Seurat") ||
      (start < 0) ||
      (end < start) ||
      (margin <= 0)) {
    message("# Argument:")
    message("#   obj: a Seurat object")
    message("#   start: the start of range of resolution (>= 0, default == 0)")
    message("#   end: the end of range of resolution (>= start, default == 1.5)")
    message("#   margin: the margin of range (positive rational number, default == 0.1)")
  } else {
    range <- seq(start, end, margin)
    # loop over range
    for(i in range) {
      obj <- FindClusters(obj, resolution = i)
      message("Added clustering under resolution: ", i)
    }
    # plot cluster tree if showPlot == T
    if (showPlot) {
      print(clustree::clustree(obj))
    }
    return(obj)
  }
}

#[END]
