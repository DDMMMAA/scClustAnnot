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
#'
#' # Example 1 using subset of pbmc data
#' # Complete pre-processing prior clustering
#' sub_pbmc <- Seurat::ScaleData(sub_pbmc)
#' sub_pbmc <- Seurat::RunPCA(sub_pbmc)
#' sub_pbmc <- Seurat::FindNeighbors(sub_pbmc, dims = 1:15, reduction = "pca")
#' # perform clustering using default parameter without plotting
#' sub_pbmc <- ClusterUnderRes(sub_pbmc, showPlot = FALSE)
#'
#' \dontrun{
#' # Example 2 using pbmc data
#' # Complete pre-processing prior clustering
#' pbmc <- Seurat::ScaleData(pbmc)
#' pbmc <- Seurat::RunPCA(pbmc)
#' pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:15, reduction = "pca")
#' # perform clustering using default parameter without plotting
#' pbmc <- ClusterUnderRes(pbmc, showPlot = FALSE)
#' }
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


ClusterUnderRes <- function(
    obj,
    start = 0,
    end = 1.5,
    margin = 0.1,
    showPlot = TRUE
    ) {
  # Arguments checking
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.", call. = FALSE)
  }

  # start
  if (!is.numeric(start) || length(start) != 1L || !is.finite(start) || start < 0) {
    stop("`start` must be a single finite numeric value >= 0.", call. = FALSE)
  }

  # end
  if (!is.numeric(end) || length(end) != 1L || !is.finite(end) || end < start) {
    stop("`end` must be a single finite numeric value >= `start`.", call. = FALSE)
  }

  # margin
  if (!is.numeric(margin) || length(margin) != 1L || !is.finite(margin) || margin <= 0) {
    stop("`margin` must be a single finite numeric value > 0.", call. = FALSE)
  }

  # showPlot
  if (!is.logical(showPlot) || length(showPlot) != 1L || is.na(showPlot)) {
    stop("`showPlot` must be a single logical (TRUE/FALSE).", call. = FALSE)
  }

  # Function logic
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

#[END]
