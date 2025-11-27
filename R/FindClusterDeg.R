#' Find, filter, and score cluster markers (Differential experssed genes)
#'
#' Runs Seurat::FindAllMarkers(), filters markers by prevalence, p-value,
#' and log2FC thresholds, then computes a p_FC score = (1 - p_val) * avg_log2FC.
#' Return a list of length two compose of splitted and combined DEG result
#'
#' @param obj Seurat object with identities set (Idents(obj)).
#' @param only_pos logical; pass to FindAllMarkers (default TRUE).
#' @param pct_min Minimum fraction expressed in each group (applied to both pct.1 and pct.2).
#' @param pval_max Maximum raw p-value to keep.
#' @param logfc_min Keep rows with |avg_log2FC| >= logfc_min.
#'
#' @return A list with:
#' \describe{
#'   \item{by_cluster}{named list of filtered marker data.frames per cluster}
#'   \item{combined}{single data.frame of filtered markers with p_FC}
#' }
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
#' # Set clustering result under resolution = 0.5
#' sub_pbmc <- Seurat::FindClusters(sub_pbmc, resolution = 0.5)
#' # Find Cluster DEGs
#' result <- FindClusterDeg(sub_pbmc)
#'
#' \dontrun{
#' # Example 2 using pbmc data
#' # Complete pre-processing prior clustering
#' pbmc <- Seurat::ScaleData(pbmc)
#' pbmc <- Seurat::RunPCA(pbmc)
#' pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:15, reduction = "pca")
#' # perform clustering using default parameter without plotting
#' pbmc <- ClusterUnderRes(pbmc, showPlot = FALSE)
#' # Set clustering result under resolution = 0.5
#' pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
#' # Find Cluster DEGs
#' result <- FindClusterDeg(pbmc)
#' }
#'
#' @references
#' Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology, 42(2), 293–304. \href{https://doi.org/10.1038/s41587-023-01767-y}{Link}
#'
#' @export
#' @import Seurat
#' @import clustree
#' @import ggplot2
#' @import ggraph


FindClusterDeg <- function(
    obj,
    only_pos = TRUE,
    pct_min = 0.1,
    pval_max = 0.05,
    logfc_min = 0.25
  ) {
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.")
  }

  # Check logical only_pos
  if (!is.logical(only_pos) || length(only_pos) != 1 || is.na(only_pos)) {
    stop("`only_pos` must be a single logical value (TRUE/FALSE).", call. = FALSE)
  }

  # Check numeric arguments
  if (!is.numeric(pct_min) || length(pct_min) != 1 || is.na(pct_min) ||
      pct_min < 0 || pct_min > 1) {
    stop("`pct_min` must be a single numeric value between 0 and 1.", call. = FALSE)
  }

  if (!is.numeric(pval_max) || length(pval_max) != 1 || is.na(pval_max) ||
      pval_max <= 0 || pval_max > 1) {
    stop("`pval_max` must be a single numeric value in (0, 1].", call. = FALSE)
  }

  if (!is.numeric(logfc_min) || length(logfc_min) != 1 || is.na(logfc_min)) {
    stop("`logfc_min` must be a single numeric value.", call. = FALSE)
  }

  # Find markers
  markers <- Seurat::FindAllMarkers(obj, only.pos = only_pos)

  # Split by cluster
  markers_split <- split(markers, markers$cluster)

  # Filter per-cluster
  filtered_list <- lapply(markers_split, function(df) {
    keep <- !is.na(df$p_val) &
      df$p_val < pval_max &
      !is.na(df$pct.1) & df$pct.1 >= pct_min &
      !is.na(df$pct.2) & df$pct.2 >= pct_min &
      !is.na(df$avg_log2FC) & abs(df$avg_log2FC) >= logfc_min
    out <- df[keep, , drop = FALSE]

    # Compute p_FC
    if (nrow(out) > 0) {
      out$p_FC <- (1 - out$p_val) * out$avg_log2FC
    } else {
      out$p_FC <- numeric(0)
    }
    out
  })

  # Combined data.frame (sorted by cluster then descending p_FC)
  combined <- do.call(rbind, filtered_list)
  if (!is.null(combined) && nrow(combined) > 0) {
    ord <- order(combined$cluster, -combined$p_FC)
    combined <- combined[ord, , drop = FALSE]
  }

  list(
    by_cluster = filtered_list,
    combined = combined
  )
}
