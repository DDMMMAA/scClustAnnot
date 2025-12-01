#' Perform Gene Ontology (GO) Enrichment Analysis for Clusters
#'
#' Iterates through a list of cluster markers (typically the output of
#' FindClusterDeg) and performs GO enrichment analysis using clusterProfiler.
#'
#' @param Deg_result A list object typically returned by `FindClusterDeg`.
#' It must contain a sub-list named `by_cluster`, where each element contains
#' a `gene` slot.
#' @param org_db The organism database object or name.
#' @param ont One of "BP", "MF", or "CC".
#' @param p_val_cutoff P-value cutoff for enrichment.
#' @param q_val_cutoff Q-value cutoff for enrichment.
#'
#' @return A named list of enrichment results (`enrichResult` objects).
#' Clusters with no significant enrichment are stored as NA.
#'
#' @examples
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
#' DegResult <- FindClusterDeg(sub_pbmc)
#' # Perform GO enrichment analysis on Cluster DEGs
#' GoResult <- ClusterGo(DegResult)
#'
#' \dontrun{
#' # Example 2 using pbmc data
#' # Complete pre-processing prior clustering
#' sub_pbmc <- Seurat::ScaleData(sub_pbmc)
#' sub_pbmc <- Seurat::RunPCA(sub_pbmc)
#' sub_pbmc <- Seurat::FindNeighbors(sub_pbmc, dims = 1:15, reduction = "pca")
#' # perform clustering using default parameter without plotting
#' sub_pbmc <- ClusterUnderRes(sub_pbmc, showPlot = FALSE)
#' # Set clustering result under resolution = 0.5
#' sub_pbmc <- Seurat::FindClusters(sub_pbmc, resolution = 0.5)
#' # Find Cluster DEGs
#' DegResult <- FindClusterDeg(sub_pbmc)
#' # Perform GO enrichment analysis on Cluster DEGs
#' GoResult <- ClusterGo(DegResult)
#' }
#'
#' @references
#' Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation, 2(3). \href{https://doi.org/10.1016/j.xinn.2021.100141}{Link}
#'
#' @export
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db

ClusterGo <- function(
    Deg_result,
    org_db = "org.Hs.eg.db",
    ont = "BP",
    p_val_cutoff = 0.05,
    q_val_cutoff = 0.02
) {

  # 1. Validate Input
  if (is.null(Deg_result$by_cluster)) {
    stop("Input `Deg_result` must contain a `by_cluster` slot.", call. = FALSE)
  }

  # Initialize output list
  go_results <- list()

  # 2. Iterate through clusters using names
  for (cluster_id in names(Deg_result$by_cluster)) {

    result_name <- paste("cluster_", cluster_id, sep = "")

    # Extract gene list using the Deg_result structure
    gene_list <- Deg_result$by_cluster[[cluster_id]]$gene

    # Skip empty gene lists to avoid errors
    if (length(gene_list) == 0) {
      message("Skipping Cluster ", cluster_id, ": No genes found.")
      go_results[[result_name]] <- NA
      next
    }

    message("Running GO analysis for Cluster: ", cluster_id, " (", length(gene_list), " genes)")

    # 3. Run enrichGO with error handling
    enrich_res <- tryCatch({
      enrichGO(
        gene          = gene_list,
        keyType       = "SYMBOL",
        OrgDb         = org_db,
        ont           = ont,
        pvalueCutoff  = p_val_cutoff,
        qvalueCutoff  = q_val_cutoff,
        readable      = FALSE
      )
    }, error = function(e) {
      warning("Error in Cluster ", cluster_id, ": ", e$message)
      return(NULL)
    })

    # 4. Store Result
    # If NULL (no terms found), store NA to preserve list structure
    if (is.null(enrich_res)) {
      go_results[[result_name]] <- NA
      message("  -> No significant enrichment found.")
    } else {
      go_results[[result_name]] <- enrich_res
      message("  -> Success.")
    }
  }

  return(go_results)
}

#[END]
