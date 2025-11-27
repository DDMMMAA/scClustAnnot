#' Annotate cell using SingleR package
#'
#' Utilize SingleR, annotate cell automatically at single cell level by map to
#' selected previously annotated scRNA-seq data. Annotation result will be stored
#' in metadata slot as SingleR.labels
#'
#' @param obj a Seurat object after pre-processing, PCA, and Neighbors constructed
#' @param name the name of reference dataset. See available dataset \href{https://www.bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html}{here}
#' @param version the version of reference dataset. See available dataset See available dataset \href{https://www.bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html}{here}
#'
#' @return Returns an Seurat object with SingleR annotation result added into
#'metadata slot as SingleR.labels
#'
#' @examples
#'
#' # Example 1 using subset of pbmc data
#' # Complete pre-processing prior clustering
#' sub_pbmc <- Seurat::ScaleData(sub_pbmc)
#' sub_pbmc <- Seurat::RunPCA(sub_pbmc)
#' sub_pbmc <- Seurat::FindNeighbors(sub_pbmc, dims = 1:15, reduction = "pca")
#' # perform annotation using bulk RNA-seq samples from https://doi.org/10.1016/j.cell.2018.10.022
#' sub_pbmc <- SingleRAnnote(sub_pbmc, "dice", "2024-02-26")
#' # PCA visualization
#' Seurat::DimPlot(sub_pbmc, reduction = "pca", label = TRUE, pt.size = 0.5,
#' group.by = "SingleR.labels", repel = TRUE)
#'
#' \dontrun{
#' # Example 2 using pbmc data
#' # Complete pre-processing prior clustering
#' pbmc <- Seurat::ScaleData(pbmc)
#' pbmc <- Seurat::RunPCA(pbmc)
#' pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:15, reduction = "pca")
#' # perform annotation using bulk RNA-seq samples from https://doi.org/10.1016/j.cell.2018.10.022
#' pbmc <- SingleRAnnote(pbmc, "dice", "2024-02-26")
#' # PCA visualization
#' Seurat::DimPlot(pbmc, reduction = "pca", label = TRUE, pt.size = 0.5,
#' group.by = "SingleR.labels", repel = TRUE)
#' }
#'
#' @references
#' Schmiedel, B. J., Singh, D., Madrigal, A., Valdovino-Gonzalez, A. G., White, B. M., Zapardiel-Gonzalo, J., Ha, B., Altay, G., Greenbaum, J. A., McVicker, G., Seumois, G., Rao, A., Kronenberg, M., Peters, B., & Vijayanand, P. (2018). Impact of Genetic Polymorphisms on Human Immune Cell Gene Expression. Cell, 175(6), 1701-1715.e16. \href{https://doi.org/10.1016/j.cell.2018.10.022}{Link}
#'
#' Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology, 42(2), 293–304. \href{https://doi.org/10.1038/s41587-023-01767-y}{Link}
#'
#' Aran, D., Looney, A. P., Liu, L., Wu, E., Fong, V., Hsu, A., Chak, S., Naikawadi, R. P., Wolters, P. J., Abate, A. R., Butte, A. J., & Bhattacharya, M. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology, 20(2), 163–172. \href{https://doi.org/10.1038/s41590-018-0276-y}{Link}
#'
#' @export
#' @import Seurat
#' @import SingleR
#' @importFrom celldex fetchReference
#' @import ggplot2
#' @import ggraph


SingleRAnnote <- function(obj, name, version) {
  # Parameter checking
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.", call. = FALSE)
  }

  # name: single non-empty character
  if (!is.character(name) || length(name) != 1L || is.na(name) || name == "") {
    stop("`name` must be a non-empty single character string (e.g., 'dice').", call. = FALSE)
  }

  # version: single non-empty character
  if (!is.character(version) || length(version) != 1L || is.na(version) || version == "") {
    stop("`version` must be a non-empty single character string (e.g., '2024-02-26').", call. = FALSE)
  }

  # Function logic
  # fetch reference scRNA-seq data via celldex package
  ref_data <- celldex::fetchReference(name, version)
  SingleR_result <- SingleR(as.data.frame(as.matrix(obj[["RNA"]]$data)),
                            ref_data, ref_data$label.main)
  obj$SingleR.labels <- SingleR_result$labels
  return(obj)
}

#[END]
