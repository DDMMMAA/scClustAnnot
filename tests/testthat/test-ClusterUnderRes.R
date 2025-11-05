test_that("ClusterUnderRes validates inputs with errors", {
  # not a Seurat object
  expect_error(
    ClusterUnderRes("not_a_seurat"),
    "`obj` must be a Seurat object\\."
  )

  # load your built-in pbmc if needed
  if (!exists("pbmc", inherits = FALSE)) {
    suppressWarnings(try(data("pbmc", package = "scClustAnnot"), silent = TRUE))
  }
  skip_if_not(exists("pbmc"), "Built-in `pbmc` not found")

  # start must be >= 0
  expect_error(
    ClusterUnderRes(pbmc, start = -0.1),
    "`start` must be a single finite numeric value >= 0\\."
  )

  # end must be >= start
  expect_error(
    ClusterUnderRes(pbmc, start = 0.5, end = 0.4),
    "`end` must be a single finite numeric value >= `start`\\."
  )

  # margin must be > 0
  expect_error(
    ClusterUnderRes(pbmc, margin = 0),
    "`margin` must be a single finite numeric value > 0\\."
  )

  # showPlot must be single logical
  expect_error(
    ClusterUnderRes(pbmc, showPlot = "TRUE"),
    "`showPlot` must be a single logical \\(TRUE/FALSE\\)\\."
  )
})

test_that("ClusterUnderRes runs and adds clustering metadata", {

  # pre-process
  pbmc <- ScaleData(pbmc)
  pbmc <- RunPCA(pbmc)
  pbmc <- FindNeighbors(pbmc, dims = 1:15, reduction = "pca")

  obj <- ClusterUnderRes(pbmc, start = 0, end = 0.2, margin = 0.2, showPlot = FALSE)

  # Ensure return is a Seurat object
  expect_s4_class(obj, "Seurat")

  # Ensure metadata columns were added
  meta_names <- colnames(obj@meta.data)
  expect_true(any(grepl("RNA_snn_res", meta_names)))
})
