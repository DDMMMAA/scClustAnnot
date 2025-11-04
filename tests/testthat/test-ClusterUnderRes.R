test_that("ClusterUnderRes handles invalid inputs properly", {

  # Invalid: not a Seurat object
  expect_message(
    ClusterUnderRes("not_a_seurat"),
    "# Argument:"
  )

  # Invalid: negative start
  expect_message(
    ClusterUnderRes(pbmc, start = -1),
    "# Argument:"
  )

  # Invalid: margin <= 0
  expect_message(
    ClusterUnderRes(pbmc, margin = 0),
    "# Argument:"
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
