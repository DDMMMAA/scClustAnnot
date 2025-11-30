test_that("ClustifyrAnnote validates inputs with errors", {
  # 1. Load data if needed
  if (!exists("pbmc", inherits = FALSE)) {
    suppressWarnings(try(data("pbmc", package = "scClustAnnot"), silent = TRUE))
  }
  skip_if_not(exists("pbmc"), "Built-in `pbmc` dataset not found")

  # Ensure pbmc has a dummy cluster column for validation tests
  # (In case the raw pbmc doesn't have it yet)
  pbmc$seurat_clusters <- factor(rep(0, ncol(pbmc)))

  # 2. Test: Wrong object type
  expect_error(
    ClustifyrAnnote("not_a_seurat", ref_mat = matrix()),
    "`obj` must be a Seurat object\\."
  )

  # 3. Test: Missing cluster column
  # We test this by asking for a column that definitely doesn't exist
  expect_error(
    ClustifyrAnnote(pbmc, ref_mat = matrix(), cluster_col = "non_existent_col"),
    "Column 'non_existent_col' not found in object metadata"
  )

  # 4. Test: Missing or invalid ref_mat
  expect_error(
    ClustifyrAnnote(pbmc, cluster_col = "seurat_clusters"), # Missing ref_mat
    "`ref_mat` must be a matrix or data frame"
  )

  expect_error(
    ClustifyrAnnote(pbmc, ref_mat = "not_a_matrix"), # Wrong type
    "`ref_mat` must be a matrix or data frame"
  )
})

