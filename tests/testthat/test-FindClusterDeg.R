test_that("FindClusterDeg validates inputs", {
  # Not a Seurat object
  expect_error(
    FindClusterDeg("not_seurat"),
    "`obj` must be a Seurat object\\."
  )

  # Load built-in pbmc from this package (ensure available in test env)
  if (!exists("pbmc", inherits = FALSE)) {
    suppressWarnings(try(data("pbmc", package = "scClustAnnot"), silent = TRUE))
  }
  skip_if_not(exists("pbmc"), "built-in `pbmc` dataset not found")

  # only_pos must be single logical
  expect_error(
    FindClusterDeg(pbmc, only_pos = "TRUE"),
    "`only_pos` must be a single logical value"
  )
  expect_error(
    FindClusterDeg(pbmc, only_pos = c(TRUE, FALSE)),
    "`only_pos` must be a single logical value"
  )

  # pct_min in [0,1]
  expect_error(
    FindClusterDeg(pbmc, pct_min = -0.1),
    "`pct_min` must be a single numeric value between 0 and 1"
  )
  expect_error(
    FindClusterDeg(pbmc, pct_min = 1.1),
    "`pct_min` must be a single numeric value between 0 and 1"
  )

  # pval_max in (0,1]
  expect_error(
    FindClusterDeg(pbmc, pval_max = 0),
    "`pval_max` must be a single numeric value in \\(0, 1\\]"
  )
  expect_error(
    FindClusterDeg(pbmc, pval_max = 2),
    "`pval_max` must be a single numeric value in \\(0, 1\\]"
  )

  # logfc_min numeric scalar
  expect_error(
    FindClusterDeg(pbmc, logfc_min = NA_real_),
    "`logfc_min` must be a single numeric value"
  )
  expect_error(
    FindClusterDeg(pbmc, logfc_min = c(0.25, 0.5)),
    "`logfc_min` must be a single numeric value"
  )
})

