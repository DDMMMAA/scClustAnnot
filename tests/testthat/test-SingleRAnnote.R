test_that("SingleRAnnote validates inputs with errors", {
  # wrong obj
  expect_error(
    SingleRAnnote("not_a_seurat", "dice", "2024-02-26"),
    "`obj` must be a Seurat object\\."
  )
  # wrong name
  if (!exists("pbmc", inherits = FALSE)) {
    suppressWarnings(try(data("pbmc", package = "scClustAnnot"), silent = TRUE))
  }
  skip_if_not(exists("pbmc"), "Built-in `pbmc` dataset not found")

  expect_error(
    SingleRAnnote(pbmc, name = NA_character_, version = "2024-02-26"),
    "`name` must be a non-empty single character string"
  )
  expect_error(
    SingleRAnnote(pbmc, name = "dice", version = ""),
    "`version` must be a non-empty single character string"
  )
})


test_that("SingleRAnnote adds SingleR.labels to metadata (mocked for speed)", {
  skip_if_not_installed("Seurat")

  # mock out SingleR and fetchReference to avoid heavy data download
  with_mocked_bindings(
    fetchReference = function(name, version) {
      list(dummy = TRUE)
    },
    SingleR = function(test, ref, labels) {
      n_cells <- ncol(as.matrix(test))
      list(labels = rep("MockType", n_cells))
    },
    {
      obj2 <- SingleRAnnote(pbmc, name = "dice", version = "2024-02-26")

      expect_s4_class(obj2, "Seurat")
      expect_true("SingleR.labels" %in% colnames(obj2@meta.data))
      expect_equal(length(obj2$SingleR.labels), ncol(obj2))
      expect_true(all(obj2$SingleR.labels == "MockType"))
    }
  )
})
