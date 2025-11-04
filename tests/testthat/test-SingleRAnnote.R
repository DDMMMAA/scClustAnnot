test_that("SingleRAnnote validates inputs and prints helpful messages", {
  # Not a Seurat object
  expect_message(
    SingleRAnnote("not_a_seurat", "dice", "2024-02-26"),
    "# Argument:"
  )

  # invalid name (not character)
  expect_message(
    SingleRAnnote(pbmc, name = 123, version = "2024-02-26"),
    "# Argument:"
  )

  # invalid version (not character)
  expect_message(
    SingleRAnnote(pbmc, name = "dice", version = 20240226),
    "# Argument:"
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
