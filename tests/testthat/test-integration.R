test_that("Integration: ClusterUnderRes -> FindClusterDeg -> ClusterGo pipeline works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("clusterProfiler")

  # --- 1. Setup Data (Tiny Subset) ---
  if (!exists("pbmc", inherits = FALSE)) {
    suppressWarnings(try(data("pbmc", package = "scClustAnnot"), silent = TRUE))
  }
  skip_if_not(exists("pbmc"), "Built-in `pbmc` dataset not found")

  # Use a subset of 50 cells to make the test fast
  mini_pbmc <- pbmc[, 1:50]

  # Ensure basic Seurat slots are ready
  mini_pbmc <- suppressWarnings(Seurat::NormalizeData(mini_pbmc, verbose = FALSE))
  mini_pbmc <- suppressWarnings(Seurat::FindVariableFeatures(mini_pbmc, verbose = FALSE))
  mini_pbmc <- suppressWarnings(Seurat::ScaleData(mini_pbmc, verbose = FALSE))
  mini_pbmc <- suppressWarnings(Seurat::RunPCA(mini_pbmc, npcs = 10, verbose = FALSE))
  mini_pbmc <- suppressWarnings(Seurat::FindNeighbors(mini_pbmc, dims = 1:5, verbose = FALSE))

  # --- 2. Mock External APIs ---
  # We mock enrichGO because we don't want to load large OrgDbs during tests
  with_mocked_bindings(
    enrichGO = function(...) { return("MockGoResult") },
    {
      # STEP A: ClusterUnderRes (Generate Clusters)
      # Run for a specific resolution (e.g., 0.5)
      obj <- ClusterUnderRes(
        mini_pbmc,
        start = 0.5,
        end = 0.5,
        margin = 0.1,
        showPlot = FALSE
      )

      # Verify the column was created
      expect_true("RNA_snn_res.0.5" %in% colnames(obj@meta.data))

      # BRIDGE STEP: Set Identity
      # FindClusterDeg requires active identities.
      # We set them to the result of Step A.
      Seurat::Idents(obj) <- obj$RNA_snn_res.0.5

      # Verify we actually have clusters (at least 1)
      expect_true(length(unique(Seurat::Idents(obj))) > 0)

      # STEP B: FindClusterDeg (Generate Markers)
      # We use very loose thresholds to ensure we find *something* # even in this tiny dataset, so the pipeline continues.
      deg_res <- FindClusterDeg(
        obj,
        pct_min = 0.01,
        logfc_min = 0.01,
        pval_max = 1.0
      )

      expect_true(!is.null(deg_res$by_cluster))

      # STEP C: ClusterGo (Functional Analysis)
      # This checks if the output of Step B is compatible with Step C
      go_res <- ClusterGo(deg_res, org_db = "dummy.db")

      # Final Assertions
      expect_type(go_res, "list")

      # Check if keys match:
      # If FindClusterDeg found markers for "0", ClusterGo should have "cluster_0"
      clusters_found <- names(deg_res$by_cluster)

      # Skip if no markers were found (possible in tiny random subsets),
      # but if markers exist, GO should exist (or be NA)
      if (length(clusters_found) > 0) {
        expected_names <- paste0("cluster_", clusters_found)
        expect_true(all(expected_names %in% names(go_res)))

        # Check that our mock was called (or NA was assigned)
        # We expect "MockGoResult" or NA
        val <- go_res[[1]]
        expect_true(identical(val, "MockGoResult") || is.na(val))
      }
    }
  )
})
