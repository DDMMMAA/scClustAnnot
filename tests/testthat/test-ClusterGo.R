test_that("ClusterGo validates inputs", {
  # 1. Test missing 'by_cluster' slot
  bad_input <- list(wrong_slot = "data")
  expect_error(
    ClusterGo(bad_input),
    "Input `Deg_result` must contain a `by_cluster` slot"
  )
})

test_that("ClusterGo handles different enrichment scenarios (Mocked)", {
  skip_if_not_installed("clusterProfiler")

  # --- Setup Dummy Input Data ---
  # We create 3 scenarios:
  # Cluster "0": Normal genes (Should succeed)
  # Cluster "1": Empty gene list (Should return NA immediately)
  # Cluster "2": Genes that produce no enrichment (Should return NA from NULL result)

  deg_data <- list(
    by_cluster = list(
      "0" = list(gene = c("GeneA", "GeneB")),
      "1" = list(gene = character(0)),
      "2" = list(gene = c("NoEnrichGene"))
    )
  )

  # --- Mocking ---
  # We mock enrichGO to avoid loading the huge OrgDb and to force specific outcomes
  with_mocked_bindings(
    enrichGO = function(gene, ...) {
      # Scenario: Genes that result in NO enrichment
      if ("NoEnrichGene" %in% gene) {
        return(NULL)
      }

      # Scenario: Success
      # We return a dummy S4-like object or list.
      # Since the function just stores it, a string is sufficient for testing structure.
      return("SuccessObject")
    },

    # --- Test Execution ---
    {
      # We pass a dummy string for org_db because the mock ignores it anyway
      res <- ClusterGo(deg_data, org_db = "dummy.db")

      # --- Assertions ---

      # 1. Check Output Structure
      expect_type(res, "list")
      expect_named(res, c("cluster_0", "cluster_1", "cluster_2"))

      # 2. Check Scenario "0" (Success)
      expect_equal(res[["cluster_0"]], "SuccessObject")

      # 3. Check Scenario "1" (Empty Input)
      # Logic: length(gene) == 0 -> assigns NA directly
      expect_true(is.na(res[["cluster_1"]]))

      # 4. Check Scenario "2" (No Enrichment found)
      # Logic: enrichGO returned NULL -> function assigns NA
      expect_true(is.na(res[["cluster_2"]]))
    }
  )
})

test_that("ClusterGo handles errors gracefully", {
  skip_if_not_installed("clusterProfiler")

  deg_data <- list(by_cluster = list("0" = list(gene = c("GeneA"))))

  # Mock enrichGO to throw an error
  with_mocked_bindings(
    enrichGO = function(...) {
      stop("Database connection failed")
    },
    {
      # The function uses tryCatch, so it should NOT crash, but return NA/NULL
      # (Depending on how your tryCatch is written. Your code returns NULL on error,
      # then checks if(is.null), so it becomes NA).

      expect_warning(
        res <- ClusterGo(deg_data, org_db = "dummy.db"),
        "Error in Cluster 0"
      )

      expect_true(is.na(res[["cluster_0"]]))
    }
  )
})
