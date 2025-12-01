library(shiny)
library(Seurat)
library(ggplot2)
library(enrichplot)
library(scClustAnnot)

# ==============================================================================
# UI DEFINITION
# ==============================================================================
ui <- fluidPage(

  titlePanel("scClustAnnot Dashboard"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # --- 1. LOAD DATA ---
      h4("1. Load Data"),
      helpText("Step 1: Start here. Load the built-in 'sub_pbmc' for a quick demo, or upload your own Seurat object (.rda)."),
      p("Note: The uploaded .rda file must contain exactly one Seurat object. Workflow based on Seurat V5 (Hao et al., 2024).",
        style = "color: #555; margin-bottom: 10px; font-style: italic;"),
      selectInput("demo_select", "Choose Built-in Data:",
                  choices = c("sub_pbmc", "pbmc"),
                  selected = "sub_pbmc"),
      actionButton("load_demo", "Load Selected Data", icon = icon("play")),
      fileInput("upload_rda", "Or Upload .rda", accept = ".rda"),
      verbatimTextOutput("data_status"),
      hr(),

      # --- 2. PRE-PROCESSING ---
      h4("2. Pre-processing"),
      helpText("Step 2: If your data is raw (counts only), run this standard pipeline: Normalize -> Find Features -> Scale -> PCA -> UMAP."),
      p(strong("Pre-processing:"), "Normalizes gene counts to correct for sequencing depth and scales data so highly expressed genes don't dominate.", br(),
        strong("Dimension Reduction:"), "PCA compresses the dataset to capture the most biological variance. UMAP projects this into 2D space for visualization (Hao et al., 2024).",
        style = "color: #555; margin-bottom: 10px;"),
      actionButton("run_preprocess", "Run Standard Workflow", class = "btn-primary"),
      hr(),

      # --- 3. CLUSTERING ---
      h4("3. Clustering"),
      helpText("Step 3: Analyze cluster stability. 1) Define a range. 2) Run the sweep. 3) Check the 'Clustering Tree' tab. 4) Select the best resolution below."),
      p("This groups cells with similar transcriptomes. Look for stable branches in the Clustering Tree (Zappia & Oshlack, 2018). A resolution where branches do not rapidly split or cross is often biologically robust.",
        style = "color: #555; margin-bottom: 10px;"),
      sliderInput("res_range", "Resolution Range:", min = 0.1, max = 1.5, value = c(0.1, 0.8), step = 0.1),
      numericInput("res_margin", "Step Size:", value = 0.1, step = 0.1),
      actionButton("run_cluster", "Run Clustering"),

      # Dynamic UI: Resolution Dropdown (Appears after clustering)
      br(), br(),
      uiOutput("cluster_res_ui"),
      hr(),

      # --- 4. MARKERS ---
      h4("4. Find Markers"),
      helpText("Step 4: Identify marker genes for your chosen resolution. These are ranked by the 'p_FC' score: (1 - p_val) * avg_log2FC."),
      p("This finds genes uniquely expressed in each cluster. High p_FC genes are both statistically significant and strongly upregulated, defining the cluster's identity.",
        style = "color: #555; margin-bottom: 10px;"),
      numericInput("deg_fc", "Min LogFC:", value = 0.25, step = 0.05),
      numericInput("deg_pval", "Max P-Value:", value = 0.05, step = 0.01),
      actionButton("run_deg", "Find DEGs"),

      # Dynamic UI: Cluster Selector for DEGs
      br(), br(),
      uiOutput("deg_cluster_ui"),
      hr(),

      # --- 5. ANNOTATION ---
      h4("5. Annotation"),
      helpText("Step 5: Predict cell types using a reference. 'SingleR' is slower but annotates every single cell. 'Clustifyr' is faster and annotates whole clusters."),
      p("This assigns biological labels (e.g., 'T-cell') by comparing your data to a reference atlas. Uses SingleR (Aran et al., 2019) or Clustifyr (Fu et al., 2020). Use this as a prediction to cross-reference with the Marker Genes found in Step 4.",
        style = "color: #555; margin-bottom: 10px;"),
      selectInput("anno_method", "Method:", choices = c("SingleR", "Clustifyr")),

      # SingleR Options
      conditionalPanel(
        condition = "input.anno_method == 'SingleR'",
        textInput("singler_ref", "Ref Name (e.g., dice):", value = "dice"),
        helpText(a("Click here for available references (celldex)",
                   href="https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html",
                   target="_blank")),
        textInput("singler_ver", "Ref Version:", value = "2024-02-26")
      ),

      # Clustifyr Options
      conditionalPanel(
        condition = "input.anno_method == 'Clustifyr'",
        selectInput("clustifyr_ref", "Reference (clustifyrdatahub):", choices = NULL),
        helpText(a("Click here for available references (clustifyrdatahub)",
                   href="https://bioconductor.org/packages/release/data/experiment/vignettes/clustifyrdatahub/inst/doc/clustifyrdatahub.html",
                   target="_blank"))
      ),

      actionButton("run_anno", "Annotate"),
      hr(),

      # --- 6. GO ENRICHMENT ---
      h4("6. GO Analysis"),
      helpText("Step 6: Discover biological pathways. This runs Gene Ontology enrichment on the markers found in Step 4. Select the correct organism."),
      p("This translates gene lists into biological functions (e.g., 'Viral Response') using clusterProfiler (Wu et al., 2021). Dots with high GeneRatio and red color (low p-adjust) represent the primary functions of the cluster.",
        style = "color: #555; margin-bottom: 10px;"),
      selectInput("go_org", "Organism:",
                  choices = c("Human (org.Hs.eg.db)" = "org.Hs.eg.db",
                              "Mouse (org.Mm.eg.db)" = "org.Mm.eg.db")),
      selectInput("go_ont", "Ontology:", choices = c("BP" = "BP", "MF" = "MF", "CC" = "CC")),
      actionButton("run_go", "Run Enrichment")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",

        # --- TAB: About & Help ---
        tabPanel("About",
                 h3("Welcome to scClustAnnot"),
                 p("This package streamlines the clustering and annotation steps in the scRNA-seq workflow"),

                 h4("Workflow Guide"),
                 tags$ul(
                   tags$li(strong("Pre-processing:"), "Cleans raw data by normalizing counts and identifying variable features."),
                   tags$li(strong("Dimension Reduction:"), "Compresses data using PCA and visualizes local structure with UMAP."),
                   tags$li(strong("Clustering:"), "Uses `ClusterUnderRes` to visualize cluster stability across resolutions (Clustree)."),
                   tags$li(strong("Markers:"), "Uses `FindClusterDeg` to rank genes by a custom score: p_FC = (1 - p_val) * avg_log2FC."),
                   tags$li(strong("Annotation:"), "Maps cells to references using `SingleR` (cell-level) or `Clustifyr` (cluster-level)."),
                   tags$li(strong("GO Analysis:"), "Performs Gene Ontology enrichment using `ClusterGo`.")
                 ),

                 br(),
                 h4("Package References"),
                 tags$ul(
                   tags$li("Jiaqi, M. (2025) scClustAnnot: An R package integrate into Seurat V5, meant to aid clustering and annotation step in scRNA-seq analysis. Unpublished. ", a(href="https://github.com/DDMMMAA/scClustAnnot", "GitHub", target="_blank"))
                 ),

                 h4("Other References"),
                 tags$ul(
                   tags$li("3k PBMCs from a Healthy Donor. (n.d.). 10x Genomics. ", a(href="https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0", "Link", target="_blank")),
                   tags$li("Analysis, visualization, and integration of Visium HD spatial datasets with Seurat. (n.d.). Satija Lab. ", a(href="https://satijalab.org/seurat/articles/pbmc3k_tutorial", "Link", target="_blank")),
                   tags$li("Aran, D., Looney, A. P., Liu, L., Wu, E., Fong, V., Hsu, A., Chak, S., Naikawadi, R. P., Wolters, P. J., Abate, A. R., Butte, A. J., & Bhattacharya, M. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology, 20(2), 163–172. ", a(href="https://doi.org/10.1038/s41590-018-0276-y", "DOI", target="_blank")),
                   tags$li("Fu, R., Gillen, A. E., Sheridan, R. M., Tian, C., Daya, M., Hao, Y., Hesselberth, J. R., & Riemondy, K. A. (2020). clustifyr: An R package for automated single-cell RNA sequencing cluster classification. F1000Research. ", a(href="https://doi.org/10.12688/f1000research.22969.2", "DOI", target="_blank")),
                   tags$li("Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology, 42(2), 293–304. ", a(href="https://doi.org/10.1038/s41587-023-01767-y", "DOI", target="_blank")),
                   tags$li("Ouyang, L., Wu, J., Jiang, X., Almeida, D., Wainwright, C. L., Mishkin, P., Zhang, C., Agarwal, S., Slama, K., Ray, A., Schulman, J., Hilton, J., Kelton, F., Miller, L., Simens, M., Askell, A., Welinder, P., Christiano, P., Leike, J., & Lowe, R. (2022). Training language models to follow instructions with human feedback (No. arXiv:2203.02155). arXiv. ", a(href="https://doi.org/10.48550/arXiv.2203.02155", "DOI", target="_blank")),
                   tags$li("Satijalab/seurat-data. (2025). [R]. satijalab. (Original work published 2019) ", a(href="https://github.com/satijalab/seurat-data", "GitHub", target="_blank")),
                   tags$li("Schmiedel, B. J., Singh, D., Madrigal, A., Valdovino-Gonzalez, A. G., White, B. M., Zapardiel-Gonzalo, J., Ha, B., Altay, G., Greenbaum, J. A., McVicker, G., Seumois, G., Rao, A., Kronenberg, M., Peters, B., & Vijayanand, P. (2018). Impact of Genetic Polymorphisms on Human Immune Cell Gene Expression. Cell, 175(6), 1701-1715.e16. ", a(href="https://doi.org/10.1016/j.cell.2018.10.022", "DOI", target="_blank")),
                   tags$li("Wasserstein, R. L., & Lazar, N. A. (2016). The ASA Statement on p-Values: Context, Process, and Purpose. The American Statistician, 70(2), 129–133. ", a(href="https://doi.org/10.1080/00031305.2016.1154108", "DOI", target="_blank")),
                   tags$li("Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T. L., Miller, E., Bache, S. M., Müller, K., Ooms, J., Robinson, D., Seidel, D. P., Spinu, V., … Yutani, H. (2019). Welcome to the Tidyverse. Journal of Open Source Software, 4(43), 1686. ", a(href="https://doi.org/10.21105/joss.01686", "DOI", target="_blank")),
                   tags$li("Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation, 2(3). ", a(href="https://doi.org/10.1016/j.xinn.2021.100141", "DOI", target="_blank")),
                   tags$li("Zappia, L., & Oshlack, A. (2018). Clustering trees: A visualization for evaluating clusterings at multiple resolutions. GigaScience, 7(7), giy083. ", a(href="https://doi.org/10.1093/gigascience/giy083", "DOI", target="_blank"))
                 ),
                 br(),
                 p("Package developed by Jiaqi Ma.")
        ),

        # --- TAB: UMAP ---
        tabPanel("1. UMAP",
                 plotOutput("umap_plot", height = "600px"),
                 helpText("Dimensionality reduction plot showing current clusters/identities.")
        ),

        # --- TAB: Clustree ---
        tabPanel("2. Clustering Tree",
                 plotOutput("clustree_plot", height = "600px"),
                 helpText("Arrows indicate how cells move between clusters as resolution increases. Stable clusters have few crossing arrows.")
        ),

        # --- TAB: Markers ---
        tabPanel("3. Marker Table",
                 h4(textOutput("table_title")),
                 tableOutput("marker_table"),
                 helpText("Markers are sorted by p_FC (significance + fold change).")
        ),

        # --- TAB: Annotation ---
        tabPanel("4. Annotation",
                 plotOutput("anno_plot", height = "600px"),
                 helpText("Predicted cell types projected onto the UMAP.")
        ),

        # --- TAB: GO ---
        tabPanel("5. GO Enrichment",
                 br(),
                 uiOutput("go_ui"),
                 helpText("Dotplot showing enriched Gene Ontology terms for the selected cluster.")
        )
      )
    )
  )
)

# ==============================================================================
# SERVER LOGIC
# ==============================================================================
server <- function(input, output, session) {

  # Reactive storage for the Seurat object and analysis results
  values <- reactiveValues(obj = NULL, deg_res = NULL, go_res = NULL)

  # --- DYNAMIC: Populate References (Clustifyr Only) ---
  # Scans the 'clustifyrdatahub' package for available reference functions
  observe({
    if (requireNamespace("clustifyrdatahub", quietly = TRUE)) {
      all_exports <- getNamespaceExports("clustifyrdatahub")
      c_options <- sort(grep("^ref_", all_exports, value = TRUE))
      updateSelectInput(session, "clustifyr_ref",
                        choices = c_options,
                        selected = "ref_hema_microarray")
    }
  })

  # --- 1. Load Data ---
  observeEvent(input$load_demo, {
    dataset_name <- input$demo_select
    showNotification(paste("Loading", dataset_name, "..."), type = "message")

    # Load data into a clean environment to avoid conflicts
    e <- new.env()
    tryCatch({
      data(list = dataset_name, package = "scClustAnnot", envir = e)
      if (!exists(dataset_name, envir = e)) stop("Dataset not found.")

      values$obj <- e[[dataset_name]]
      output$data_status <- renderText(paste("Loaded:", dataset_name, "\nCells:", ncol(values$obj)))

      # Warn user if data is raw (no PCA calculated yet)
      if (!"pca" %in% names(values$obj@reductions)) {
        showNotification("Raw data loaded. Please run Pre-processing.", type = "warning")
      }
    }, error = function(err) showNotification(paste("Error:", err$message), type = "error"))
  })

  # UPDATED: Simplified .rda loading
  observeEvent(input$upload_rda, {
    req(input$upload_rda)
    tryCatch({
      # Load into environment
      e <- new.env()
      loaded_objects <- load(input$upload_rda$datapath, envir = e)

      if (length(loaded_objects) == 0) {
        stop("The .rda file contains no objects.")
      }

      # Assumption: The file contains a single Seurat object (take the first one)
      obj_name <- loaded_objects[1]
      obj <- e[[obj_name]]

      # Check class
      if (!inherits(obj, "Seurat")) {
        stop(paste("The loaded object '", obj_name, "' is not a Seurat object."))
      }

      values$obj <- obj
      output$data_status <- renderText(paste("Loaded:", ncol(values$obj), "Cells"))

    }, error = function(e) showNotification(paste("Upload Error:", e$message), type = "error"))
  })

  # --- 2. Pre-processing ---
  # Runs standard Seurat workflow: Norm -> VarFeatures -> Scale -> PCA -> Neighbors -> UMAP
  observeEvent(input$run_preprocess, {
    req(values$obj)
    showNotification("Running Standard Workflow...", type = "message")
    tryCatch({
      obj <- values$obj
      obj <- Seurat::NormalizeData(obj)
      obj <- Seurat::FindVariableFeatures(obj)
      obj <- Seurat::ScaleData(obj)
      obj <- Seurat::RunPCA(obj, verbose = FALSE)
      obj <- Seurat::FindNeighbors(obj, dims = 1:10)
      obj <- Seurat::RunUMAP(obj, dims = 1:10)
      values$obj <- obj
      showNotification("Pre-processing Complete!", type = "message")
      updateTabsetPanel(session, "main_tabs", selected = "1. UMAP")
    }, error = function(e) showNotification(paste("Error:", e$message), type = "error"))
  })

  # --- 3. Clustering (ClusterUnderRes) ---
  # Iteratively clusters data across the requested resolution range
  observeEvent(input$run_cluster, {
    req(values$obj)
    showNotification("Running Clustering...", type = "message")
    values$obj <- scClustAnnot::ClusterUnderRes(
      values$obj, start = input$res_range[1], end = input$res_range[2],
      margin = input$res_margin, showPlot = FALSE
    )
    showNotification("Clustering Complete! Check the 'Clustering Tree' tab.", type = "message")
  })

  # Dynamic UI: Shows resolution options ONLY after clustering is run
  output$cluster_res_ui <- renderUI({
    req(values$obj)
    cols <- colnames(values$obj@meta.data)
    res_cols <- grep("^RNA_snn_res", cols, value = TRUE)
    if (length(res_cols) == 0) return(helpText("Run Clustering to see resolution options."))

    tagList(
      selectInput("selected_res", "Display Resolution:", choices = res_cols, selected = res_cols[length(res_cols)]),
      actionButton("set_ident", "Use for Downstream Analysis", class = "btn-success", style = "width: 100%; margin-top: 5px;"),
      helpText("Click to lock this resolution for Marker/Annotation steps.", style = "font-size: 0.8em; color: #666;")
    )
  })

  # Logic to set the active identity based on user selection
  observeEvent(input$set_ident, {
    req(values$obj, input$selected_res)
    target_res <- input$selected_res
    Idents(values$obj) <- values$obj@meta.data[[target_res]]
    showNotification(paste("Active Identity updated to:", target_res), type = "message")
  })

  output$clustree_plot <- renderPlot({
    req(values$obj)
    clustree::clustree(values$obj, prefix = "RNA_snn_res.")
  })

  output$umap_plot <- renderPlot({
    req(values$obj)
    target_col <- input$selected_res
    # Plot using selected resolution column if available, otherwise default
    if (!is.null(target_col) && target_col %in% colnames(values$obj@meta.data)) {
      DimPlot(values$obj, reduction = "umap", group.by = target_col, label = TRUE) +
        ggtitle(paste("Visualizing:", target_col))
    } else {
      DimPlot(values$obj, reduction = "umap", label = TRUE) +
        ggtitle("Active Identity")
    }
  })

  # --- 4. Markers (FindClusterDeg) ---
  observeEvent(input$run_deg, {
    req(values$obj)
    showNotification("Finding Markers based on Active Identity...", type = "message")
    # Calculates p_FC score for every marker
    values$deg_res <- scClustAnnot::FindClusterDeg(
      values$obj, logfc_min = input$deg_fc, pval_max = input$deg_pval
    )
    showNotification("Markers Found!", type = "message")
  })

  # Dynamic UI to pick which cluster to view in the table
  output$deg_cluster_ui <- renderUI({
    req(values$deg_res)
    clusters <- names(values$deg_res$by_cluster)
    if (length(clusters) == 0) return(NULL)
    selectInput("deg_view_cluster", "View Cluster Markers:", choices = c("All Clusters", clusters), selected = "All Clusters")
  })

  output$table_title <- renderText({
    req(values$deg_res)
    if (is.null(input$deg_view_cluster) || input$deg_view_cluster == "All Clusters") {
      return("Top 50 Markers (All Clusters Combined)")
    } else {
      return(paste("Markers for Cluster", input$deg_view_cluster))
    }
  })

  output$marker_table <- renderTable({
    req(values$deg_res)
    view_choice <- input$deg_view_cluster
    if (is.null(view_choice) || view_choice == "All Clusters") {
      head(values$deg_res$combined, 50)
    } else {
      clus_df <- values$deg_res$by_cluster[[view_choice]]
      if (is.null(clus_df) || nrow(clus_df) == 0) {
        data.frame(Message = "No markers found for this cluster.")
      } else {
        head(clus_df, 50)
      }
    }
  })

  # --- 5. Annotation ---
  observeEvent(input$run_anno, {
    req(values$obj)
    showNotification("Annotating...", type = "message")

    if (input$anno_method == "SingleR") {
      # Cell-level annotation via SingleR
      if (requireNamespace("celldex", quietly = TRUE)) {
        tryCatch({
          values$obj <- scClustAnnot::SinglerAnnote(
            values$obj,
            name = input$singler_ref,
            version = input$singler_ver
          )
          showNotification("SingleR Annotation Complete!", type = "message")
        }, error = function(e) {
          showNotification(paste("SingleR Error:", e$message), type = "error")
        })
      } else { showNotification("Package 'celldex' missing.", type = "error") }
    } else {
      # Cluster-level annotation via Clustifyr
      if (requireNamespace("clustifyrdatahub", quietly = TRUE)) {
        ref_name <- input$clustifyr_ref
        tryCatch({
          # Dynamically call the selected reference function
          ref_func <- get(ref_name, envir = asNamespace("clustifyrdatahub"))
          ref <- ref_func()
          values$obj <- scClustAnnot::ClustifyrAnnote(values$obj, ref_mat = ref)
          showNotification("Clustifyr Annotation Complete!", type = "message")
        }, error = function(e) {
          showNotification(paste("Clustifyr Error:", e$message), type = "error")
        })
      } else { showNotification("Package 'clustifyrdatahub' missing.", type = "error") }
    }
  })

  output$anno_plot <- renderPlot({
    req(values$obj)
    group_col <- if(input$anno_method == "SingleR") "SingleR.labels" else "Clustifyr.labels"
    if (group_col %in% colnames(values$obj@meta.data)) {
      DimPlot(values$obj, group.by = group_col, label = TRUE) + ggtitle(input$anno_method)
    }
  })

  # --- 6. GO Analysis (ClusterGo) ---
  observeEvent(input$run_go, {
    req(values$deg_res)
    selected_org <- input$go_org
    showNotification(paste("Running GO using", selected_org, "..."), type = "message")

    # Check if the organism database is installed
    if (requireNamespace(selected_org, quietly = TRUE)) {
      tryCatch({
        values$go_res <- scClustAnnot::ClusterGo(
          values$deg_res,
          org_db = selected_org,
          ont = input$go_ont
        )
        showNotification("GO Analysis Complete!", type = "message")
      }, error = function(e) {
        showNotification(paste("GO Error:", e$message), type = "error")
      })
    } else { showNotification(paste("Package", selected_org, "missing."), type = "error") }
  })

  output$go_ui <- renderUI({
    req(values$go_res)
    valid_clusters <- names(values$go_res)[!is.na(values$go_res)]
    if (length(valid_clusters) == 0) return(h4("No significant GO terms."))

    tagList(
      selectInput("go_cluster_select", "Select Cluster:", choices = valid_clusters),
      plotOutput("go_dotplot", height = "600px")
    )
  })

  output$go_dotplot <- renderPlot({
    req(values$go_res, input$go_cluster_select)
    res <- values$go_res[[input$go_cluster_select]]
    if (!is.na(res)) {
      enrichplot::dotplot(res, showCategory = 15) + ggtitle(input$go_cluster_select)
    }
  })
}

# Run the App
shinyApp(ui, server)
