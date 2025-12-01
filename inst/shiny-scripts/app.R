library(shiny)
library(Seurat)
library(ggplot2)
library(enrichplot)
library(scClustAnnot)

# Define UI
ui <- fluidPage(

  titlePanel("scClustAnnot Dashboard"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # --- 1. LOAD DATA ---
      h4("1. Load Data"),
      selectInput("demo_select", "Choose Built-in Data:",
                  choices = c("sub_pbmc", "pbmc"),
                  selected = "sub_pbmc"),
      actionButton("load_demo", "Load Selected Data", icon = icon("play")),
      fileInput("upload_rds", "Or Upload .rds", accept = ".rds"),
      verbatimTextOutput("data_status"),
      hr(),

      # --- 2. PRE-PROCESSING ---
      h4("2. Pre-processing"),
      helpText("Normalize, Find Features, Scale, PCA, UMAP"),
      actionButton("run_preprocess", "Run Standard Workflow", class = "btn-primary"),
      hr(),

      # --- 3. CLUSTERING ---
      h4("3. Clustering"),
      sliderInput("res_range", "Resolution Range:", min = 0.1, max = 1.5, value = c(0.1, 0.8), step = 0.1),
      numericInput("res_margin", "Step Size:", value = 0.1, step = 0.1),
      actionButton("run_cluster", "Run Clustering"),
      br(), br(),
      uiOutput("cluster_res_ui"),
      hr(),

      # --- 4. MARKERS ---
      h4("4. Find Markers"),
      numericInput("deg_fc", "Min LogFC:", value = 0.25, step = 0.05),
      numericInput("deg_pval", "Max P-Value:", value = 0.05, step = 0.01),
      actionButton("run_deg", "Find DEGs"),
      br(), br(),
      uiOutput("deg_cluster_ui"),
      hr(),

      # --- 5. ANNOTATION ---
      h4("5. Annotation"),
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

      # --- 6. GO ENRICHMENT (UPDATED) ---
      h4("6. GO Analysis"),
      selectInput("go_org", "Organism:",
                  choices = c("Human (org.Hs.eg.db)" = "org.Hs.eg.db",
                              "Mouse (org.Mm.eg.db)" = "org.Mm.eg.db")),
      selectInput("go_ont", "Ontology:", choices = c("BP", "MF", "CC")),
      actionButton("run_go", "Run Enrichment")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("1. UMAP", plotOutput("umap_plot", height = "600px")),
        tabPanel("2. Clustering Tree", plotOutput("clustree_plot", height = "600px")),

        tabPanel("3. Marker Table",
                 h4(textOutput("table_title")),
                 tableOutput("marker_table")),

        tabPanel("4. Annotation", plotOutput("anno_plot", height = "600px")),
        tabPanel("5. GO Enrichment", br(), uiOutput("go_ui"))
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {

  values <- reactiveValues(obj = NULL, deg_res = NULL, go_res = NULL)

  # --- DYNAMIC: Populate References ---
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
    e <- new.env()
    tryCatch({
      data(list = dataset_name, package = "scClustAnnot", envir = e)
      if (!exists(dataset_name, envir = e)) stop("Dataset not found.")
      values$obj <- e[[dataset_name]]
      output$data_status <- renderText(paste("Loaded:", dataset_name, "\nCells:", ncol(values$obj)))
      if (!"pca" %in% names(values$obj@reductions)) {
        showNotification("Raw data loaded. Please run Pre-processing.", type = "warning")
      }
    }, error = function(err) showNotification(paste("Error:", err$message), type = "error"))
  })

  observeEvent(input$upload_rds, {
    req(input$upload_rds)
    tryCatch({
      obj <- readRDS(input$upload_rds$datapath)
      if (!inherits(obj, "Seurat")) stop("File is not a Seurat object")
      values$obj <- obj
      output$data_status <- renderText(paste("Loaded:", ncol(values$obj), "Cells"))
    }, error = function(e) showNotification(paste("Upload Error:", e$message), type = "error"))
  })

  # --- 2. Pre-processing ---
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

  # --- 3. Clustering ---
  observeEvent(input$run_cluster, {
    req(values$obj)
    showNotification("Running Clustering...", type = "message")
    values$obj <- scClustAnnot::ClusterUnderRes(
      values$obj, start = input$res_range[1], end = input$res_range[2],
      margin = input$res_margin, showPlot = FALSE
    )
    showNotification("Clustering Complete!", type = "message")
  })

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
    if (!is.null(target_col) && target_col %in% colnames(values$obj@meta.data)) {
      DimPlot(values$obj, reduction = "umap", group.by = target_col, label = TRUE) +
        ggtitle(paste("Visualizing:", target_col))
    } else {
      DimPlot(values$obj, reduction = "umap", label = TRUE) +
        ggtitle("Active Identity")
    }
  })

  # --- 4. Markers ---
  observeEvent(input$run_deg, {
    req(values$obj)
    showNotification("Finding Markers based on Active Identity...", type = "message")
    values$deg_res <- scClustAnnot::FindClusterDeg(
      values$obj, logfc_min = input$deg_fc, pval_max = input$deg_pval
    )
    showNotification("Markers Found!", type = "message")
  })

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
      if (requireNamespace("clustifyrdatahub", quietly = TRUE)) {
        ref_name <- input$clustifyr_ref
        tryCatch({
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

  # --- 6. GO Analysis (UPDATED) ---
  observeEvent(input$run_go, {
    req(values$deg_res)
    selected_org <- input$go_org # Get user selection (e.g. "org.Hs.eg.db")

    showNotification(paste("Running GO using", selected_org, "..."), type = "message")

    # Dynamically check for the selected organism package
    if (requireNamespace(selected_org, quietly = TRUE)) {
      tryCatch({
        values$go_res <- scClustAnnot::ClusterGo(
          values$deg_res,
          org_db = selected_org, # Pass selected DB
          ont = input$go_ont
        )
        showNotification("GO Analysis Complete!", type = "message")
      }, error = function(e) {
        showNotification(paste("GO Error:", e$message), type = "error")
      })
    } else {
      showNotification(paste("Package", selected_org, "is missing. Please install it."), type = "error")
    }
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

shinyApp(ui, server)
