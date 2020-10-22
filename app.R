#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(data.table)
library(ggplot2)
library(BiocFileCache)
library(digest)
library(shinyjs)
library(shinycssloaders)

source("setup_gene_lists.R")
source("setup_clinical.R")
source("setup_tcga_expression.R")
source("app_module_expression_matrix.R")
source("app_module_upload.R")
source("app_module_tsne.R")
source("app_module_gene_xy_vis.R")
source("functions.R")
# based on analyze_BRCA_tsne_allGenes_plusCells.R

code2type = c("01" = "tumor", "06" = "metastasis", "11" = "normal")
FACET_VAR = list(NONE = "none", SAMPLE_TYPE = "sample type", PAM50 = "PAM50")
clean_list = function(l){
    tmp = unlist(l)
    names(tmp) = NULL
    tmp
}





# Define UI for application that draws a histogram
ui <- fluidPage(theme = "bootstrap.css",
                # Application title
                useShinyjs(),
                titlePanel("TCGA t-sne"),
                tabsetPanel(
                    tabPanel("Main",
                             # Sidebar with a slider input for number of bins 
                             sidebarLayout(
                                 sidebarPanel(
                                     selectInput("sel_data", label = "TCGA data", choices = c("BRCA")),
                                     checkboxGroupInput("sel_sample_type_filter", label = "Samples Included", 
                                                        choices = c("normal", "tumor", "metastasis"), 
                                                        selected = c("normal", "tumor", "metastasis")),
                                     radioButtons("sel_facet_var", label = "Facet By", choices = clean_list(FACET_VAR)), #c("none", "sample type", "PAM50")),
                                     radioButtons("sel_gene_list", label = "Gene List", choices = c(names(gene_lists), "custom"), inline = TRUE, selected = "PAM50"),
                                     disabled((selectInput("sel_custom_gene_set", label = "Custom gene lists", choices = ""))),
                                     radioButtons("sel_color_by", label = "Color By", choices = c("sample type", "PAM50")),
                                     # radioButtons("sel_facet_by", label = "Facet By", choices = c("sample type", "PAM50")),
                                     
                                     selectInput("txtGene", label = "Select Gene To Plot", choices = NULL),
                                 ),
                                 mainPanel(
                                     withSpinner(plotOutput("plot_tsne", width = "600px", height = "600px")),
                                     withSpinner(plotOutput("plot_tsne_gene", width = "600px", height = "600px"))
                                 )
                             )
                    ),
                    tabPanel("Add Gene Set",
                             ui_upload()
                    ),
                    tabPanel("DE",
                             sidebarLayout(
                                 sidebarPanel(
                                     disabled(selectInput("sel_A_clust", label = "A group clusters", choices = "")),
                                     disabled(selectInput("sel_B_clust", label = "B group clusters", choices = "")),
                                     tabsetPanel(id = "tabs_cluster_method",
                                                 tabPanel("knn", #id = "knn", 
                                                          numericInput("num_nn", label = "Nearest neighbors", value = 5, min = 2, max = Inf)
                                                 ), 
                                                 tabPanel("kmeans", #id = "kmeans", 
                                                          numericInput("num_kmeans", label = "k", value = 5, min = 2, max = Inf)
                                                 ), 
                                                 tabPanel("hierarchical", #id = "hierarchical", 
                                                          numericInput("num_clust", label = "n_clust", value = 5, min = 2, max = Inf)
                                                 )
                                     )
                                 ),
                                 mainPanel(
                                     withSpinner(plotOutput("plot_tsne_clusters", width = "600px", height = "600px"))
                                 )
                             ))
                )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
    ### TCGA expression data
    tcga_data = reactiveVal()
    tsne_input = reactiveVal()
    ### Genes to use in t-sne
    input_genes = reactiveVal()
    #subset of parsed genes in tcga expression data
    valid_genes = reactiveVal()
    #table shown in Add Gene Set main panel, includes parsed genes and other data used for selection
    gene_table = reactiveVal()
    #all added custom gene sets
    custom_gene_sets = reactiveVal(list())
    
    #metadata for patient entries
    meta_data = reactiveVal()
    #result of tsne
    tsne_res = reactiveVal()
    #tsne_res with clustering applied
    tsne_clust = reactiveVal()
    ## watch gene inputs
    observeEvent({
        input$sel_gene_list
        # input$txt_genes
    }, {
        sel = input$sel_gene_list
        if(sel == "custom"){
            gl = custom_gene_sets()[[input$sel_custom_gene_set]]
        }else{
            gl = NULL
            gl = gene_lists[[sel]]
            # message(paste(gl, collapse = ", "))
            gl = sort(unique(gl))
        }
        input_genes(gl)
    })
    
    
    ##
    vis_gene = reactiveVal()
    observeEvent({
        tcga_data()
    }, {
        req(tcga_data())
        def = vis_gene()
        if(is.null(def)) def = ""
        all_genes = rownames(tcga_data())
        if(def %in% all_genes){
            updateSelectizeInput(session, 'txtGene', choices = all_genes, selected = def, server = TRUE)
        }else{
            updateSelectizeInput(session, 'txtGene', choices = all_genes, server = TRUE)
        }
        
    })
    observeEvent({
        input$txtGene   
        tcga_data()
    }, {
        req(tcga_data())
        if(input$txtGene %in% rownames(tcga_data())){
            vis_gene(input$txtGene)
        }
    })
    
    ## compare input genes with available expression data
    observeEvent({
        tcga_data()
        input_genes()
    }, {
        req(tcga_data())
        gl = input_genes()
        
        missed = setdiff(gl, rownames(tcga_data()))
        if(length(missed) > 0){
            showNotification(paste("genes not present in TCGA:", paste(missed, collapse = ", ")), type = "warning")
        }
        valid_genes(setdiff(gl, missed))
    })
    ## notify about genes loaded
    observeEvent({
        valid_genes()
    }, {
        showNotification(paste0(length(valid_genes()), " genes loaded from ", input$sel_gene_list, "."))
    })
    
    observeEvent({
        custom_gene_sets()
    }, {
        gls = custom_gene_sets()
        if(length(gls) > 0) enable("sel_custom_gene_set")
        curr_sel = input$sel_custom_gene_set
        if(curr_sel %in% names(gls)){
            updateSelectInput(session, "sel_custom_gene_set", choices = names(gls), selected = curr_sel)
        }else{
            updateSelectInput(session, "sel_custom_gene_set", choices = names(gls))
        }
        showNotification(paste(length(gls), " custom gene sets"))
    })
    
    output$plot_tsne_clusters = renderPlot({
        req(tsne_clust())
        tsne_dt = tsne_clust()
        tsne_dt = merge(tsne_dt, meta_data(), by = "submitter_id")
        
        p = ggplot(tsne_dt, aes(x = x, y = y, color = cluster_id)) + 
            geom_point() + 
            coord_fixed() +
            labs(x = "", y = "", title = "t-sne of TCGA samples", subtitle = paste(input$sel_gene_list, "gene list"))
        p
    })
    
    observe({
        if(is.null(tsne_res())){
            tsne_clust(NULL)
        }
        req(tsne_res())
        tsne_dt = tsne_res()
        clust_method = input$tabs_cluster_method
        showNotification(paste("clustering method is", clust_method))
        if(clust_method == "knn"){
            tsne_dt.clust = nn_clust(tsne_dt, nn = input$num_nn)
        }else if(clust_method == "kmeans"){
            tsne_dt.clust = km_clust(tsne_dt, k = input$num_kmeans)
        }else if(clust_method == "hierarchical"){
            tsne_dt.clust = h_clust(tsne_dt, n_clust = input$num_clust)
        }else{
            stop("Unrecognized clustering method, ", clust_method)
        }
        tsne_clust(tsne_dt.clust)
        
    })
    
    server_tsne(input, output, session, tsne_res, tsne_input, valid_genes, meta_data, code2type, FACET_VAR)
    server_expression_matrix(input, output, session,
                             expression_files,
                             expression_loaded,
                             tcga_data,
                             clinical_loaded,
                             meta_data,
                             code2type,
                             tsne_input,
                             tsne_res)
    server_gene_xy(input, output, session, tsne_res, tcga_data, vis_gene)
    server_upload(input, output, session, gene_table, tcga_data, custom_gene_sets)
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
