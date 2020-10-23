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
ui <- fluidPage(
    theme = "bootstrap.css",
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
                         tags$h5("TODO: clean up some weird A/B interaction, B trumps A currently."),
                         disabled(selectInput("sel_A_clust", label = "A group clusters", choices = "", multiple = TRUE)),
                         disabled(selectInput("sel_B_clust", label = "B group clusters", choices = "", multiple = TRUE)),
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
                         ),
                         withSpinner(plotOutput("plot_tsne_clusters", width = "600px", height = "600px"))
                     ),
                     mainPanel(
                         withSpinner(plotOutput("plot_A_group", width = "600px", height = "600px",
                                                brush = brushOpts(id = "plot_A_brush"))),
                         tags$h3("Refine selection"),
                         fluidRow(
                             column(width = 4,
                                    actionButton("btn_add_A", "Add A"),
                                    actionButton("btn_rm_A", "Remove A"),
                                    actionButton("btn_limit_A", "Limit A"),
                             ), 
                             column(width = 4,
                                    actionButton("btn_add_B", "Add B"),
                                    actionButton("btn_rm_B", "Remove B"),
                                    actionButton("btn_limit_B", "Limit B")
                             )
                         )
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
        tsne_dt.clust$group = "bg"
        tsne_clust(tsne_dt.clust)
        # updateSelectInput(session, "sel_A_clust", choices = cl, selected = cl[1])
        # updateSelectInput(session, "sel_B_clust", choices = cl, selected = cl[2])
    })
    
    observeEvent({
        tsne_clust()
    },
    {
        if(is.null(tsne_clust())){
            disable("sel_A_clust")
            disable("sel_B_clust")
        }
        cl = sort(unique(tsne_clust()$cluster_id))
        cl = cl[order(as.numeric(sub("[a-zA-Z ]+", "", cl)))]
        enable("sel_A_clust")
        enable("sel_B_clust")
        updateSelectInput(session, "sel_A_clust", choices = cl, selected = cl[1])
        updateSelectInput(session, "sel_B_clust", choices = cl, selected = cl[2])
        prev_A_clust = cl[1]
        prev_B_clust = cl[2]
        sel_A_ids(tsne_clust()[cluster_id %in% cl[1]]$bcr_patient_barcode)
        sel_B_ids(tsne_clust()[cluster_id %in% cl[2]]$bcr_patient_barcode)
        
    })
    
    
    sel_A_ids = reactiveVal(character())
    sel_B_ids = reactiveVal(character())
    prev_A_clust = character()
    prev_B_clust = character()
    
    update_ids = function(current_clust, previous_clust, current_ids){
        added_clust = setdiff(current_clust, previous_clust)
        removed_clust = setdiff(previous_clust, current_clust)
        if(length(added_clust) > 0 & length(removed_clust) > 0){
            stop("inconceivable! new and missing longer than 0.")
        }
        tsne_dt = tsne_clust()
        new_ids = current_ids
        if(length(added_clust) > 0){
            new_ids = union(
                current_ids,
                tsne_dt[cluster_id %in% added_clust]$bcr_patient_barcode          
            )
        }
        if(length(removed_clust) > 0){
            new_ids = setdiff(
                current_ids,
                tsne_dt[cluster_id %in% removed_clust]$bcr_patient_barcode          
            )
        }
        new_ids
    }
    
    observeEvent({
        input$sel_A_clust
    }, {
        new_ids = update_ids(input$sel_A_clust, 
                             prev_A_clust,
                             isolate(sel_A_ids()))
        prev_A_clust = input$sel_A_clust
        sel_A_ids(new_ids)
    })
    
    observeEvent({
        input$sel_B_clust
    }, {
        new_ids = update_ids(input$sel_B_clust, 
                             prev_B_clust,
                             isolate(sel_B_ids()))
        prev_B_clust = input$sel_B_clust
        sel_B_ids(new_ids)
    })
    
    output$plot_A_group = renderPlot({
        tsne_dt = tsne_clust()
        # tsne_dt[, sel := cluster_id %in% input$sel_A_clust]
        # ggplot(tsne_dt, aes(x = x, y = y)) +
        #     geom_point(data= tsne_dt[sel == FALSE], color = 'gray', size = .3) +
        #     geom_point(data= tsne_dt[sel == TRUE], color = 'blue', size = .8) 
        tsne_dt$group = "bg"
        # tsne_dt[cluster_id %in% input$sel_A_clust, group := "A"]
        # tsne_dt[cluster_id %in% input$sel_B_clust, group := "B"]
        
        tsne_dt[bcr_patient_barcode %in% sel_A_ids(), group := "A"]
        tsne_dt[bcr_patient_barcode %in% sel_B_ids(), group := "B"]
        
        ggplot(tsne_dt, aes(x = x, y = y)) +
            geom_point(data= tsne_dt[group == "bg"], color = 'gray', size = .3) +
            geom_point(data= tsne_dt[group != "bg"], aes(color = group), size = .8) +
            scale_color_manual(values = c("A" = "red", "B" = "blue")) +
            theme(panel.background = element_blank(),
                  panel.grid = element_blank())
        
    })
    
    # actionButton("btn_add_A", "Add A"),
    # actionButton("btn_rm_A", "Remove A"),
    # actionButton("btn_limit_A", "Limit A"),
    # actionButton("btn_add_B", "Add B"),
    # actionButton("btn_rm_B", "Remove B"),
    # actionButton("btn_limit_B", "Limit B")
    #plot A buttons
    observeEvent({
        input$btn_add_A
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_A_ids(
            union(isolate(sel_A_ids()), 
                  ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_rm_A
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_A_ids(
            setdiff(isolate(sel_A_ids()), 
                  ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_limit_A
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_A_ids(
            intersect(isolate(sel_A_ids()), 
                  ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_add_B
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_B_ids(
            union(isolate(sel_B_ids()), 
                  ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_rm_B
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_B_ids(
            setdiff(isolate(sel_B_ids()), 
                    ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_limit_B
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        sel_B_ids(
            intersect(isolate(sel_B_ids()), 
                      ids_in_rng)    
        )
    })
    
    output$plot_B_group = renderPlot({
        tsne_dt = tsne_clust()
        tsne_dt[, sel := cluster_id %in% input$sel_B_clust]
        ggplot(tsne_dt, aes(x = x, y = y)) +
            geom_point(data= tsne_dt[sel == FALSE], color = 'gray', size = .3) +
            geom_point(data= tsne_dt[sel == TRUE], color = 'blue', size = .8) 
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
