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
    
    # Application title
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
                         radioButtons("sel_color_by", label = "Color By", choices = c("sample type", "PAM50")),
                         # radioButtons("sel_facet_by", label = "Facet By", choices = c("sample type", "PAM50")),
                         
                         selectizeInput("txtGene", label = "Select Gene To Plot", choices = NULL),
                         
                     ),
                     
                     # Show a plot of the generated distribution
                     mainPanel(
                         withSpinner(plotOutput("plot_tsne", width = "600px", height = "600px")),
                         withSpinner(plotOutput("plot_tsne_gene", width = "600px", height = "600px"))
                     )
                 )
        ),
        tabPanel("Add Gene Set",
                 
                 sidebarLayout(
                     sidebarPanel(
                         tabsetPanel(
                             id = "gene_list_method",
                             tabPanel("Paste", value = "paste",
                                      textAreaInput("txt_genes", label = "Custom Genes", value = "paste genes here")#,
                             ),
                             ui_upload()
                         )
                     ),
                     mainPanel(
                         DT::dataTableOutput("DT_PasteGenes_DataFrame")
                     )
                 )
        )
        
    )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
    ### TCGA expression data
    tcga_data = reactiveVal()
    tsne_input = reactiveVal()
    ### Genes to use in t-sne
    input_genes = reactiveVal()
    #genes parsed for paste/upload
    parsed_genes = reactiveVal()
    #subset of parsed genes in tcga expression data
    valid_genes = reactiveVal()
    #table shown in Add Gene Set main panel, includes parsed genes and other data used for selection
    gene_table = reactiveVal()
    #metadata for patient entries
    meta_data = reactiveVal()
    tsne_res = reactiveVal()
    
    
    ## watch gene inputs
    observeEvent({
        input$sel_gene_list
        # input$txt_genes
    }, {
        sel = input$sel_gene_list
        gl = NULL
        gl = gene_lists[[sel]]
        # message(paste(gl, collapse = ", "))
        gl = sort(unique(gl))
        input_genes(gl)
    })
    
    observeEvent({
        input$txt_genes
    }, {
        gl = parse_gl(input$txt_genes)
        parsed_genes(gl)
    })
    
    observe({
        if(input$gene_list_method == "paste"){
            showNotification("DEBUG list by paste")
            gl = parsed_genes()
            if(is.null(gl)){
                gene_table(data.frame())  
            }else if(length(gl) == 0){
                gene_table(data.frame())  
            }else{
                gene_table(data.frame(gene_name = gl))
            }    
        }
        
    })
    
    output$DT_PasteGenes_DataFrame = DT::renderDataTable({
        DT::datatable(gene_table())    
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
    server_upload(input, output, session, gene_table)
    
}

# Run the application 
shinyApp(ui = ui, server = server)
