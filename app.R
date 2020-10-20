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
    
    observeEvent({
        input$sel_data
    }, {
        sel = input$sel_data
        if(!sel %in% names(expression_files)){
            stop(sel, " not found in expression_files.")
        }
        if(is.null(expression_loaded[[sel]])){
            expression_loaded[[sel]] = load_expression(expression_files[[sel]])
        }
        tcga_data(expression_loaded[[sel]])
    })
    observe({
        showNotification(paste0("expression: ", nrow(tcga_data()), " rows x ", ncol(tcga_data()), " columns loaded."))
    })
    ### Sample Type
    observeEvent({
        input$sel_sample_type_filter
        tcga_data()
    }, {
        req(tcga_data())
        sample_codes = sub("[A-Z]", "", sapply(strsplit(colnames(tcga_data()), "-"), function(x)x[4]))
        sample_types = code2type[sample_codes]
        k = toupper(sample_types) %in% toupper(input$sel_sample_type_filter)
        # tsne_dt[, sample_code := tstrsplit(bcr_patient_barcode, "-", keep = 4)]
        # tsne_dt[, sample_code := sub("[A-Z]", "", sample_code)]
        # tsne_dt[, sample_type := code2type[sample_code]]
        tsne_input(tcga_data()[, k])
        tsne_res(NULL)
    })
    
    
    ### patient metadata
    observeEvent({
        input$sel_data
    }, {
        sel = input$sel_data
        if(!sel %in% names(clinical_loaded)){
            stop(sel, " not found in loaded metadata.")
        }
        # if(is.null(clinical_loaded[[sel]])){
        #     clinical_loaded[[sel]] = load_metadata(clinical_files[[sel]])
        # }
        meta_data(clinical_loaded[[sel]])
    })
    observe({
        showNotification(paste0("metadata: ", nrow(meta_data()), " rows x ", ncol(meta_data()), " columns loaded."))
    })
    
    
    
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
    
    
    tsne_res = server_tsne(input, output, session, tsne_input, valid_genes, meta_data, code2type, FACET_VAR)
    server_gene_xy(input, output, session, tsne_res, tcga_data, vis_gene)
    
    server_upload(input, output, session, gene_table)
    
}

# Run the application 
shinyApp(ui = ui, server = server)
