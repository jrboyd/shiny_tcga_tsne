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
source("setup_gene_lists.R")
source("setup_tcga_expression.R")
# based on analyze_BRCA_tsne_allGenes_plusCells.R

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("TCGA t-sne"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("sel_data", label = "TCGA data", choices = c("BRCA")),
            radioButtons("sel_gene_list", label = "Gene List", choices = c("PAM50", "custom")),
            textInput("txt_genes", label = "Custom Genes", value = "paste genes here")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("plot_tsne")
        )
    )
)

parse_gl = function(txt){
    gl = strsplit(txt, "[, \n]")[[1]]
    gl[gl != ""]
}

load_expression = function(f){
    dat = fread(f)
    mat = as.matrix(dat[,-1])
    rownames(mat) = dat$gene_name
    mat
}

# Define server logic required to draw a histogram
server <- function(input, output) {
    ### TCGA expression data
    tcga_data = reactiveVal()
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
        showNotification(paste0(nrow(tcga_data()), " rows x ", ncol(tcga_data()), " columns loaded."))
    })
    
    ### Genes to use in t-sne
    input_genes = reactiveVal()
    tsne_genes = reactiveVal()
    
    observeEvent({
        input$sel_gene_list
        input$txt_genes
    }, {
        sel = input$sel_gene_list
        gl = NULL
        if(!sel %in% names(gene_lists)){
            gl = parse_gl(input$txt_genes)
            # stop(sel, " not found in gene lists.")
        }else{
            gl = gene_lists[[sel]]
        }
        # message(paste(gl, collapse = ", "))
        gl = sort(unique(gl))
        input_genes(gl)
    })
    
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
        tsne_genes(setdiff(gl, missed))
    })
    
    observeEvent({
        tsne_genes()
    }, {
        showNotification(paste0(length(tsne_genes()), " genes loaded from ", input$sel_gene_list, "."))
    })
    
    ### Running t-sne
    tsne_res = reactiveVal()
    observeEvent({
        tcga_data()
        tsne_genes()
    }, {
        
    })
    
    ### Plot t-sne
    output$plot_tsne <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = 10 + 1)
        
        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
