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
source("app_module_point_selection.R")
source("app_module_diff_expression.R")
source("functions.R")
# based on analyze_BRCA_tsne_allGenes_plusCells.R

code2type = c("01" = "tumor", "06" = "metastasis", "11" = "normal")
FACET_VAR = list(NONE = "none", SAMPLE_TYPE = "sample type", PAM50 = "PAM50")
clean_list = function(l){
    tmp = unlist(l)
    names(tmp) = NULL
    tmp
}

reinit_data = function(sel){
    sel = input$sel_data
    #meta data
    if(!sel %in% names(clinical_loaded)){
        stop(sel, " not found in loaded metadata.")
    }
    # if(is.null(clinical_loaded[[sel]])){
    #     clinical_loaded[[sel]] = load_metadata(clinical_files[[sel]])
    # }
    meta_data(clinical_loaded[[sel]])
    
    #expression data
    if(!sel %in% names(expression_files)){
        stop(sel, " not found in expression_files.")
    }
    if(is.null(expression_loaded[[sel]])){
        expression_loaded[[sel]] = load_expression(expression_files[[sel]])
    }
    dat = expression_loaded[[sel]]
    exp_dt = as.data.table(dat)
    exp_dt$gene_name = rownames(dat)
    exp_dt = melt(exp_dt, id.vars = "gene_name")
    exp_dt = exp_dt[, .(value = max(value)), .(gene_name, variable)]
    # exp_dt = unique(exp_mat$gene_name[duplicated(exp_mat$gene_name)])
    exp_dt = dcast(exp_dt, gene_name~variable, value.var = "value")
    exp_mat = as.matrix(exp_dt[, -1])
    rownames(exp_mat) = exp_dt$gene_name
    
    stopifnot(!any(duplicated(rownames(exp_mat))))
    tcga_data(exp_mat)
    
    #reset downstream
    
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
                         selectInput("sel_data", label = "TCGA data", choices = names(expression_files)),
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
                 ui_point_selection(),
                 tabsetPanel(
                     tabPanel("simple",
                              withSpinner(plotOutput("plot_volcano", width = "400px", height = "400px")),
                              withSpinner(plotOutput("plot_boxes", width = "400px", height = "400px"))
                              
                     ),
                     tabPanel("DESeq2",
                              actionButton("btn_runDE", "Run DE"),
                              withSpinner(DT::dataTableOutput("dt_DE_res"))
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
    
    
    
    
    #running tsne
    server_expression_matrix(input, output, session,
                             expression_files,
                             expression_loaded,
                             tcga_data,
                             clinical_loaded,
                             meta_data,
                             code2type,
                             tsne_input,
                             tsne_res)
    
    server_tsne(input, output, session, tsne_res, tsne_input, valid_genes, meta_data, code2type, FACET_VAR)
    
    #the gene expression mapped to tsne space
    server_gene_xy(input, output, session, tsne_res, tcga_data, vis_gene)
    #upload user data for gene lists
    server_upload(input, output, session, gene_table, tcga_data, custom_gene_sets)
    #interface to select A and B set of points from scatterplot
    sample_groups = server_point_selection(input, output, session, tsne_clust = tsne_clust, meta_data = meta_data, tsne_res = tsne_res)
    
    DE_res = reactiveVal()
    
    observeEvent({
        sample_groups()
    },{
        req(sample_groups())
        showNotification(paste0("A ", length(sample_groups()$A), "\n",
                                "B ", length(sample_groups()$B)))
    })
    
    observeEvent({
        input$btn_runDE
    }, {
        req(tsne_input())
        req(sample_groups())
        diff_res = run_DE(tsne_input(), sample_groups()$A(), sample_groups()$B())
        DE_res(diff_res)
    })
    
    observeEvent({
        DE_res()
    }, {
        req(DE_res())
        showNotification(paste("diff gene count:", nrow(DE_res())))
    })
    
    output$dt_DE_res = DT::renderDataTable({
        if(is.null(DE_res())){
            DT::datatable(data.frame(waiting = "", for_ = "", data = ""))    
        }else{
            DT::datatable(DE_res())    
        }
        
    })
    
    DE_fast_raw = reactiveVal()
    DE_fast_res = reactiveVal()
    
    observeEvent({
        tsne_input()
        sample_groups()
        sample_groups()$A()
        sample_groups()$B()
    }, {
        req(tsne_input())
        req(sample_groups())
        if(length(sample_groups()$A()) == 0 | length(sample_groups()$B()) == 0){
            DE_fast_raw(NULL)
            DE_fast_res(NULL)
        }else{
            dt = run_group.fast(tsne_input(), sample_groups()$A(), sample_groups()$B())
            DE_fast_raw(dt)
            p_dt = run_DE.fast(dt)   
            DE_fast_res(p_dt)    
        }
        
    })
    
    output$plot_boxes = renderPlot({
        req( DE_fast_raw())
        req(DE_fast_res())
        
        dt = DE_fast_raw()
        p_dt = DE_fast_res()
        
        p_dt
        high_dt = p_dt[lg2_min > 10][order(-abs(lg2_fc))][1:10][order(lg2_fc)]
        high_genes = as.character(high_dt$gene_name)
        
        sel_dt = dt[gene_name %in% high_genes]
        sel_dt$gene_name = factor(sel_dt$gene_name, levels = high_genes)
        # browser()
        ggplot(sel_dt, aes(x = gene_name, y = log2(expression+1), color = group)) +
            geom_boxplot(position = "dodge", width = .6) +
            labs(y = "log2 expression", x= "") +
            scale_color_manual(values = c("A" = "red", "B" = "blue", "." = "gray"), drop = FALSE) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    })
    
    output$plot_volcano = renderPlot({
        req(DE_fast_res())
        p_dt = DE_fast_res()
        
        p_dt
        high_dt = p_dt[lg2_min > 10][order(-abs(lg2_fc))][1:10]
        
        ggplot(p_dt, aes(x = lg2_fc, y = lg2_min)) +
            geom_point(alpha = .1) +
            geom_point(data= high_dt) +
            ggrepel::geom_text_repel(data = high_dt, aes(label = gene_name)) +
            labs(x = "log2 fold-change(B / A)", y = "log2 min(A, B)")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
