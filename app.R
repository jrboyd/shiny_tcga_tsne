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
bfc = BiocFileCache()
bfcif = ssvRecipes::bfcif
source("setup_gene_lists.R")
source("setup_clinical.R")
source("setup_tcga_expression.R")
source("functions.R")
# based on analyze_BRCA_tsne_allGenes_plusCells.R

code2type = c("01" = "tumor", "06" = "metastasis", "11" = "normal")
FACET_VAR = list(NONE = "none", SAMPLE_TYPE = "sample type", PAM50 = "PAM50")
clean_list = function(l){
    tmp = unlist(l)
    names(tmp) = NULL
    tmp
}

ex_files = dir("example_data", full.names = TRUE)
names(ex_files) = basename(ex_files)

run_tsne = function(expression_matrix, perplexity = 30){
    if(perplexity > ncol(expression_matrix)/4){
        warning("auto reducing perplexity")
        perplexity = round(ncol(expression_matrix)/4)
    }
    tsne_patient = bfcif(bfc, digest(list(expression_matrix, perplexity)), function(){
        Rtsne::Rtsne(t(expression_matrix), 
                     num_threads = 20, 
                     check_duplicates = FALSE,
                     perplexity = perplexity)    
    })
    
    tsne_df = as.data.table(tsne_patient$Y)
    colnames(tsne_df) = c("x", "y")
    tsne_df$bcr_patient_barcode = colnames(expression_matrix)
    tsne_df
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
                             tabPanel("Paste",
                                      textAreaInput("txt_genes", label = "Custom Genes", value = "paste genes here")#,
                             ),
                             tabPanel("Upload",
                                      fileInput(inputId = "BtnUploadPeakfile", label = "Browse Local Files"),
                                      actionButton("btnExampleData", label = "use example data"),
                                      selectInput("selExampleData", label = "select example data", choices = ex_files)    
                                      
                             )
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
    meta_data = reactiveVal()
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
    
    ### Genes to use in t-sne
    input_genes = reactiveVal()
    parsed_genes = reactiveVal()
    valid_genes = reactiveVal()
    
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
    
    gene_table = reactiveVal()
    
    observe({
        gl = parsed_genes()
        if(is.null(gl)){
            gene_table(data.frame())  
        }else if(length(gl) == 0){
            gene_table(data.frame())  
        }else{
            gene_table(data.frame(gene_name = gl))
        }
    })
    
    observe({
        df = PreviewSet_DataFrame()
        if(is.null(df)){
            gene_table(data.frame())  
        }else if(nrow(df) == 0){
            gene_table(data.frame())
        }else{
            gene_table(df)
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
    
    ### Running t-sne
    tsne_res = reactiveVal()
    observeEvent({
        # tcga_data()
        tsne_input()
        valid_genes()
    }, {
        expr_mat = tsne_input()
        gl = valid_genes()
        if(length(gl) > 0){
            tsne_worked = tryCatch({
                tsne_dt = run_tsne(expr_mat[gl,])    
                TRUE
            }, error = function(e){
                
                FALSE
            })
            if(tsne_worked){
                regx = regexpr("^.{4}-.{2}-.{4}", tsne_dt$bcr_patient_barcode)
                tsne_dt$submitter_id = regmatches(tsne_dt$bcr_patient_barcode, regx)
                tsne_dt[, sample_code := tstrsplit(bcr_patient_barcode, "-", keep = 4)]
                tsne_dt[, sample_code := sub("[A-Z]", "", sample_code)]
                tsne_dt[, sample_type := code2type[sample_code]]
                # tsne_dt$sample_code %>% table
                # tsne_dt$sample_type %>% table
                
                tsne_dt[, x := scales::rescale(x, c(-.5, .5))]
                tsne_dt[, y := scales::rescale(y, c(-.5, .5))]
                
                tsne_res(tsne_dt) 
            }else{
                showNotification("Need more valid genes/samples to run t-sne.", type = "error")
                tsne_res(NULL)
            }
            
        }else{
            tsne_res(NULL)
        }
    })
    
    ### Plot t-sne
    output$plot_tsne <- renderPlot({
        req(tsne_res())
        req(input$sel_facet_var)
        tsne_dt = tsne_res()
        
        tsne_dt = merge(tsne_dt, meta_data(), by = "submitter_id")
        # browser()
        if(input$sel_color_by == "sample type"){ #c("sample type", "PAM50")
            p = ggplot(tsne_dt, aes(x = x, y = y, color = sample_type)) 
        }else if(input$sel_color_by == "PAM50"){
            p = ggplot(tsne_dt, aes(x = x, y = y, color = pam_call)) 
        }else{
            stop("unrecognized input$sel_color_by: ", input$sel_color_by)
        }
        if(input$sel_facet_var != FACET_VAR$NONE){
            p = p + annotate("point", x= tsne_dt$x, y = tsne_dt$y, color = 'gray70', size = .3)
        }
        p = p + 
            geom_point() + 
            coord_fixed() +
            labs(x = "", y = "", title = "t-sne of TCGA samples", subtitle = paste(input$sel_gene_list, "gene list"))
        if(input$sel_facet_var == FACET_VAR$NONE){
            p
        }else if(input$sel_facet_var == FACET_VAR$SAMPLE_TYPE){
            p + facet_wrap(~sample_type)
        }else if(input$sel_facet_var == FACET_VAR$PAM50){
            p + facet_wrap(~pam_call)
        }else{
            stop("unrecognized input$sel_facet_var: ", input$sel_facet_var)
        }
    })
    
    output$plot_tsne_gene = renderPlot({
        req(tsne_res())
        tsne_dt = tsne_res()
        req(vis_gene())
        # browser()
        gene_vals = tcga_data()[vis_gene(),]
        tsne_dt$gene_val = gene_vals[tsne_dt$bcr_patient_barcode]
        ggplot(tsne_dt, aes(x = x, y = y, color = log10(gene_val + 1))) + 
            geom_point() + 
            coord_fixed() +
            labs(x = "", y = "", title = paste(vis_gene(), "expression"), subtitle = "log10 scale") +
            scale_color_viridis_c()
        
    })
    
    ### File Upload
    #Set Preview reactives
    PreviewSet_Filepath = reactiveVal(value = "", label = "PreviewSet_Filepath")
    PreviewSet_Name = reactiveVal(value = "", label = "PreviewSet_Name")
    PreviewSet_DataFrame = reactiveVal()
    # PreviewSet_DataFrame = reactive({
    #     showNotification(PreviewSet_Filepath())
    #     if(is.null(PreviewSet_Filepath())) return(NULL)
    #     if(PreviewSet_Filepath() == "") return(NULL)
    #     browser()
    #     load_peak_wValidation(PreviewSet_Filepath(), with_notes = T)
    # })
    observeEvent({
        PreviewSet_Filepath()
        PreviewSet_Name()
    }, {
        if(PreviewSet_Filepath() == ""){
            PreviewSet_DataFrame(NULL)
        }else{
            # browser()
            # fread(PreviewSet_Filepath())
            ##TODO load a DESEQ file or other gene list
            showNotification(paste0("loading ", PreviewSet_Name()))
            out = decide_parse_FUN(PreviewSet_Filepath(), PreviewSet_Name())
            PreviewSet_DataFrame(out)
        }
    })
    
    observeEvent(input$BtnUploadPeakfile, {
        PreviewSet_Filepath(input$BtnUploadPeakfile$datapath)
        PreviewSet_Name(input$BtnUploadPeakfile$name)
    })
    
    observeEvent(input$BtnCancelFile, {
        PreviewSet_Filepath("")
        PreviewSet_Name("")
    })
    
    observeEvent(input$btnExampleData,
                 {
                     PreviewSet_Filepath(input$selExampleData)
                     PreviewSet_Name(names(ex_files)[which(input$selExampleData == ex_files)])
                 })
    
    observeEvent(
        PreviewSet_DataFrame(),
        {
            req(PreviewSet_DataFrame())
            # showModal(modalDialog( size = "l",
            #     DT::dataTableOutput("DT_UploadGenes_DataFrame")
            #     
            # ))
        }
    )
    
    output$DT_UploadGenes_DataFrame = DT::renderDataTable({
        DT::datatable(PreviewSet_DataFrame()[[1]], options = list(scrollX = TRUE))  
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
