
ex_files = dir("example_data", full.names = TRUE)
names(ex_files) = basename(ex_files)

ui_upload = function(){
    
    sidebarLayout(
        sidebarPanel(

            fluidRow(
                column(textInput("txt_gene_set_name", label = "Gene set name", value = "custom"), width = 9),
                column(actionButton("btn_add_gene_set", label = "Add"), width = 3, style = "margin-top: 24px;",)
            ),
            tags$hr(),
            tabsetPanel(
                id = "gene_list_method",
                tabPanel("Paste", value = "paste",
                         textAreaInput("txt_genes", label = "Custom Genes", value = "paste genes here")#,
                ),
                tabPanel("Upload", value = "upload",
                         fileInput(inputId = "BtnUploadPeakfile", label = "Browse Local Files"),
                         actionButton("btnExampleData", label = "use example data"),
                         selectInput("selExampleData", label = "select example data", choices = ex_files)    
                         
                )
                
            ),
            uiOutput("rpt_custom_gene_sets"),
            tags$h5("TODO: gene set import/export"),
            tags$h5("TODO: filtering after upload"),
            tags$h5("TODO: hide example data")
        ),
        mainPanel(
            DT::dataTableOutput("DT_PasteGenes_DataFrame")
        )
    )
    
}

server_upload = function(input, output, session, gene_table, tcga_data, custom_gene_sets){
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
        }
    )
    
    observe({
        if(input$gene_list_method == "upload"){
            showNotification("DEBUG list by upload")
            df = PreviewSet_DataFrame()[[1]]
            if(is.null(df)){
                gene_table(data.frame())  
            }else if(nrow(df) == 0){
                gene_table(data.frame())
            }else{
                gene_table(df)
            }
        }
    })
    
    ### Pasting genes
    #genes parsed for paste/upload
    parsed_genes = reactiveVal()
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
    display_table = reactiveVal()
    observeEvent({
        gene_table()
    }, {
        df = gene_table()
        df = as.data.frame(df)
        valid_genes = rownames(tcga_data())
        if(nrow(df) > 0){
            df = locate_genes_in_df(df, valid_genes)
            df = validate_genes_in_df(df, valid_genes)    
        }else{
            df = data.table(A = "please upload data")
        }
        display_table(df)
    })
    
    ### Show table
    output$DT_PasteGenes_DataFrame = DT::renderDataTable({
        DT::datatable(display_table())    
    })
    
    ### Report valid
    
    ### Add gene set
    observeEvent({
        input$btn_add_gene_set
    }, {
        current = custom_gene_sets()
        df = display_table()
        gl = subset(df, valid == TRUE)$gene_name
        current[[input$txt_gene_set_name]] = gl
        custom_gene_sets(current)
    })
    
    # observeEvent({
    #     custom_gene_sets()
    # }, {
    output$rpt_custom_gene_sets = renderUI({
        gls =  custom_gene_sets()
        hl = lapply(names(gls), function(nam){
            x = gls[[nam]]
            tags$tr(
                tags$th(nam),
                tags$th(length(x))
            )
        })
        tags$table(
            tags$tr(
                tags$th("Name"),
                tags$th("Length")
            ),
            hl
        )
    })
    # })
}