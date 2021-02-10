as_id = function(Name){
    gsub(" ", "_", tolower(Name))
}

metadataServer = function(rv_metadata, Name = "metadata", var_ignore = c("patient_id", "sample_id")){
    id = as_id(Name)
    moduleServer(
        id,
        function(input, output, session){
            sel_metadata_var = reactiveVal()
            
            observeEvent({
                input$sel_attrib
            },{
                sel_metadata_var(input$sel_attrib)
            })
            
            output$dynamic_sel_attrib = renderUI({
                ns <- session$ns
                req(rv_metadata())
                print(rv_metadata())
                clin_var = colnames(rv_metadata())
                clin_var = setdiff(clin_var, var_ignore)
                sapply(rv_metadata()[, clin_var, with = FALSE], class)
                selectInput(ns("sel_attrib"), label = Name, choices = clin_var)
            })
            
            output$plot_attrib = renderPlot({
                input$sel_data
                req(rv_metadata())
                sel_var = sel_metadata_var()
                req(sel_var)
                clin_dt = rv_metadata()
                if(is.numeric(clin_dt[[sel_var]])){
                    ggplot(clin_dt, aes_string(x = sel_var)) +
                        geom_histogram() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
                }else{
                    ggplot(clin_dt, aes_string(x = sel_var)) +
                        geom_bar() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
                }
            })
            
            output$ui_filter = renderUI({
                input$sel_data
                req(rv_metadata())
                sel_var = sel_metadata_var()
                req(sel_var)
                clin_dt = rv_metadata()
                if(is.numeric(clin_dt[[sel_var]])){
                    numericFilterUI("num_filter")
                }else{
                    characterFilterUI("char_filter")
                }
            })
            
            
            return(sel_metadata_var)
        }
    )
}

metadataUI = function(Name = "metadata"){
    id = as_id(Name)
    ns = NS(id)
    tagList(
        tabsetPanel(
            tabPanel("Overview", 
                     uiOutput(ns("dynamic_sel_attrib")),
                     shinycssloaders::withSpinner(plotOutput(ns("plot_attrib")))
            ),
            tabPanel("Filter", uiOutput(ns("ui_filter"))),
            tabPanel("Group"),
            tabPanel("Transform")
        ),
        
    )
    
}

numericFilterServer = function(Name = "Numeric Filter", rv_input_table, rv_var, rv_output_table){
    id = as_id(Name)
    moduleServer(
        id,
        function(input, output, session){
        }
    )
}

numericFilterUI = function(Name = "Numeric Filter"){
    id = as_id(Name)
    ns = NS(id)
    actionButton(ns("NumNYI"), label = "NumNYI")
    
}

characterFilterServer = function(Name = "Character Filter"){
    id = as_id(Name)
    moduleServer(
        id,
        function(input, output, session){
        }
    )
}
characterFilterUI = function(Name = "Character Filter"){
    id = as_id(Name)
    ns = NS(id)
    actionButton(ns("CharNYI"), label = "CharNYI")
}
