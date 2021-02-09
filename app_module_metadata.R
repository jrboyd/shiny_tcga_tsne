metadataServer = function(rv_metadata, Name = "metadata", var_ignore = c("patient_id", "sample_id")){
    name = tolower(Name)
    moduleServer(
        name,
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
            
            
            return(sel_metadata_var)
        }
    )
}

server_metadata = function(input, output, session, rv_metadata, Name = "metadata", var_ignore = c("patient_id", "sample_id")){
    name = tolower(Name)
    message(Name)
    sel_metadata_var = reactiveVal()
    
    observeEvent({
        input[[paste0("sel_", name)]]
    },{
        sel_metadata_var(input[[paste0("sel_", name)]])
    })
    
    output[[paste0("dynamic_sel_", name)]] = renderUI({
        # message("dynamic_sel_", name)
        req(rv_metadata())
        print(rv_metadata())
        clin_var = colnames(rv_metadata())
        clin_var = setdiff(clin_var, var_ignore)
        sapply(rv_metadata()[, clin_var, with = FALSE], class)
        selectInput(paste0("sel_", name), label = Name, choices = clin_var)
    })
    
    output[[paste0("plot_", name)]] = renderPlot({
        # message("plot_", name)
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
    
    
    return(sel_metadata_var)
}

metadataUI = function(Name = "metadata"){
    name = tolower(Name)
    ns = NS(name)
    tagList(
        uiOutput(ns("dynamic_sel_attrib")),
        shinycssloaders::withSpinner(plotOutput(ns("plot_attrib"))))
}

ui_metadata = function(name = "metadata"){
    tags$span(
        uiOutput(paste0("dynamic_sel_", name)),
        shinycssloaders::withSpinner(plotOutput(paste0("plot_", name))))
}