
server_expression_matrix = function(input, output, session,
                                    app_datasets,
                                    active_dataset,
                                    dataset_downstream,
                                    tcga_data,
                                    meta_data,
                                    sample_data,
                                    code2type,
                                    tsne_input,
                                    tsne_res
){
    observeEvent({
        input$sel_data
    }, {
        if(!isLoaded(app_datasets[[input$sel_data]])){
            app_datasets[[input$sel_data]] = LoadAppDataset(app_datasets[[input$sel_data]])
        }
        active_dataset(app_datasets[[input$sel_data]])
    })
    
    #handle change in dataset selection
    observeEvent({
        active_dataset()
    }, {
        dat = active_dataset()
        tcga_data(ExpressionMatrix(dat))
        meta_data(PatientInfo(dat))
        sample_data(SampleInfo(dat))
        
        #reset downstream
        # browser()
        for(rv in dataset_downstream){
            rv(NULL)
        }
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
        samp_dt = sample_data()
        sel_types = input$sel_sample_type_filter
        
        samp_dt = samp_dt[sample_type %in% sel_types]
        k = colnames(tcga_data()) %in% samp_dt$sample_id
        tsne_input(tcga_data()[, k])
        tsne_res(NULL)
    })
}