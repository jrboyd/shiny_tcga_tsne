
server_expression_matrix = function(input, output, session,
                                    expression_files,
                                    expression_loaded,
                                    tcga_data,
                                    clinical_loaded,
                                    meta_data,
                                    code2type,
                                    tsne_input,
                                    tsne_res){
    
    
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
    
}