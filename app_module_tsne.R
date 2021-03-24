library(Rtsne)
bfc = BiocFileCache()
bfcif = ssvRecipes::bfcif

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
    tsne_df$sample_id = colnames(expression_matrix)
    tsne_df
}

server_tsne = function(input, output, session, tsne_res, tsne_input, valid_genes, meta_data, sample_data, code2type, FACET_VAR){
    ### Running t-sne
    
    observeEvent({
        # tcga_data()
        tsne_input()
        valid_genes()
    }, {
        expr_mat = tsne_input()
        gl = valid_genes()
        
        #choose dimensional reduction method, if possible
        if(ncol(expr_mat) < 3){
            showNotification("Too few samples for dimensional reduction.", type = "error")
            tsne_worked = FALSE
        }else if(ncol(expr_mat) < 20){
            browser()
            pc = prcomp(expr_mat[gl,])
            tsne_dt = as.data.table(pc$rotation[,1:2], keep.rownames = TRUE)[, c(2:3, 1)]
            setnames(tsne_dt, c("rn", "PC1", "PC2"), c("sample_id", "x", "y"))
            tsne_worked = TRUE
        }else if(length(gl) > 0){
            tsne_worked = tryCatch({
                tsne_dt = run_tsne(expr_mat[gl,])    
                TRUE
            }, error = function(e){
                FALSE
            })
        }else{
            tsne_worked = FALSE
        }
        
        if(tsne_worked){
            samp_dt = sample_data()
            tsne_dt = merge(tsne_dt, samp_dt, by = "sample_id")
            tsne_dt[, x := scales::rescale(x, c(-.5, .5))]
            tsne_dt[, y := scales::rescale(y, c(-.5, .5))]
            tsne_res(tsne_dt) 
            
        }else{
            showNotification("Need more valid genes/samples to run t-sne.", type = "error")
            tsne_res(NULL)
        }
    })
    
    ### Plot t-sne
    output$plot_tsne <- renderPlot({
        req(tsne_res())
        req(input$sel_facet_var)
        tsne_dt = tsne_res()
        tsne_dt = merge(tsne_dt, meta_data(), by = "patient_id")
        # browser()
        # if(input$sel_color_by == "sample type"){ #c("sample type", "PAM50")
        #     p = ggplot(tsne_dt, aes(x = x, y = y, color = sample_type)) 
        # }else if(input$sel_color_by == "PAM50"){
        #     p = ggplot(tsne_dt, aes(x = x, y = y, color = pam_call)) 
        # }else{
        #     stop("unrecognized input$sel_color_by: ", input$sel_color_by)
        # }
        if(input$sel_color_by ==  FACET_VAR$NONE){
            p = ggplot(tsne_dt, aes_string(x = "x", y = "y")) 
        }else{
            p = ggplot(tsne_dt, aes_string(x = "x", y = "y", color = input$sel_color_by))     
        }
        
        if(input$sel_facet_var != FACET_VAR$NONE){
            p = p + annotate("point", x= tsne_dt$x, y = tsne_dt$y, color = 'gray70', size = .3)
        }
        p = p + 
            geom_point() + 
            coord_fixed() +
            labs(x = "", y = "", title = "t-sne of TCGA samples", subtitle = paste(input$sel_gene_list, "gene list")) +
            theme(panel.background = element_blank(),
                  panel.grid = element_blank())
        if(input$sel_facet_var != FACET_VAR$NONE){
            p = p + facet_wrap(paste0("~", input$sel_facet_var))
        }
        # }else if(input$sel_facet_var == FACET_VAR$SAMPLE_TYPE){
        #     p + facet_wrap(~sample_type)
        # }else if(input$sel_facet_var == FACET_VAR$PAM50){
        #     p + facet_wrap(~pam_call)
        # }else{
        #     stop("unrecognized input$sel_facet_var: ", input$sel_facet_var)
        # }
        p
    })
    
    
    
    tsne_res
}