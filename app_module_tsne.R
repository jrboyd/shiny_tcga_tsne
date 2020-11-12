
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
    tsne_df$bcr_patient_barcode = colnames(expression_matrix)
    tsne_df
}
 
server_tsne = function(input, output, session, tsne_res, tsne_input, valid_genes, meta_data, code2type, FACET_VAR){
    ### Running t-sne
    
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
            labs(x = "", y = "", title = "t-sne of TCGA samples", subtitle = paste(input$sel_gene_list, "gene list")) +
            theme(panel.background = element_blank(),
                  panel.grid = element_blank())
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
    

    
    tsne_res
}