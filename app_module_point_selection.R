point_selection_plot_sizes = "400px"
ps_size = point_selection_plot_sizes

ui_point_selection = function(){
    sidebarLayout(
        sidebarPanel(
            tags$h5("TODO: clean up some weird A/B interaction, B trumps A currently."),
            disabled(selectInput("sel_A_clust", label = "A group clusters", choices = "", multiple = TRUE)),
            disabled(selectInput("sel_B_clust", label = "B group clusters", choices = "", multiple = TRUE)),
            tags$label("Patient Clustering Methods"),
            tabsetPanel(id = "tabs_cluster_method",
                        tabPanel("knn", #id = "knn", 
                                 numericInput("num_nn", label = "Nearest neighbors", value = 5, min = 2, max = Inf)
                        ), 
                        tabPanel("kmeans", #id = "kmeans", 
                                 numericInput("num_kmeans", label = "k", value = 5, min = 2, max = Inf)
                        ), 
                        tabPanel("hierarchical", #id = "hierarchical", 
                                 numericInput("num_clust", label = "n_clust", value = 5, min = 2, max = Inf)
                        )
            ),
            withSpinner(plotOutput("plot_tsne_clusters", width = ps_size, height = ps_size))
        ),
        mainPanel(
            withSpinner(plotOutput("plot_A_group", width = ps_size, height = ps_size,
                                   brush = brushOpts(id = "plot_A_brush"))),
            tags$h3("Refine selection"),
            fluidRow(
                column(width = 4,
                       actionButton("btn_add_A", "Add A"),
                       actionButton("btn_rm_A", "Remove A"),
                       actionButton("btn_limit_A", "Limit A"),
                ), 
                column(width = 4,
                       actionButton("btn_add_B", "Add B"),
                       actionButton("btn_rm_B", "Remove B"),
                       actionButton("btn_limit_B", "Limit B")
                )
            )
        )
    )
}

server_point_selection = function(input, output, session, tsne_clust, meta_data, tsne_res){
    
    output$plot_tsne_clusters = renderPlot({
        req(tsne_clust())
        tsne_dt = tsne_clust()
        tsne_dt = merge(tsne_dt, meta_data(), by = "patient_id")
        
        p = ggplot(tsne_dt, aes(x = x, y = y, color = cluster_id)) + 
            geom_point() + 
            coord_fixed() +
            labs(x = "", y = "", title = "patient clustering", subtitle = paste(input$sel_gene_list, "gene list")) +
            theme(panel.background = element_blank(),
                  panel.grid = element_blank())
        p
    })
    
    observe({
        if(is.null(tsne_res())){
            tsne_clust(NULL)
        }
        req(tsne_res())
        tsne_dt = tsne_res()
        clust_method = input$tabs_cluster_method
        showNotification(paste("clustering method is", clust_method))
        if(clust_method == "knn"){
            tsne_dt.clust = nn_clust(tsne_dt, nn = input$num_nn)
        }else if(clust_method == "kmeans"){
            tsne_dt.clust = km_clust(tsne_dt, k = input$num_kmeans)
        }else if(clust_method == "hierarchical"){
            tsne_dt.clust = h_clust(tsne_dt, n_clust = input$num_clust)
        }else{
            stop("Unrecognized clustering method, ", clust_method)
        }
        tsne_dt.clust$group = "bg"
        tsne_clust(tsne_dt.clust)
        # updateSelectInput(session, "sel_A_clust", choices = cl, selected = cl[1])
        # updateSelectInput(session, "sel_B_clust", choices = cl, selected = cl[2])
    })
    
    observeEvent({
        tsne_clust()
    },
    {
        if(is.null(tsne_clust())){
            disable("sel_A_clust")
            disable("sel_B_clust")
        }
        cl = sort(unique(tsne_clust()$cluster_id))
        cl = cl[order(as.numeric(sub("[a-zA-Z ]+", "", cl)))]
        enable("sel_A_clust")
        enable("sel_B_clust")
        updateSelectInput(session, "sel_A_clust", choices = cl, selected = cl[1])
        updateSelectInput(session, "sel_B_clust", choices = cl, selected = cl[2])
        prev_A_clust = cl[1]
        prev_B_clust = cl[2]
        sel_A_ids(tsne_clust()[cluster_id %in% cl[1]]$bcr_patient_barcode)
        sel_B_ids(tsne_clust()[cluster_id %in% cl[2]]$bcr_patient_barcode)
        
    })
    
    
    sel_A_ids = reactiveVal(character())
    sel_B_ids = reactiveVal(character())
    prev_A_clust = character()
    prev_B_clust = character()
    
    update_ids = function(current_clust, previous_clust, current_ids){
        added_clust = setdiff(current_clust, previous_clust)
        removed_clust = setdiff(previous_clust, current_clust)
        if(length(added_clust) > 0 & length(removed_clust) > 0){
            stop("inconceivable! new and missing longer than 0.")
        }
        tsne_dt = tsne_clust()
        new_ids = current_ids
        if(length(added_clust) > 0){
            new_ids = union(
                current_ids,
                tsne_dt[cluster_id %in% added_clust]$bcr_patient_barcode          
            )
        }
        if(length(removed_clust) > 0){
            new_ids = setdiff(
                current_ids,
                tsne_dt[cluster_id %in% removed_clust]$bcr_patient_barcode          
            )
        }
        new_ids
    }
    
    observeEvent({
        input$sel_A_clust
    }, {
        new_ids = update_ids(input$sel_A_clust, 
                             prev_A_clust,
                             isolate(sel_A_ids()))
        prev_A_clust = input$sel_A_clust
        sel_A_ids(new_ids)
    })
    
    observeEvent({
        input$sel_B_clust
    }, {
        new_ids = update_ids(input$sel_B_clust, 
                             prev_B_clust,
                             isolate(sel_B_ids()))
        prev_B_clust = input$sel_B_clust
        sel_B_ids(new_ids)
    })
    
    output$plot_A_group = renderPlot({
        tsne_dt = tsne_clust()
        tsne_dt$group = "."
        
        tsne_dt[bcr_patient_barcode %in% sel_A_ids(), group := "A"]
        tsne_dt[bcr_patient_barcode %in% sel_B_ids(), group := "B"]
        
        tsne_dt$group = factor(tsne_dt$group, )
        
        p = ggplot(tsne_dt, aes(x = x, y = y, color = group)) +
            geom_point(data= tsne_dt[group == "."], size = .3) +
            geom_point(data= tsne_dt[group != "."], size = .8) +
            scale_color_manual(values = c("A" = "red", "B" = "blue", "." = "gray"), drop = FALSE) +
            coord_fixed() +
            theme(panel.background = element_blank(),
                  panel.grid = element_blank()) +
            labs(x = "", y = "", title = "A and B selection")
        p
    })
    
    #plot A buttons
    observeEvent({
        input$btn_add_A
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_A_ids(
            union(isolate(sel_A_ids()), 
                  ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_rm_A
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_A_ids(
            setdiff(isolate(sel_A_ids()), 
                    ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_limit_A
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_A_ids(
            intersect(isolate(sel_A_ids()), 
                      ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_add_B
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_B_ids(
            union(isolate(sel_B_ids()), 
                  ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_rm_B
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        
        sel_B_ids(
            setdiff(isolate(sel_B_ids()), 
                    ids_in_rng)    
        )
    })
    
    observeEvent({
        input$btn_limit_B
    }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = tsne_dt[x >= brsh$xmin & x <= brsh$xmax & y >= brsh$ymin & y <= brsh$ymax]$bcr_patient_barcode
        sel_B_ids(
            intersect(isolate(sel_B_ids()), 
                      ids_in_rng)    
        )
    })
    
    # id_groups = reactiveValues(A = list(), B = list())
    
    A = reactiveVal()
    B = reactiveVal()
    
    observeEvent({
        sel_A_ids()
        sel_B_ids()
    }, {
        A(sel_A_ids)
        B(sel_B_ids)
    })
    
    list(A = A, B = B)
    
}