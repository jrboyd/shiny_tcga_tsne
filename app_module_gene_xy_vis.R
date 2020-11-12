server_gene_xy = function(input, output, session, xy_df, color_df, rownames_to_vis){
    
    output$plot_tsne_gene = renderPlot({
        req(xy_df())
        xy = xy_df()
        req(rownames_to_vis())
        # browser()
        color_vals = color_df()[rownames_to_vis(),]
        xy$color_val = color_vals[xy$bcr_patient_barcode]
        ggplot(xy, aes(x = x, y = y, color = log10(color_val + 1))) + 
            geom_point() + 
            coord_fixed() +
            labs(x = "", y = "", title = paste(rownames_to_vis(), "expression"), subtitle = "log10 scale") +
            scale_color_viridis_c() +
            theme(panel.background = element_blank(),
                  panel.grid = element_blank())
        
    })
    
}