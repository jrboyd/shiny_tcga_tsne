library(shiny)
library(data.table)
library(ggplot2)
library(BiocFileCache)
library(digest)
library(shinyjs)
library(shinycssloaders)

source("AppDataset.R")
source("setup_gene_lists.R")
source("setup_datasets.R")
source("app_module_expression_matrix.R")
source("app_module_upload.R")
source("app_module_tsne.R")
source("app_module_gene_xy_vis.R")
source("app_module_point_selection.R")
source("app_module_diff_expression.R")
source("functions.R")
# based on analyze_BRCA_tsne_allGenes_plusCells.R



# code2type = c("01" = "tumor", "06" = "metastasis", "11" = "normal")
code2type = code_dt$long_name
names(code2type) = code_dt$code

FACET_VAR = list(NONE = "none", SAMPLE_TYPE = "sample type", PAM50 = "PAM50")
clean_list = function(l){
    tmp = unlist(l)
    names(tmp) = NULL
    tmp
}

# Define UI for application that draws a histogram
ui2 <- fluidPage(
    theme = "bootstrap.css",
    # Application title
    useShinyjs(),
    titlePanel("development for clinical module"),
    selectInput("sel_data", label = "TCGA data", choices = dataset_names),
    uiOutput("dynamic_sel_clinical"),
    shinycssloaders::withSpinner(plotOutput("plot_attribute"))
)

server2 <- function(input, output, session) {
    ### TCGA expression data
    tcga_data = reactiveVal()
    tsne_input = reactiveVal()
    ### Genes to use in t-sne
    input_genes = reactiveVal()
    #subset of parsed genes in tcga expression data
    valid_genes = reactiveVal()
    #table shown in Add Gene Set main panel, includes parsed genes and other data used for selection
    gene_table = reactiveVal()
    #all added custom gene sets
    custom_gene_sets = reactiveVal(list())
    
    #metadata for patient entries
    meta_data = reactiveVal()
    sel_clinical_var = reactiveVal()
    #metadata for sample entries
    sample_data = reactiveVal()
    active_dataset = reactiveVal()
    
    #result of tsne
    tsne_res = reactiveVal()
    #tsne_res with clustering applied
    tsne_clust = reactiveVal()
    
    dataset_downstream = list(
        sel_clinical_var
        # sample_groups_A, 
        # sample_groups_B, 
        # DE_res, 
        # DE_fast_raw, 
        # DE_fast_res,
        # DE_res
    )
    
    server_expression_matrix(input, output, session,
                             app_datasets,
                             active_dataset,
                             dataset_downstream,
                             tcga_data,
                             meta_data,
                             sample_data,
                             code2type,
                             tsne_input,
                             tsne_res)
    
    output$dynamic_sel_clinical = renderUI({
        req(meta_data())
        clin_var = colnames(meta_data())
        clin_var = setdiff(clin_var, "patient_id")
        sapply(meta_data()[, clin_var, with = FALSE], class)
        selectInput("sel_clinical", label = "Attributes", choices = clin_var)
    })
    
    observeEvent({
        input$sel_clinical
    },{
        sel_clinical_var(input$sel_clinical)
    })
    
    output$plot_attribute = renderPlot({
        input$sel_data
        req(meta_data())
        sel_var = sel_clinical_var()
        req(sel_var)
        clin_dt = meta_data()#[, c(sel_var, "patient_id"), with = FALSE]
        if(is.numeric(clin_dt[[sel_var]])){
            ggplot(clin_dt, aes_string(x = sel_var)) +
                geom_histogram() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
        }else{
            ggplot(clin_dt, aes_string(x = sel_var)) +
                geom_bar()
        }
    })
}

# Run the application 
shinyApp(ui = ui2, server = server2)
