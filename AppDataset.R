
setClass("AppDataset", representation(
    name = "character",
    expression_file = "character", 
    sample_file = "character",
    clinical_file = "character", 
    info_file = "character",
    expression_data = "matrix",
    sample_data = "data.frame",
    clinical_data = "data.frame",
    info_data = "data.frame",
    is_loaded = "logical"))

CreateAppDataset = function(name, expression_file, sample_file, clinical_file, info_file){
    new("AppDataset", 
        name = name,
        expression_file = expression_file, 
        sample_file = sample_file, 
        clinical_file = clinical_file, 
        info_file = info_file,
        # expression_data = data.table(), 
        # sample_data = data.table(),
        # clinical_data = data.table(),
        # info_data = data.table(),
        is_loaded = FALSE)
}

CreateAppDataset.installed = function(install_dir){
    CreateAppDataset(
        name = basename(install_dir),
        file.path(install_dir, "expression.csv"),
        file.path(install_dir, "samples.csv"),
        file.path(install_dir, "clinical.csv"),
        file.path(install_dir, "info.txt"))
}

setMethod("show", "AppDataset", function(object){
    message("This is an AppDataset named ", object@name, ".")
    if(object@is_loaded){
        message("It has been loaded.")
    }else{
        message("It has NOT been loaded.")
        message("Use LoadAppDataset(this) to Load the relevevant data.")
    }
})

LoadAppDataset = function(object){
    dat = data.table::fread(object@expression_file)
    mat = as.matrix(dat[,-1])
    rownames(mat) = dat$gene_name
    dat = mat
    exp_dt = data.table::as.data.table(dat)
    exp_dt$gene_name = rownames(dat)
    exp_dt = data.table::melt(exp_dt, id.vars = "gene_name")
    exp_dt = exp_dt[, .(value = max(value)), .(gene_name, variable)]
    # exp_dt = unique(exp_mat$gene_name[duplicated(exp_mat$gene_name)])
    exp_dt = data.table::dcast(exp_dt, gene_name~variable, value.var = "value")
    exp_mat = as.matrix(exp_dt[, -1])
    rownames(exp_mat) = exp_dt$gene_name
    stopifnot(!any(duplicated(rownames(exp_mat))))
    
    object@expression_data = exp_mat
    object@sample_data = data.table::fread(object@sample_file)
    object@clinical_data = data.table::fread(object@clinical_file)
    object@info_data = data.table::fread(object@info_file)
    object@is_loaded = TRUE
    object
}

object = CreateAppDataset.installed("installed_datasets/BRCA_tiny")
object = LoadAppDataset(object)

ExpressionMatrix = function(object){
    object@expression_data    
}
SampleInfo = function(object){
    object@sample_data    
}
PatientInfo = function(object){
    object@clinical_data    
}
DatasetName = function(object){
    object@name
}

