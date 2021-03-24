


check_AppDataset = function(object){
    probs = c()
    #at a minimum, files should exist
    if(!file.exists(object@expression_file)){
        probs = c(probs, paste("expression_file does not exist:", object@expression_file))
    }
    if(!file.exists(object@sample_file)){
        probs = c(probs, paste("sample_file does not exist:", object@sample_file))
    }
    if(!file.exists(object@clinical_file)){
        probs = c(probs, paste("clinical_file does not exist:", object@clinical_file))
    }
    if(!file.exists(object@info_file)){
        probs = c(probs, paste("info_file does not exist:", object@info_file))
    }
    if(!is.null(object@attribute_dir))
        if(!dir.exists(object@attribute_dir)){
            probs = c(probs, paste("attribute_dir does not exist:", object@attribute_dir))
        }
    
    #additional checks on loaded data
    if(object@is_loaded){
        if(!all(colnames(object@expression_data) %in% object@sample_data$sample_id)){
            probs = c(probs, "samples missing from sample info that are present in expression data.")
            probs = c(probs, paste(setdiff(colnames(object@expression_data), object@sample_data$sample_id), collapse = ", "))
        }
        if(!all(object@sample_data$sample_id %in% colnames(object@expression_data))){
            probs = c(probs, "samples missing from expression data that are present in sample info.")
            probs = c(probs, paste(setdiff(object@sample_data$sample_id, colnames(object@expression_data)), collapse = ", "))
        }
        if(!all(object@sample_data$patient_id %in% object@clinical_data$patient_id)){
            probs = c(probs, "patients missing from clinical info that are present in sample data.")
            probs = c(probs, paste(setdiff(object@sample_data$patient_id, object@clinical_data$patient_id), collapse = ", "))
        }
        if(!all(object@clinical_data$patient_id %in% object@sample_data$patient_id)){
            probs = c(probs, "patients missing from sample data that are present in clinical info.")
            probs = c(probs, paste(setdiff(object@clinical_data$patient_id, object@sample_data$patient_id), collapse = ", "))
        }
    }
    
    if(length(probs) == 0){
        return(TRUE)
    }else{
        smax = 80
        probs = substr(probs, 0, smax)
        probs[nchar(probs) >= smax] = paste0(probs[nchar(probs) >= smax], "...etc")
        return(probs)
    }
}

setClass("AppDataset", representation(
    name = "character",
    #files
    expression_file = "character", 
    sample_file = "character",
    clinical_file = "character", 
    info_file = "character",
    attribute_dir = "character",
    #raw data
    expression_data = "matrix",
    sample_data = "data.table",
    clinical_data = "data.table",
    info_data = "data.table",
    #processed data
    expression_processed = "matrix",
    sample_processed = "data.table",
    clinical_processed = "data.table",
    #color mappings
    clinical_color_mappings = "list",
    sample_color_mappings = "list",
    default_color_function = "function",
    
    is_loaded = "logical"), 
    validity = check_AppDataset)

setMethod("show", "AppDataset", function(object){
    message("This is an AppDataset named ", object@name, ".")
    if(object@is_loaded){
        message("It has been loaded.")
    }else{
        message("It has NOT been loaded.")
        message("Use LoadAppDataset(this) to Load the relevevant data.")
    }
})

# object = LoadAppDataset(object)
# probs = check_AppDataset(object)

CreateAppDataset = function(name, expression_file, sample_file, clinical_file, info_file, attribute_dir = NULL){
    ad = new("AppDataset", 
             name = name,
             expression_file = expression_file, 
             sample_file = sample_file, 
             clinical_file = clinical_file, 
             info_file = info_file,
             attribute_dir = attribute_dir,
             expression_data = matrix(),
             sample_data = data.table(),
             clinical_data = data.table(),
             info_data = data.table(),
             sample_processed = data.table(),
             clinical_processed = data.table(),
             is_loaded = FALSE)
    ad@default_color_function = function(items){
        if(is.factor(items)){
            item_names = levels(items)    
        }else if(is.character(items)){
            item_names = unique(items)
        }else{
            stop("items must be character or factor")
        }
        
        n = length(item_names)
        stopifnot(is.numeric(n))
        pal = "Dark2"
        stopifnot(is.character(pal))
        if (n < 1) 
            stop("n must be at least 1")
        
        pal_info = RColorBrewer::brewer.pal.info
        pal_info$brewName = rownames(pal_info)
        rownames(pal_info) = tolower(rownames(pal_info))
        pal = tolower(pal)
        if (!any(pal == rownames(pal_info))) 
            stop("Palette", pal, "not a valid RColorBrewer palette, ", 
                 "see RColorBrewer::brewer.pal.info")
        maxColors = pal_info[pal, ]$maxcolors
        nbrew = min(max(n, 3), maxColors)
        cols = RColorBrewer::brewer.pal(n = nbrew, name = pal_info[pal, ]$brewName)[(seq_len(n) - 1)%%maxColors + 1]
        names(cols) = item_names
        cols
    }
    ad
}

apply_lev_file = function(appDat, lev_file){
    lev_root = basename(lev_file)
    if(!appDat@is_loaded){
        stop("appDat should be loaded prior to applying level files.")
    }
    if(grepl("^clinical", lev_root)){
        data_slot = "clinical_data"
    }else if(grepl("^sample", lev_root)){
        data_slot = "sample_data"
    }else{
        stop("could not determine slot names from file name")
    }
    df = slot(appDat, data_slot)
    var = sub(".levels", "", sub("(^clinical_)|(^sample_)", "", lev_root))
    if(!var %in% colnames(df)){
        stop("variable ", var, " was not found in ", data_slot)
    }
    lev = readLines(lev_file)
    if(!all(df[[var]] %in% lev)){
        missed = setdiff(df[[var]], lev)
        if(length(missed) > 1 | missed[1] != "")
            warning(paste(collapse = "\n", c("all levels of loaded data were not present in provided levels. missing:", missed)))
        lev = c(lev, missed)
    }
    df[[var]] = factor(df[[var]], levels = lev)
    slot(appDat, data_slot) = df
    appDat
}

apply_color_file = function(appDat, color_file){
    color_root = basename(color_file)
    if(!appDat@is_loaded){
        stop("appDat should be loaded prior to applying color files.")
    }
    if(grepl("^clinical", color_root)){
        data_slot = "clinical_data"
        map_slot = "clinical_color_mappings"
    }else if(grepl("^sample", color_root)){
        data_slot = "sample_data"
        map_slot = "sample_color_mappings"
    }else{
        stop("could not determine slot names from file name")
    }
    df = slot(appDat, data_slot)
    var = sub(".colors", "", sub("(^clinical_)|(^sample_)", "", color_root))
    if(!var %in% colnames(df)){
        stop("variable ", var, " was not found in ", data_slot)
    }
    color_tbl = read.table(color_file, col.names = c("value", "color"), sep = ",")
    colors = color_tbl$color
    color_names = color_tbl$value
    if(!all(df[[var]] %in% color_names)){
        missed = setdiff(df[[var]], color_names)
        if(length(missed) > 1 | missed[1] != "")
            warning(paste(collapse = "\n", c("all values of loaded data were not present in provided colors mapping. missing:", missed)))
        colors = c(colors, gray.colors(length(missed)))
        color_names = c(color_names, missed)
    }
    names(colors) = color_names
    new_col_list = list(colors)
    names(new_col_list) = var
    old_col_list = slot(appDat, map_slot)
    #remove possible duplicate
    old_col_list = old_col_list[!names(old_col_list) %in% names(new_col_list)]
    slot(appDat, map_slot) = c(old_col_list, new_col_list)
    appDat  
}

CreateAppDataset.installed = function(install_dir){
    ad = CreateAppDataset(
        name = basename(install_dir),
        file.path(install_dir, "expression.csv"),
        file.path(install_dir, "samples.csv"),
        file.path(install_dir, "clinical.csv"),
        file.path(install_dir, "info.txt"),
        file.path(install_dir, "attribute_info"))
    ad
}

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
    
    if(!is.null(object@attribute_dir)){
        lev_files = dir(object@attribute_dir, pattern = "levels$", full.names = TRUE)
        for(f in lev_files){
            object = apply_lev_file(object, f)
        }
        col_files = dir(object@attribute_dir, pattern = "colors$", full.names = TRUE)
        for(f in col_files){
            object = apply_color_file(object, f)
        }
    }
    
    object = InitLoadedAppDataset(object)
    object
}

#reduce sample info to expression sample_id, and clinical to sample patient_id
#error if sample_id missing from sample, or patient_id missing from clinical
InitLoadedAppDataset = function(object){
    object@sample_data = object@sample_data[object@sample_data$sample_id %in% colnames(object@expression_data),]
    object@clinical_data = object@clinical_data[object@clinical_data$patient_id %in% object@sample_data$patient_id,]
    object
}

UnloadAppDataset = function(object){
    object@expression_data = matrix()
    object@sample_data = data.table()
    object@clinical_data = data.table()
    object@info_data = data.table()
    object@is_loaded = FALSE
    object
}

object = CreateAppDataset.installed("installed_datasets/TARGET")
object = CreateAppDataset.installed("installed_datasets/BRCA_tiny/")
object

check_AppDataset(object)

object = LoadAppDataset(object)
object

object = UnloadAppDataset(object)
object

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

isLoaded = function(object){
    object@is_loaded
}

ValidateInstallDataset = function(dest){
    object = CreateAppDataset.installed(dest)
    LoadAppDataset(object)
    check_AppDataset(object)
}

InstallDataset.describe = function(){
    message("Clinical data must have 'patient_id'")
    message("Clinical data supports survival analysis if it has 'days_to_death', 'days_to_last_follow_up', and 'vital_status'")
    message("Sample data must have 'patient_id', 'sample_id', 'sample_code', 'sample_type', and 'sample_type_short'")
    message("Expression data has colnames that are 'sample_id' in sample data;  rownames are unique gene ids or names.")
}

InstallDataset = function(dataset_dirname, 
                          expression_data, 
                          clinical_data, 
                          sample_data, 
                          dataset_name = dataset_dirname, 
                          install_location = "installed_datasets", 
                          survival_ids = c("days_to_death", "days_to_last_follow_up", "vital_status")){
    #validity checks
    if(colnames(expression_data)[1] != "gene_name"){
        stop("column 1 of expression_data must be gene_name and should contain Gene Symbols.")
    }
    req_clinical = "patient_id"
    if(!all(req_clinical %in% colnames(clinical_data))){
        stop(paste(setdiff(req_clinical, colnames(clinical_data)), collapse = ", "), " are missing from clinical_data")
    }
    req_sample = c("sample_code", "sample_id", "patient_id", "sample_type", "sample_type_short")
    if(!all(req_sample %in% colnames(sample_data))){
        stop(paste(setdiff(req_sample, colnames(sample_data)), collapse = ", "), " are missing from sample_data")
    }
    
    if(any(duplicated(clinical_data$patient_id))) stop("duplicated patient_id in clinical_data")
    if(any(duplicated(sample_data$sample_id))) stop("duplicated sample_id in sample_data")
    if(!all(colnames(expression_data)[-1] %in% sample_data$sample_id)) stop("some sample info missing for expression_data")
    if(!all(sample_data$patient_id %in% clinical_data$patient_id)) stop("some patient info missing for sample_data")
    
    if(any(duplicated(clinical_data$patient_id))) stop("patient_id not all unique in clinical_data")
    if(any(duplicated(sample_data$sample_id))) stop("sample_id not all unique in sample_data")
    
    factors = list()
    #column types
    for(var in colnames(clinical_data)){
        message(var)
        if(is.factor(clinical_data[[var]])){
            k = levels(clinical_data[[var]]) == "--"
            levels(clinical_data[[var]])[k] = NA
            factors[[paste0("clinical_",var)]] = levels(clinical_data[[var]])
        }else{
            k = clinical_data[[var]] == "--"
            k[is.na(k)] = FALSE
            if(any(k)){
                clinical_data[[var]][k] = NA
                as_num = suppressWarnings({
                    as.numeric(clinical_data[[var]])
                })
                if(all(!is.na(as_num[!k]))){
                    clinical_data[[var]] = as_num
                }
            }
        }
    }
    
    for(var in colnames(sample_data)){
        message(var)
        if(is.factor(sample_data[[var]])){
            k = levels(sample_data[[var]]) == "--"
            levels(sample_data[[var]])[k] = NA
            factors[[paste0("sample_",var)]] = levels(sample_data[[var]])
        }else{
            k = sample_data[[var]] == "--"
            k[is.na(k)] = FALSE
            if(any(k)){
                sample_data[[var]][k] = NA
                as_num = suppressWarnings({
                    as.numeric(sample_data[[var]])
                })
                if(all(!is.na(as_num[!k]))){
                    sample_data[[var]] = as_num
                }
            }
        }
    }
    
    
    
    #setup dir
    dest = file.path(install_location, dataset_dirname)
    dir.create(dest, showWarnings = FALSE, recursive = TRUE)
    #write files
    fwrite(expression_data, file.path(dest, "expression.csv"))
    fwrite(clinical_data, file.path(dest, "clinical.csv"))
    fwrite(sample_data, file.path(dest, "samples.csv"))
    
    att_dir = file.path(dest, "attribute_info")
    dir.create(att_dir, showWarnings = FALSE)
    for(var in names(factors)){
        lev = factors[[var]]
        writeLines(lev, file.path(att_dir, paste0(var, ".levels")))
    }
    
    #determine info
    surv_ready = all(survival_ids %in% colnames(clinical_data))
    attrib = rbind(
        paste0("name:", dataset_name),
        paste0("survival_ready:", ifelse(surv_ready, "yes", "no"))
    )
    writeLines(attrib, file.path(dest, "info.txt"))
    
    ValidateInstallDataset(dest)
    
}

appDat = CreateAppDataset.installed("installed_datasets/BRCA_tiny/")
appDat = LoadAppDataset(appDat)

as_id = function(Name){
    gsub(" ", "_", tolower(Name))
}

UI.GenClinicalFilter = function(appDat, Name = "clinical_filter"){
    id = as_id(Name)
    ns = NS(id)
    cn = setdiff(colnames(SampleInfo(appDat)), "sample_id")
    tagList(
        selectizeInput(label = "Clinical Filter", ns("sel_var"), choices = cn),
        uiOutput(ns("filter_ui"))
    )
    
}

Server.GenClinicalFilter = function(appDat, Name = "clinical_filter"){
    id = as_id(Name)
    moduleServer(
        id,
        function(input, output, session){
            
            # observeEvent({
            #     input$sel_var
            # }, {
            #     
            # })
            
            output$filter_ui = renderUI({
                var = input$sel_var
                if(is.null(var)){
                    msg = "waiting"
                }else{
                    msg = paste("UI for", var)
                }
                tags$h3(msg)
            })
            
        }
    )
}

# UI.GenSampleFilter
# 
# UI.GenClinicalTransform
# UI.GenSampleTransform
# 
# UI.GenSampleAddExpression


