


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
    expression_file = "character", 
    sample_file = "character",
    clinical_file = "character", 
    info_file = "character",
    expression_data = "matrix",
    sample_data = "data.frame",
    clinical_data = "data.frame",
    info_data = "data.frame",
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
    object = InitLoadedAppDataset(object)
    object
}

#reduce sample info to expression sample_id, and clinical to sample patient_id
#error if sample_id missing from sample, or patient_id missing from clinical
InitLoadedAppDataset = function(object){
    object@sample_data = object@sample_data[sample_id %in% colnames(object@expression_data)]
    object@clinical_data = object@clinical_data[patient_id %in% object@sample_data$patient_id]
    object
}

UnloadAppDataset = function(object){
    object@expression_data = matrix()
    object@sample_data = data.frame()
    object@clinical_data = data.frame()
    object@info_data = data.frame()
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
    
    #column types
    for(var in colnames(clinical_data)){
        message(var)
        k = clinical_data[[var]] == "--"
        k[is.na(k)] = FALSE
        any(k)
        which(k)
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
    
    #setup dir
    dest = file.path(install_location, dataset_dirname)
    dir.create(dest, showWarnings = FALSE, recursive = TRUE)
    # browser()
    #write files
    fwrite(expression_data, file.path(dest, "expression.csv"))
    fwrite(clinical_data, file.path(dest, "clinical.csv"))
    fwrite(sample_data, file.path(dest, "samples.csv"))
    
    #determine info
    surv_ready = all(survival_ids %in% colnames(clinical_data))
    attrib = rbind(
        paste0("name:", dataset_name),
        paste0("survival_ready:", ifelse(surv_ready, "yes", "no"))
    )
    writeLines(attrib, file.path(dest, "info.txt"))
    
    ValidateInstallDataset(dest)
    
}



#this installs TCGA and TARGET datasets
if(FALSE){
    {
        library(data.table)
        code_dt = fread("sample_codes.txt", colClasses = rep("character", 3), col.names = c("sample_code", "sample_type", "sample_type_short"))
        
        expression_files = list(
            BRCA_tiny = "data/BRCA_TCGA_expression.tiny.csv",
            BRCA = "data/BRCA_TCGA_expression.csv"
        )
        
        expression_loaded = lapply(expression_files, fread)
        
        clinical_files = list(
            BRCA = "data/BRCA_TCGA_clinical.csv"
        )
        
        brca_pam50_class = fread("data/BRCA_TCGA_pam50.csv")
        
        brca_clin = fread("data/BRCA_TCGA_clinical.csv")
        # brca_clin[, .(vital_status, days_to_last_follow_up, days_to_death)][vital_status == "alive"][order(days_to_death)]
        # brca_clin[, .(vital_status, days_to_last_follow_up, days_to_death)][vital_status != "alive"][order(days_to_last_follow_up)]
        brca_clin = brca_clin[, .(submitter_id, gender, year_of_birth, race, ethnicity, primary_diagnosis, tumor_stage, age_at_diagnosis, days_to_birth, days_to_death, days_to_last_follow_up, vital_status)]
        # brca_clin = merge(brca_clin, brca_pam50_class, by = "submitter_id")
        brca_clin
        
        sample_codes = function(sample_ids){
            sdt = data.table(sample_id = sample_ids)
            
            # if(grepl("\\.", sdt$sample_id))
            sdt$cleaner = sub(".+\\.", "", sdt$sample_id)
            sdt[, c("p1", "p2", "p3", "p4", "p5") := tstrsplit(cleaner, "-", keep = 1:5)]
            sdt[, patient_id := paste(p1, p2, p3, sep = "-")]
            sdt[, sample_code := sub("[A-Z]", "", p4)]
            
            
            # pat_ids = 
            # sdt[, rm := sub("(?:[^-]*-){2}([^-]*)", "", sample_id, perl = TRUE)]
            # sdt[, patient_id := sub(rm, "", sample_id), .(sample_id)]
            # sdt[, sample_code := sub("[A-Z]", "", sub("-.+", "", sub("-", "", rm)))]
            # sdt = merge(sdt, code_dt, by = "sample_code")
            # sdt$rm = NULL
            merge(sdt[, .(sample_id, patient_id, sample_code)], code_dt, by = "sample_code")
        }
        
        sample_dts = lapply(expression_loaded, function(x){
            ids = colnames(x)[-1]
            sample_codes(ids)
        })
        
        setnames(brca_clin, "submitter_id", "patient_id")
        
        # setnames(clin_dt, c("Vital.Status", "Overall.Survival.Time.in.Days"), c("vital_status", "days_to_last_follow_up"))
        # clin_dt[, days_to_death := ifelse(vital_status == "Dead", days_to_last_follow_up, NA) ]
        InstallDataset(dataset_dirname = "BRCA_tiny", 
                        expression_data = expression_loaded$BRCA_tiny, 
                        clinical_data = brca_clin, 
                        sample_data = sample_dts$BRCA_tiny)
        InstallDataset(dataset_dirname = "BRCA", 
                        expression_data = expression_loaded$BRCA, 
                        clinical_data = brca_clin, 
                        sample_data = sample_dts$BRCA)
        
        target_exp_dt = fread("../SF_target_RNAseq/TARGET_expression_wide.csv")
        target_exp_dt = target_exp_dt[, -"RNA-Seq.PBDC-PB.TARGET-15-SJMPAL011911-03A.1-01D"]
        target_exp_dt[1:5, 1:5]
        
        target_surv_dt = fread("../SF_target_RNAseq/clinical_merged.survival.csv")
        target_clin_dt = fread("../SF_target_RNAseq/clinical_merged.csv")    
        target_clin_dt$i = NULL
        cn = colnames(target_clin_dt)
        cn = cn[!grepl("V[0-9]+$", cn)]
        target_clin_dt = target_clin_dt[, cn, with = FALSE]
        target_clin_dt = unique(target_clin_dt)    
        target_clinical_dt = merge(target_clin_dt, target_surv_dt, by = "patient_id", all = TRUE)
        
        ref_gr = rtracklayer::import.gff("~/gencode.v35.annotation.gtf.gz", feature.type = "gene")
        names(ref_gr) = ref_gr$gene_id
        target_exp_dt[, gene_id := ref_gr[gene_id]$gene_name]
        setnames(target_exp_dt, "gene_id", "gene_name")
        sum(duplicated(target_exp_dt$gene_name))
        
        target_exp_dt[duplicated(gene_name),]$gene_name
        tmp = melt(target_exp_dt, id.vars = "gene_name")
        tmp = tmp[, .(value = round(mean(value))), .(gene_name, variable)]
        dim(tmp)
        tmp
        tmp2 = dcast(tmp, gene_name~variable, value.var = "value")
        dim(tmp2)
        dim(target_exp_dt)
        
        target_sample_dt = sample_codes(colnames(target_exp_dt)[-1])
        target_sample_dt = target_sample_dt[patient_id %in% target_clinical_dt$patient_id]
        
        tmp3 = tmp2[, colnames(tmp2) %in% c("gene_name", target_sample_dt$sample_id), with = FALSE]
        
        InstallDataset(dataset_dirname = "TARGET", 
                        expression_data = tmp3, 
                        clinical_data = target_clinical_dt, 
                        sample_data = target_sample_dt)
    }
    
}

