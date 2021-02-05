install_dataset = function(dataset_dirname, 
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
    
    #setup dir
    dest = file.path(install_location, dataset_dirname)
    dir.create(dest, showWarnings = FALSE, recursive = TRUE)
    
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
}

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
    brca_clin = merge(brca_clin, brca_pam50_class, by = "submitter_id")
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
    install_dataset(dataset_dirname = "BRCA_tiny", 
                    expression_data = expression_loaded$BRCA_tiny, 
                    clinical_data = brca_clin, 
                    sample_data = sample_dts$BRCA_tiny)
    install_dataset(dataset_dirname = "BRCA", 
                    expression_data = expression_loaded$BRCA, 
                    clinical_data = brca_clin, 
                    sample_data = sample_dts$BRCA)
    
    target_exp_dt = fread("../SF_target_RNAseq/TARGET_expression_wide.csv")
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
    
    install_dataset(dataset_dirname = "TARGET", 
                    expression_data = tmp2, 
                    clinical_data = target_clinical_dt, 
                    sample_data = sample_codes(colnames(target_exp_dt)[-1]))
    }
    
}

