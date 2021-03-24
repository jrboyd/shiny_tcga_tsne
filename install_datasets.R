source("AppDataset.R")

#this installs TCGA and TARGET datasets
if(FALSE){
    {
        library(data.table)
        code_dt = fread("sample_codes.txt", colClasses = rep("character", 3), col.names = c("sample_code", "sample_type", "sample_type_short"))
        
        expression_files = list(
            BRCA_tiny = "data/BRCA_TCGA_expression.tiny2.csv",
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
        
        pam_dt = fread("data/tcga_brca_metadata_pam50.csv")
        pam_dt = pam_dt[, .(sample_id = bcr_patient_barcode, PAM50_Basal, PAM50_Her2, PAM50_LumA, PAM50_LumB, PAM50_Normal, pam_call, pam_max)]
        
        sample_dts = lapply(sample_dts, function(x){
            merge(x, pam_dt, by = "sample_id")
        })
        
        brca_clin$gender  = factor(brca_clin$gender, levels = c("female", "male", "--"))
        
        stage_lev = c("stage i", "stage ia", "stage ib", "stage ii", 
                      "stage iia", "stage iib", "stage iii", 
                      "stage iiia", "stage iiib", "stage iiic", 
                      "stage iv", "stage x", "not reported", "--")
        
        brca_clin$tumor_stage = factor(brca_clin$tumor_stage, levels = stage_lev)
        
        sample_dts$BRCA$sample_type = factor(sample_dts$BRCA$sample_type, levels = c("Solid Tissue Normal", "Primary Solid Tumor", "Metastatic"))
        sample_dts$BRCA_tiny$sample_type = factor(sample_dts$BRCA_tiny$sample_type, levels = c("Solid Tissue Normal", "Primary Solid Tumor", "Metastatic"))
        
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
        
        snp_dt = fread("ik_snps_annotated.csv")
        del_dt = snp_dt[group %in% c("SIFT", "polyphen")]
        traits = del_dt[, .N, .(Trait, group)]$Trait
        traits2simp = c("Deleterious", "Mixed", "Mixed", "Tolerated", "Tolerated", "Deleterious", "Mixed")
        names(traits2simp) = traits
        del_dt[, simple_trait := traits2simp[Trait]]
        table(del_dt$simple_trait)
        any(duplicated(del_dt$sample))
        del_dt$simple_trait = factor(del_dt$simple_trait, levels = c("Deleterious", "Mixed", "Tolerated"))
        del_dt[, sample_id := sub("\\.SRR.+", "", sample)]
        simple_del_dt =del_dt[order(simple_trait)][!duplicated(sample_id)][, .(sample_id, simple_snp = simple_trait)]
        
        target_sample_dt =merge(target_sample_dt, simple_del_dt, by = "sample_id", all.x = TRUE)
        target_sample_dt[is.na(simple_snp), simple_snp := "No AA Impact"]
        table(target_sample_dt$simple_snp)
        
        InstallDataset(dataset_dirname = "TARGET", 
                       expression_data = tmp3, 
                       clinical_data = target_clinical_dt, 
                       sample_data = target_sample_dt)
    }
    
}