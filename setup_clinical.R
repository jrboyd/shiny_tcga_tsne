clinical_files = list(
    BRCA = "data/BRCA_TCGA_clinical.csv"
)

BRCA_clinical = fread("data/BRCA_TCGA_pam50.csv")

clinical_loaded = lapply(clinical_files, function(x)NULL)
clinical_loaded$BRCA_tiny = BRCA_clinical
clinical_loaded$BRCA = BRCA_clinical

process_clinical = function(name, clinical_dt){
    if(name == "BRCA" | name == "BRCA_tiny"){
        dt = clinical_dt[, .(submitter_id, gender, year_of_birth, race, ethnicity, primary_diagnosis, tumor_stage, age_at_diagnosis, days_to_death, days_to_birth)]
    }
}