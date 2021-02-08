# clinical_files = file.path(dir("installed_datasets/", full.names = TRUE), "clinical.csv")
# names(clinical_files) = basename(dirname(clinical_files))
# stopifnot(file.exists(clinical_files))
# clinical_files = as.list(clinical_files)
# clinical_loaded = lapply(clinical_files, function(x)NULL)
# 
# sample_files = file.path(dir("installed_datasets/", full.names = TRUE), "samples.csv")
# names(sample_files) = basename(dirname(sample_files))
# stopifnot(file.exists(sample_files))
# sample_files = as.list(sample_files)
# sample_loaded = lapply(sample_files, function(x)NULL)
# 
# #TCGA seems to use gencode.v22
# #expression files should have a first column of 
# #gene_names named gene_name and expression after that
# 
# expression_files = file.path(dir("installed_datasets/", full.names = TRUE), "expression.csv")
# names(expression_files) = basename(dirname(expression_files))
# stopifnot(file.exists(expression_files))
# expression_files = as.list(expression_files)
# expression_loaded = lapply(expression_files, function(x)NULL)

installed_datasets = dir("installed_datasets/", full.names = TRUE)
names(installed_datasets) = basename(installed_datasets)

app_datasets = lapply(installed_datasets, CreateAppDataset.installed)

dataset_names = names(app_datasets)
# lapply(list(expression_files, expression_loaded, clinical_files, clinical_loaded, sample_files, sample_loaded), function(x){
#     eq = setequal(names(x), dataset_names)
#     stopifnot(eq)
#     eq
# })

#want tiny first so it's loaded initially
dataset_names = dataset_names[order(grepl("tiny", dataset_names), decreasing = TRUE)]

code_dt = fread("sample_codes.txt", colClasses = rep("character", 3), col.names = c("code", "long_name", "short_name"))
