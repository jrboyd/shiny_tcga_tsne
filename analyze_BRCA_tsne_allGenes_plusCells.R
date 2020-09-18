library(jsonlite)
library(XML)
require(data.table)
require(magrittr)
require(pbapply)
library(ggplot2)
setwd("~/R/GDC-TCGA_data/")
source("functions_gdc_dl_and_verify.R")

cell_dt = merge(
  merge(
    fread("~/R/RNAseq_AT1_sorted/DE_res/DESeq2_normalized_counts.csv"),
    fread("~/R/RNAseq_pfizer/DE_res/DESeq2_normalized_counts.csv")),
  fread("~/R/RNAseq_BRCA/DE_res/DESeq2_normalized_counts.csv"))
cell_dt[, gene_id := sub("\\..+", "", gene_id)]
cell_mat = as.matrix(cell_dt[, -1])
rownames(cell_mat) = cell_dt$gene_id

# source("functions_gdc_dl_and_verify.R")
# gdc_dl("jsons/BRCA_metadata.cart.2018-12-13.json")
# gdc_dl("jsons/BRCA_biospecimen.cart.2018-12-13.json")
c_file = "analyze_BRCA_tsne_allGenes_plusCells.save"
if(file.exists(c_file)){
  load(c_file)
}else{
  library(genefu)
  data("pam50.robust")
  my_model = pam50.robust
  pam50 = my_model$centroids
  subGN = function(match, replace, df = pam50){
    rownames(df)[rownames(df) == match] = replace
    df
  }
  pam50 = subGN("CDCA1", "NUF2")
  pam50 = subGN("KNTC2", "NDC80")
  pam50 = subGN("ORC6L", "ORC6")
  
  ref_gr = rtracklayer::import.gff("~/gencode.v22.annotation.gtf.gz", format = "gtf", feature.type = "gene")
  names(ref_gr) = sub("\\..+", "", ref_gr$gene_id)
  
  
  json2info = function(json_file){
    json_content = fromJSON(json_file)
    info = list()
    info$case_ids = sapply(json_content$cases, function(x)x$case_id)
    info$file_names = json_content$file_name
    info$file_ids = json_content$file_id
    info$entity_id = sapply(json_content$associated_entities, function(x)x$entity_id)
    info$content = json_content
    return(info)
  }
  
  df = read.table("sample_codes", sep = "\t")
  codes = sub(" ", "0", format(df$V1, width = 2))
  sample_code2type = as.character(df$V2)
  names(sample_code2type)  = codes
  
  json2df = function(json_file){
    json_content = fromJSON(json_file)
    df = data.frame(case_ids = sapply(json_content$cases, function(x)x$case_id))
    df$file_names = json_content$file_name
    df$file_ids = json_content$file_id
    df$entity_id = sapply(json_content$associated_entities, function(x)x$entity_id)
    df$entity_submitter_id = sapply(json_content$associated_entities, function(x)x$entity_submitter_id)
    df$sample_type = sample_code2type[strsplit(df$entity_submitter_id, "-") %>% 
                                        sapply(., function(x)x[4]) %>% sub("[A-Z]$", "", .)]
    df$sample_source = strsplit(df$entity_submitter_id, "-") %>% 
      sapply(., function(x)paste(x[1:3], collapse = "-"))
    # info$content = json_content
    
    return(df)
  }
  # clin_info = json2info("jsons/prostate_clinical_metadata.cart.2017-10-20T18-50-23.087644.json")
  # mir_info = json2info(json_file = "jsons/prostate_miR_metadata.cart.2017-10-20T18-45-58.687440.json")
  # bio_info = json2info("jsons/prostate_biospecimen_metadata.cart.2017-10-23T17-15-47.468588.json")
  json_file = "jsons/BRCA_metadata.cart.2018-12-13.json"
  json_content = fromJSON(json_file)
  df = data.frame(case_ids = sapply(json_content$associated_entities, function(x)x$case_id))
  df$file_names = json_content$file_name
  df$file_ids = json_content$file_id
  df$entity_id = sapply(json_content$associated_entities, function(x)x$entity_id)
  df$entity_submitter_id = sapply(json_content$associated_entities, function(x)x$entity_submitter_id)
  df$sample_type = sample_code2type[strsplit(df$entity_submitter_id, "-") %>% 
                                      sapply(., function(x)x[4]) %>% sub("[A-Z]$", "", .)]
  df$sample_source = strsplit(df$entity_submitter_id, "-") %>% 
    sapply(., function(x)paste(x[1:3], collapse = "-"))
  
  
  master_df = df#json2df(json_file = "jsons/BRCA_metadata.cart.2018-12-13.json")
  htseq_df =  master_df[grepl("FPKM", master_df$file_names),]
  df_codes = matrix(unlist(strsplit(htseq_df$entity_submitter_id, '-')), ncol = 7, byrow = T)
  patient_ids = apply(df_codes[,1:3], 1, function(x)paste(x, collapse = "-"))
  # sub("[A-Z]", "", df_codes[,4])
  htseq_df$sample_type = gsub(" ", "_", htseq_df$sample_type)
  
  htseq_files = paste0("download/", htseq_df$file_ids, "/", htseq_df$file_names)
  tcga_file = "tcga_expression.save"
  if(file.exists(tcga_file)){
    load(tcga_file)
  }else{
    all_df = pbapply::pblapply(htseq_files, fread)
    names(all_df) = htseq_df$entity_submitter_id# paste(htseq_df$sample_source, htseq_df$sample_type)
    all_df = rbindlist(all_df, use.names = TRUE, idcol = "patient_id")
    colnames(all_df)[2:3] = c("gene_id", "count")
    all_df = all_df[!grepl("^__", gene_id)]
    
    tmp = all_df[, .N, by = .(patient_id, gene_id)]
    table(tmp$N)
    
    all_df_wide = dcast(all_df, "gene_id~patient_id", value.var = "count")
    save(all_df, all_df_wide, file = tcga_file)
  }
  
  mat = as.matrix(all_df_wide[,-1])
  rownames(mat) = all_df_wide$gene_id
  rownames(mat) = sub("\\..+", "", rownames(mat))
  mat[1:5, 1:5]
  
  common = intersect(rownames(mat), names(ref_gr))
  common = intersect(common, cell_dt$gene_id)
  mat = mat[common,]
  mat = cbind(mat, cell_mat[common,])
  train_k = ref_gr[rownames(mat)]$gene_name %in% rownames(pam50)
  sum(train_k)
  
  
  library(Rtsne)
  options(mc.cores = 12)
  dim(mat)
  mat = normalizeQuantiles(mat)
  # train_mat = mat[train_k,]
  train_mat= mat
  
  sum(duplicated(ref_gr[rownames(mat)]$gene_name))
  # rsum = rowSums(mat)
  # table(rsum > 5000)
  
  tsne_gene = Rtsne::Rtsne(train_mat, num_threads = 20, check_duplicates = FALSE, perplexity = 10)
  plot(tsne_gene$Y, pch = 16, col = rgb(0,0,0,1))
  
  tsne_patient = Rtsne::Rtsne(t(train_mat), num_threads = 20)
  
  tsne_df = as.data.table(tsne_patient$Y)
  colnames(tsne_df) = c("x", "y")
  tsne_df$bcr_patient_barcode = colnames(mat)
  
  json_file = "jsons/BRCA_clinical_metadata.json"
  json_content = fromJSON(json_file)
  df = data.frame(case_ids = sapply(json_content$associated_entities, function(x)x$case_id))
  df$file_names = json_content$file_name
  df$file_ids = json_content$file_id
  df$entity_id = sapply(json_content$associated_entities, function(x)x$entity_id)
  df$entity_submitter_id = sapply(json_content$associated_entities, function(x)x$entity_submitter_id)
  df$sample_type = sample_code2type[strsplit(df$entity_submitter_id, "-") %>% 
                                      sapply(., function(x)x[4]) %>% sub("[A-Z]$", "", .)]
  df$sample_source = strsplit(df$entity_submitter_id, "-") %>% 
    sapply(., function(x)paste(x[1:3], collapse = "-"))
  
  clinical_dat = fread("BRCA_clinical.tsv")
  clinical_dat = clinical_dat[,which(!duplicated(colnames(clinical_dat))), with = FALSE]

    tsne_df[, submitter_id := paste(tstrsplit(bcr_patient_barcode, "-", keep = 1)[[1]], 
                                  tstrsplit(bcr_patient_barcode, "-", keep = 2)[[1]],
                                  tstrsplit(bcr_patient_barcode, "-", keep = 3)[[1]], sep = "-")]
  tsne_df[grepl(" ", bcr_patient_barcode), submitter_id := sub(" .+", "", submitter_id)]
  tsne_df = merge(tsne_df, clinical_dat, all.x = TRUE)
  
  save(tsne_df, tsne_gene, tsne_patient, ref_gr, clinical_dat, mat, pam50, file = c_file)
}

# tsne_df = tsne_df[-5,]
tsne_df = tsne_df[!bcr_patient_barcode == "MCF10A 2"]
mat = mat[, tsne_df$bcr_patient_barcode]
tsne_df[, sample_type_code := tstrsplit(bcr_patient_barcode, "-", keep = 4)]
tsne_df[, sample_type_code := sub("[A-Z]", "", sample_type_code)]

code2type = c("01" = "tumor", "11" = "normal", "06" = "metastasis")
tsne_df[, sample_type := code2type[sample_type_code]]
tsne_df[is.na(sample_type_code), sample_type := "cell line"]
table(tsne_df$sample_type)
tsne_df[, sample_source := 'TCGA']
tsne_df[is.na(sample_type_code), sample_source := "cell line"]

ggplot(tsne_df[sample_source == "TCGA"], aes(x = x, y = y, color = sample_type)) + geom_point()
ggplot(tsne_df[sample_source != "TCGA"], aes(x = x, y = y, color = submitter_id )) + 
  annotate("point", x = tsne_df[sample_source == "TCGA"]$x, y = tsne_df[sample_source == "TCGA"]$y, color = 'gray') +
  geom_point() + facet_wrap("submitter_id")
###start
make_plot = function(gn, cap = 3, split = 1){
  goi = subset(ref_gr, gene_name == gn) %>% names
  
  tsne_df$goi_val = mat[goi,]
  # plot(tsne_patient$Y, pch = 16, col = rgb(0,0,0,1))
  
  to_num = function(dt, colname){
    dt = copy(dt)
    suppressWarnings({
      dt = dt[!is.na(as.numeric(get(colname)))]
    })
    dt[[colname]] = as.numeric(dt[[colname]])
    dt
    # dt[is.na(as.numeric(get(colname)))]
  }
  
  sub_df = tsne_df[, .N, by = .(primary_diagnosis)][N > 10]
  tsne_df[primary_diagnosis %in% sub_df$primary_diagnosis]
  
  val = log2(tsne_df$goi_val)
  val = (val - mean(val)) / sd(val)
  
  val = ifelse(val > cap, cap, val)
  val = ifelse(val < -cap, -cap, val)
  tsne_df$goi_z = val
  tsne_df$goi_grp = "mid"
  
  tsne_df$goi_grp[val < -split] = "low"
  tsne_df$goi_grp[val > split] = "high"
  tsne_df$goi_grp = factor(tsne_df$goi_grp, levels = c("low", "mid", "high"))
  # hist(val, breaks = 20)
  p = ggplot(tsne_df, 
             aes_string(x = "x", y = "y", color = "goi_z")) + 
    geom_point() + 
    scale_color_gradientn(colours = c("blue", "darkgray", "red"), 
                          limits = c(-cap, cap)) +
    facet_wrap("goi_grp") +
    theme_classic() + 
    labs(title = gn, #subtitle = "patient clustering by gene expression",
         x = "tsne-1", y = "tsne-2", color = "z-score") +
    theme(panel.background = element_rect(fill = "white"))
  
  p2 = ggplot(tsne_df, 
              aes_string(x = "x", y = "y", color = "goi_grp")) + 
    geom_density2d(bins = 2) + 
    scale_color_manual(values = c("low" = "blue", "mid" = "darkgray", "high" = "red")) +
    theme_classic() + 
    guides(color = "none") +
    labs(title = gn, x = "tsne-1", y = "tsne-2") +
    theme(panel.background = element_rect(fill = "white"))
  
  plots[[gn]] <<- p
  sides[[gn]] <<- p2
}

setwd("~/R/GDC-TCGA_data/")
dir.create(sub(".save", "", c_file))
setwd(sub(".save", "", c_file))
plots = list()
sides = list()
make_plot("RUNX1")
make_plot("RUNX2")
make_plot("ESR1")
make_plot("GREB1")
pg1 = cowplot::plot_grid(
  cowplot::plot_grid(plotlist = plots),
  cowplot::plot_grid(plotlist = sides), 
  rel_widths = c(2.5, 1)
)

plots = list()
sides = list()
make_plot("EZH1")
make_plot("EZH2")
make_plot("PRC1")
make_plot("SUZ12")
make_plot("RYBP")
make_plot("WDR5")
#jarad2
#rng1
#rybp
# make_plot("VIM")
# make_plot("BMP6")
pg2 = cowplot::plot_grid(
  cowplot::plot_grid(plotlist = plots),
  cowplot::plot_grid(plotlist = sides), 
  rel_widths = c(2.5, 1)
)

pdf("tsne_patients.pdf", width = 16, height = 5)
pg1
pg2
dev.off()


gn = "RUNX1"

tsne_add_gene_info = function(gn, 
                              breaks = c(-.5, .5), 
                              bin_names = c("low", "mid", "high"),
                              cap = 2){
  if(!all(breaks == sort(breaks))){
    stop("breaks must be sorted ascending")
  }
  goi = subset(ref_gr, gene_name == gn) %>% names
  val = mat[goi,]
  val = log2(val + 1)
  
  val = (val - mean(val)) / sd(val)
  
  dt = tsne_df[, c("bcr_patient_barcode", "x", "y"), with = F]
  dt$goi_raw_val = val
  
  dt$goi_val = dt$goi_raw_val
  dt$goi_val = ifelse(dt$goi_val > cap, cap, dt$goi_val)
  dt$goi_val = ifelse(dt$goi_val < -cap, -cap, dt$goi_val)
  
  dt$goi_grp = bin_names[1]
  for(i in seq_along(breaks)){
    dt[goi_val > breaks[i], goi_grp := bin_names[i + 1]]
  }
  dt$goi_grp = factor(dt$goi_grp, levels = bin_names)
  hres = hist(dt$goi_raw_val, breaks = 30)
  for(b in breaks){
    lines(rep(b, 2), c(0, max(hres$counts)), col = "red")
  }
  
  dt[]
}

tsne_bin_gene_info = function(dt, 
                              xn = 10,
                              yn = 10){
  if(is.null(dt$goi_grp)){
    stop("dt must have goi_grp set. did you run tsne_add_gene_info()?")
  }
  if(is.null(dt$goi_val)){
    stop("dt must have goi_val set. did you run tsne_add_gene_info()?")
  }
  xrng = range(dt$x)
  yrng = range(dt$y)
  # xstep = round((xrng[2] - xrng[1]) / (xn-1) + .005, 2)
  # ystep = round((yrng[2] - yrng[1]) / (yn-1) + .005, 2)
  xstep = (xrng[2] - xrng[1]) / (xn-1)
  ystep = (yrng[2] - yrng[1]) / (yn-1)
  xbreaks = xstep*(seq(xn)-xn/2-.5) + mean(xrng) - .0000001
  ybreaks = ystep*(seq(yn)-yn/2-.5) + mean(yrng) - .0000001
  
  dt[, xbin := max(which(x >= xbreaks)), by = .(bcr_patient_barcode)]
  dt[, ybin := max(which(y >= ybreaks)), by = .(bcr_patient_barcode)]
  
  dtm = dt[, .(weight = .N, goi_val = mean(goi_val)), by = .(xbin, ybin, goi_grp)]
  dtm$x = xrng[1] + diff(xrng) * (dtm$xbin - .5) / xn
  dtm$y = yrng[1] + diff(yrng) * (dtm$ybin - .5) / yn
  dtm[, weight_norm := weight / sum(weight), by = .(goi_grp)]
  
  dtm[, exp_weight := sum(weight) / nrow(dt), by = .(xbin, ybin) ]
  dtm[, weight_over_exp := weight_norm / exp_weight]
  dtm[]
}

"ESR1" %>% 
  tsne_add_gene_info(., breaks = -.8 + c(-.2, .4), bin_names = c("neg", "?", "pos")) %>% 
  tsne_bin_gene_info()
"GREB1" %>% tsne_add_gene_info() %>% tsne_bin_gene_info()

plot_tsne_gene_circs = function(gn, 
                                breaks = c(-.5, .5), 
                                bin_names = c("low", "mid", "high"),
                                cap = 2){
  dtm = gn %>% 
    tsne_add_gene_info(., 
                       breaks = breaks, 
                       bin_names = bin_names, 
                       cap = cap) %>% tsne_bin_gene_info(.)
  ggplot(dtm, aes(x = x, y = y, size = weight, color = goi_val)) + 
    geom_point() + 
    scale_color_gradientn(colours = c("blue", "white", "red"), 
                          limits = c(-cap, cap)) + 
    theme_classic() + labs(x = "tsne-1", y = "tsne-2") +
    theme(panel.background = element_rect(fill = "black"), legend.position = "bottom") + 
    facet_wrap("goi_grp", nrow = 1) + labs(title = gn)
}

pdf("tsne_gene_circs.pdf")
plot_tsne_gene_circs("ESR1", cap = 3, breaks = -.8 + c(-.2, .4), bin_names = c("neg", "?", "pos"))
plot_tsne_gene_circs("GREB1", cap = 3, breaks = -.3+ c(-.1, .3), bin_names = c("neg", "?", "pos"))
plot_tsne_gene_circs("RUNX1", cap = 3, breaks = 1.5*c(-1, 1))
plot_tsne_gene_circs("RUNX2", cap = 3, breaks = 1.5*c(-1, 1))
plot_tsne_gene_circs("VIM", cap = 3, breaks = 1.5*c(-1, 1))
plot_tsne_gene_circs("EZH1", cap = 3, breaks = 1.5*c(-1, 1))
plot_tsne_gene_circs("PRC1", cap = 3, breaks = 1.5*c(-1, 1))
dev.off()

# pam50 = subGN("KRT14", "NDC80")
# KRT14
# KNTC2 NDC80

pdf("tsne_pam50.pdf", width = 6, height = 3.7)
for(gn in rownames(pam50)){
  p = plot_tsne_gene_circs(gn, cap = 3, breaks = 0, bin_names = c("low", "high"))
  print(p)
}
dev.off()

mat[1:5, 1:5]
common = intersect(rownames(mat), names(ref_gr))
pam_mat = mat[common,]
rownames(pam_mat) = ref_gr[common]$gene_name

pam_mat = log2(pam_mat[rownames(pam50), ] + 1)
boxplot(cbind(pam_mat, pam50), xlab = "", xaxt='n')

pam_join = limma::normalizeQuantiles(cbind(pam_mat, pam50))
boxplot(pam_join)
pam_cor = cor(pam_join[, seq_len(ncol(pam_mat))], pam_join[, ncol(pam_mat) + seq_len(ncol(pam50))], method = "spearman")
hcol = rgb(colorRamp(c("blue", "white", "red"))(0:50/50)/255)
gplots::heatmap.2(pam_cor, col = hcol, 
                  density.info = NULL, key.title = "", margins = c(10,10), 
                  labRow = "", trace = "n")
colnames(pam_cor) = paste0("PAM50_", colnames(pam_cor))
pam_cor = as.data.table(pam_cor, keep.rownames = TRUE)
colnames(pam_cor)[1] = "bcr_patient_barcode"
tsne_df = merge(tsne_df, pam_cor, by = "bcr_patient_barcode")
library(data.table)
# tmp = data.table(rn = rownames(pam_cor), id = 1:nrow(pam_cor))
# tmp = tmp[, .(tstrsplit(rn, "-", keep = 4), id), by = .(id)]
# table(unlist(tmp$V1))

pdf("tsne_pam50_corr.pdf", width = 5.5, height = 4.1)
ggplot(tsne_df, aes(x = x, y = y, color = PAM50_Basal)) + geom_point() + 
  scale_color_gradientn(colours = c("blue", "white", "red"), limits = c(-1,1)) +
  theme_classic() + theme(panel.background = element_rect(fill = "black")) +
  labs(x = "tsne-1", y = "tsne-2")
ggplot(tsne_df, aes(x = x, y = y, color = PAM50_Her2)) + geom_point() +
  scale_color_gradientn(colours = c("blue", "white", "red"), limits = c(-1,1))+
  theme_classic() + theme(panel.background = element_rect(fill = "black")) +
  labs(x = "tsne-1", y = "tsne-2")
ggplot(tsne_df, aes(x = x, y = y, color = PAM50_LumA)) + geom_point() +
  scale_color_gradientn(colours = c("blue", "white", "red"), limits = c(-1,1))+
  theme_classic() + theme(panel.background = element_rect(fill = "black")) +
  labs(x = "tsne-1", y = "tsne-2")
ggplot(tsne_df, aes(x = x, y = y, color =  PAM50_LumB)) + geom_point() +
  scale_color_gradientn(colours = c("blue", "white", "red"), limits = c(-1,1))+
  theme_classic() + theme(panel.background = element_rect(fill = "black")) +
  labs(x = "tsne-1", y = "tsne-2")
ggplot(tsne_df, aes(x = x, y = y, color =  PAM50_Normal)) + geom_point() +
  scale_color_gradientn(colours = c("blue", "white", "red"), limits = c(-1,1))+
  theme_classic() + theme(panel.background = element_rect(fill = "black")) +
  labs(x = "tsne-1", y = "tsne-2")
dev.off()

pam_col = colnames(pam_cor)[-1]
tsne_df[, pam_call := pam_col[which.max(c(PAM50_Basal,  PAM50_Her2,  PAM50_LumA,   PAM50_LumB, PAM50_Normal))], by = .(bcr_patient_barcode)]
tsne_df[, pam_max := max(c(PAM50_Basal,  PAM50_Her2,  PAM50_LumA,   PAM50_LumB, PAM50_Normal)), by = .(bcr_patient_barcode)]
lum_cut = .1
tsne_df[PAM50_LumA > lum_cut & PAM50_LumB > lum_cut, pam_call := "PAM50_Lum?"]
tsne_df[pam_max < .4, pam_call := "?"]
# ggplot(tsne_df[sample_type == "tumor"], 
# aes(x = x, y = y, color =  pam_call)) + geom_point(alpha = .2, shape = 16) + facet_wrap("pam_call")
# ggplot(tsne_df[sample_type == "tumor"][!pam_call %in% c("?", "PAM50_Normal")], 
# aes(x = x, y = y, color =  pam_call)) +geom_point()+ geom_density2d(bins = 2)

tsne_df$cell = "tcga"
tsne_df[sample_type == "cell line", cell := submitter_id]

ann_df = tsne_df[!pam_call %in% c("?")]
ann_df$sample_type = "tumor"
ann_df = rbind(ann_df, 
               copy(ann_df)[, sample_type := "normal"],
               copy(ann_df)[, sample_type := "metastasis"],
               copy(ann_df)[, sample_type := "cell line"])

ggplot() + 
  annotate("point", 
           x = tsne_df[cell == "tcga"]$x, 
           y = tsne_df[cell == "tcga"]$y,
           color = "gray") +
  geom_density2d(data = ann_df,
                 aes(x = x, y = y, color =  pam_call),
                 size = 1.5,
                 bins = 2) +
  geom_point(data = tsne_df,#[!pam_call %in% c("?")], 
             aes(x = x, y = y, color =  pam_call)) +
  labs(x = "tsne-1", y = "tsne-2") +
  scale_color_brewer(palette = "Dark2") + 
  facet_wrap("sample_type") +
  theme_classic()
ggsave("tsne_sample_type_by_PAM.pdf", width = 7.5, height = 6)
# + 
#   annotate("density2d", x = ann_df$x, y = ann_df$y, color = ann_df$pam_call, bins = 2)


ann_df = copy(tsne_df[!pam_call %in% c("?")])
ann_df$cell = "DCIS"
ann_df = rbind(ann_df, 
               copy(ann_df)[, cell := "MCF10A"],
               copy(ann_df)[, cell := "MCF10AT1"],
               copy(ann_df)[, cell := "MCF10CA1"],
               copy(ann_df)[, cell := "MCF7"],
               copy(ann_df)[, cell := "MDA231"])

ggplot(tsne_df[sample_type == "cell line"], aes(x = x, y = y)) + 
  annotate("point", 
           x = tsne_df[cell == "tcga"]$x, 
           y = tsne_df[cell == "tcga"]$y,
           color = "gray") +
  geom_density2d(data = ann_df,
                 aes(x = x, y = y, color =  pam_call),
                 size = 1.5,
                 bins = 2) +
  geom_point(fill = "black", color = "white", size = 3, shape = 21) + 
  facet_wrap("cell") + 
  labs(x = "tsne-1", y = "tsne-2") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic()
ggsave("tsne_cell_lines_by_PAM.pdf", width = 10, height = 6)


tsne_df[, is_drugged := FALSE]
tsne_df[cell == "MCF7" & grepl(" .+ ", bcr_patient_barcode), is_drugged := TRUE]
tsne_df[is_drugged == TRUE, TM_drug := tstrsplit(bcr_patient_barcode, " ", keep = 2)]



ggplot(tsne_df[is_drugged == TRUE], aes(x = x, y = y)) + 
  annotate("point", 
           x = tsne_df[cell == "tcga"]$x, 
           y = tsne_df[cell == "tcga"]$y,
           color = "gray") +
  geom_density2d(data = ann_df,
                 aes(x = x, y = y, color =  pam_call),
                 size = 1.5,
                 bins = 2) +
  geom_jitter(aes(fill = TM_drug), color = "black", size = 3, shape = 21, width = 1, height = 1) +
  scale_fill_brewer(palette = "Set1")

setwd("~/R/GDC-TCGA_data/")


