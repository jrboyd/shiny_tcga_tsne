library(data.table)
library(magrittr)
source("functions.R")
exp_mat = load_expression("installed_datasets/TARGET_IKZF1/expression.csv")
total_reads = apply(exp_mat, 2, function(x)sum(x))
# assay_dt = data.table(assay_id = colnames(exp_mat))
# assay_dt[, sample_id := sapply(strsplit(assay_id, "-"), function(x)paste(x[1:5], collapse = "-"))]

#splice 1 and 2, deleterious SNPs, or low expression
sample_dt = fread("installed_datasets/TARGET_IKZF1/samples.csv")

low_exp = sample_dt[IKZF1_bin == "(0,171]"]$patient_id %>% unique
bad_splice = sample_dt[splicing %in% c("cluster 1", "cluster 2")]$patient_id %>% unique
del_SNP = sample_dt[SNP_deleterious == TRUE]$patient_id %>% unique

seqsetvis::ssvFeatureVenn(list(expression = low_exp, splicing = bad_splice, SNP = del_SNP), circle_colors = seqsetvis::safeBrew(6, "reds")[c(2, 4, 6)])

ik_malign = unique(c(low_exp, bad_splice, del_SNP))
# ik_malign = unique(c(bad_splice, del_SNP))
length(ik_malign)

ik_benefit = sample_dt[intron1_SNPs == 6]$patient_id %>% unique
length(ik_benefit)

seqsetvis::ssvFeatureVenn(list(malign = ik_malign, beneficial = ik_benefit), circle_colors = c("red", "lightblue"))

surv_dt = fread("installed_datasets/TARGET_IKZF1/clinical.csv")

surv_dt$Cell.of.Origin %>% table

res_impacted = run_survival_by_patients(surv_dt, ik_malign)
res_benefit = run_survival_by_patients(surv_dt, ik_benefit)

res_impacted.preb = run_survival_by_patients(surv_dt[Cell.of.Origin == "B-Precursor"], ik_malign)
res_impacted.preb$plot = res_impacted.preb$plot + labs(title = "B-Precursor - malign")
res_benefit.preb = run_survival_by_patients(surv_dt[Cell.of.Origin == "B-Precursor"], ik_benefit)
res_benefit.preb$plot = res_benefit.preb$plot + labs(title = "B-Precursor - benefit")

res_impacted.tcell = run_survival_by_patients(surv_dt[Cell.of.Origin == "T Cell ALL"], ik_malign)
res_impacted.tcell$plot = res_impacted.tcell$plot + labs(title = "T Cell ALL - malign")
res_benefit.tcell = run_survival_by_patients(surv_dt[Cell.of.Origin == "T Cell ALL"], ik_benefit)
res_benefit.tcell$plot = res_benefit.tcell$plot + labs(title = "T Cell ALL - benefit")

pg_surv_origin = cowplot::plot_grid(res_impacted$plot + labs(title = "Any origin - malign"), 
                                    res_benefit$plot + labs(title = "Any origin - benefit"), 
                                    res_impacted.preb$plot, res_benefit.preb$plot, 
                                    res_impacted.tcell$plot, res_benefit.tcell$plot, ncol = 2)

ggsave("survival_origin.png", pg_surv_origin, width = 8, height = 12)

res_phase3 = run_survival_by_patients(surv_dt, surv_dt[file == "target_data/metadata/TARGET_ALL_ClinicalData_Phase_III_20181213.xlsx"]$patient_id, legend_title = "PhaseIII")
res_phase3$plot

res_cell_t = run_survival_by_patients(surv_dt, surv_dt[`Cell.of.Origin` == "T Cell ALL"]$patient_id, legend_title = "T Cell")


res_cell_b = run_survival_by_patients(surv_dt, surv_dt[`Cell.of.Origin` == "B-Precursor"]$patient_id, legend_title = "B-Precursor")
cowplot::plot_grid(res_cell_t$plot, res_cell_b$plot)

"B-Precursor"


surv_dt$Cell.of.Origin



res_impacted$plot
res_benefit$plot

ik_malign.strict = setdiff(ik_malign, ik_benefit)
ik_benefit.strict = setdiff(ik_benefit, ik_malign)
ik_ambiguous = intersect(ik_benefit, ik_malign)

res_impacted2 = run_survival_by_patients(surv_dt, ik_malign.strict)
res_benefit2 = run_survival_by_patients(surv_dt, ik_benefit.strict)
res_ambiguous = run_survival_by_patients(surv_dt, ik_ambiguous)

pg_status_surv = cowplot::plot_grid(
    res_impacted2$plot + labs(title = "Malign"),
    res_benefit2$plot + labs(title = "Benefit"),
    res_ambiguous$plot + labs(title = "Ambiguous"), nrow = 1)
ggsave("survival_by_status.png", pg_status_surv, width = 12, height = 4 )

ik_malign.strict.samples = sample_dt[patient_id %in% ik_malign.strict & sample_code == "9"]$sample_id
ik_benefit.strict.samples = sample_dt[patient_id %in% ik_benefit.strict & sample_code == "9"]$sample_id
ik_ambiguous.samples = sample_dt[patient_id %in% ik_ambiguous & sample_code == "9"]$sample_id

intersect(colnames(exp_mat), ik_malign.strict.samples)
intersect(colnames(exp_mat), ik_benefit.strict.samples)

stopifnot(ik_malign.strict.samples %in% colnames(exp_mat))
stopifnot(ik_benefit.strict.samples %in% colnames(exp_mat))

other = setdiff(colnames(exp_mat), c(ik_malign.strict.samples, ik_benefit.strict.samples))

v1 = exp_mat["MT-CYB", ik_malign.strict.samples]
v2 = exp_mat["MT-CYB", ik_benefit.strict.samples]
tmp = data.frame(value = c(v1, v2), group = c(rep("malign", length(v1)), rep("benefit", length(v2))))
boxplot(value~group, tmp)

v1 = exp_mat["MT-CYB", other]
v2 = exp_mat["MT-CYB", ik_benefit.strict.samples]
tmp = data.frame(value = c(v1, v2), group = c(rep("other", length(v1)), rep("benefit", length(v2))))
boxplot(value~group, tmp)

res_de_imp_vs_ben = run_DE(exp_mat = exp_mat, a = ik_malign.strict.samples, b = ik_benefit.strict.samples, on_windows = TRUE, min_log2FoldChange = .5, max_padj = .05)
res_de_imp_vs_other = run_DE(exp_mat = exp_mat, a = ik_malign.strict.samples, b = other, on_windows = TRUE, min_log2FoldChange = .5, max_padj = .05)
res_de_ben_vs_other = run_DE(exp_mat = exp_mat, a = ik_benefit.strict.samples, b = other, on_windows = TRUE, min_log2FoldChange = .5, max_padj = .05)
res_de_amb_vs_other = run_DE(exp_mat = exp_mat, a = ik_ambiguous.samples, b = setdiff(other, ik_ambiguous.samples), on_windows = TRUE, min_log2FoldChange = .5, max_padj = .05)

all_res = list(amb_vs_other = res_de_amb_vs_other, 
               imp_vs_other = res_de_imp_vs_other, 
               ben_vs_other = res_de_ben_vs_other,
               imp_vs_ben = res_de_imp_vs_ben
)

all_plots = lapply(names(all_res), function(nam){
    df = all_res[[nam]]
    df$log2BaseMean = log2(df$baseMean+1)
    df$is_sig = df$log2BaseMean > 3 & df$padj < 10e-20
    df$padj_neglog = -log10(df$padj)
    df$gene_name = rownames(df)
    ggplot(df, aes(x = log2FoldChange, y = log2BaseMean, color = padj_neglog)) + 
        geom_point() +
        ggrepel::geom_label_repel(data = df[order(df$padj),][1:15,], aes(label = gene_name), color = "black", max.overlaps = 20) +
        labs(title = nam)
})

pg_de_plots = cowplot::plot_grid(plotlist = all_plots)
ggsave("DE_plots.png", pg_de_plots, width = 12, height = 10)

all_res.sig = lapply(all_res, function(df){
    df$log2BaseMean = log2(df$baseMean+1)
    df$is_sig = df$log2BaseMean > 3 & df$padj < 10e-20
    subset(df, is_sig == TRUE)
})

all_res.top200 = lapply(all_res, function(df){
    rownames(df[order(df$padj),][1:200,])
})
seqsetvis::ssvFeatureUpset(all_res.top200
)

lapply(all_res.sig, nrow)
gl = unique(unlist(all_res.top200))

run_tsne

df = res_de_imp_vs_ben
df$log2BaseMean = log2(df$baseMean+1)
df$is_sig = df$log2BaseMean > 3 & df$padj < 10e-20
ggplot(df, aes(x = log2FoldChange, y = log2BaseMean, color = is_sig)) + 
    geom_point()

res_de_imp_vs_other



tsne_res = run_tsne(log10(exp_mat+1), gl)
tsne_res$group = "other"
tsne_res[tsne_res$sample_id %in% ik_malign.strict.samples]$group = "malign"
tsne_res[tsne_res$sample_id %in% ik_benefit.strict.samples]$group = "benefit"
tsne_res[tsne_res$sample_id %in% ik_ambiguous.samples]$group = "ambiguous"

clin_dt = surv_dt[, .(patient_id, cell = Cell.of.Origin, phase = file)]
clin_dt[grepl("Phase_I_", phase), phase := "Phase_I"]
clin_dt[grepl("Phase_II_", phase), phase := "Phase_II"]
clin_dt[grepl("Phase_III_", phase), phase := "Phase_III"]
clin_dt[grepl("Dicentric", phase), phase := "Dicentric"]

tsne_res[, patient_id := sub(".+\\.", "", sample_id)]
tsne_res[, patient_id := sapply(strsplit(patient_id, "-"), function(x)paste(x[1:3], collapse = "-"))]


tsne_res.clin = merge(tsne_res, clin_dt, by = "patient_id")
tsne_res.clin$total_reads = total_reads[tsne_res$sample_id]

p_cell = ggplot(tsne_res.clin, aes(x = x, y = y, color = cell)) + 
    annotate("point", x= tsne_res$x, y = tsne_res$y, color = "gray85", size = .8) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    facet_wrap(~cell)+ theme(panel.background = element_blank()) +
    theme(legend.position = "bottom")

p_phase = ggplot(tsne_res.clin, aes(x = x, y = y, color = phase)) + 
    annotate("point", x= tsne_res$x, y = tsne_res$y, color = "gray85", size = .8) +
    geom_point() +
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~phase)+ theme(panel.background = element_blank())+
    theme(legend.position = "bottom")

ggplot(tsne_res.clin, aes(x = phase, y = total_reads)) +
    geom_boxplot()

ggplot(tsne_res.clin, aes(x = cell, y = total_reads)) +
    geom_boxplot()

ggplot(tsne_res.clin, aes(x = paste(phase, cell), y = total_reads)) +
    # geom_boxplot() +
    ggbeeswarm::geom_beeswarm() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


pg_cell_phase = cowplot::plot_grid(p_cell, p_phase, rel_widths = c(3, 2))


p_tsne_groups = ggplot(tsne_res, aes(x = x, y = y, color = group)) + 
    annotate("point", x= tsne_res$x, y = tsne_res$y, color = "gray85", size = .8) +
    geom_point() +
    scale_color_manual(values = c("ambiguous" = "purple", benefit = "blue", malign = "red", other = "black")) +
    facet_wrap(~group)+ theme(panel.background = element_blank())
#intron 1 DE

goi_a = c("CD79B", "CTNND1", "BLNK", "IKZF1", "SUZ12", "JAG1")
goi_b = c("RAG1", "CD34", "KLF2", "JAG1", "IL6R", "S1PR1", "SPIB")
goi = union(goi_a, goi_b)
stopifnot(goi %in% rownames(exp_mat))



norm_mat = apply(exp_mat, 2, function(x)log2( (x / sum(x)*1e6) + 1))
range(norm_mat)

exp_dt.full = as.data.table(reshape2::melt(exp_mat))
setnames(exp_dt.full, c("gene_name", "sample_id", "count"))
exp_dt.full = merge(exp_dt.full, tsne_res, by = "sample_id")
agg_dt.full = exp_dt.full[, .(count = mean(count)), .(group, gene_name)]
ggplot(agg_dt.full, aes(x = group, y = count)) +
    geom_boxplot() +
    scale_y_log10()

agg_dt.full[, .(Million_Counts = sum(count)/1e6), .(group)]
qtile = agg_dt.full[, quantile(count, 0:20/20), .(group)]
qtile$q = rep(0:20/20, 4)
ggplot(qtile[q < 1], aes(x = q, y = V1, color = group)) +
    geom_path() +
    scale_y_log10()


# exp_dt = as.data.table(reshape2::melt(exp_mat[goi,]))
exp_dt = as.data.table(reshape2::melt(norm_mat[goi,]))
setnames(exp_dt, c("gene_name", "sample_id", "count"))
exp_dt$total_reads = total_reads[exp_dt$sample_id]
# 
# 
# boxplot(exp_mat)


exp_dt = merge(exp_dt, tsne_res, by = "sample_id")

exp_dt$group = factor(exp_dt$group, levels = c("other", "benefit", "ambiguous", "malign"))

p_status_bees = ggplot(exp_dt, aes(x = group, y = count)) + 
    ggbeeswarm::geom_beeswarm() +
    facet_wrap(~gene_name, scales = "free_y")

p_status_violin = ggplot(exp_dt, aes(x = group, y = count, fill = group)) + 
    geom_violin() +
    facet_wrap(~gene_name, scales = "free_y") +
    scale_fill_manual(values = c("ambiguous" = "plum", benefit = "lightblue", malign = "red", other = "gray")) +
    theme(panel.background = element_blank())

box_dt = exp_dt[, boxplot.stats(count)$stats, .(gene_name, group)]
box_dt$stat = rep(c("ymin", "lower", "middle", "upper", "ymax"), nrow(box_dt)/5)
box_dt = dcast(box_dt, gene_name+group~stat, value.var = "V1")

p_status_box = ggplot(box_dt, aes(x = group, ymin = ymin, ymax = ymax, lower = lower, middle = middle, upper = upper, fill = group)) + 
    geom_boxplot(stat = "identity") +
    facet_wrap(~gene_name, scales = "free_y") +
    scale_fill_manual(values = c("ambiguous" = "plum", benefit = "lightblue", malign = "red", other = "gray")) +
    theme(panel.background = element_blank())
p_status_box

ggplot(exp_dt[gene_name == "IKZF1"], aes(x = group, y = total_reads)) + geom_boxplot() +
    # facet_wrap(scales = "free_y") +
    labs(title = 'total reads')

lapply(all_res.top200, function(x){intersect(x, goi)})
lapply(all_res, function(x)x[goi,])


exp_dt[, norm_count := scales::rescale(count), .(gene_name)]
p_tsne_expression = ggplot(exp_dt, aes(x = x, y = y, color = norm_count)) + 
    geom_point(size = .6) +
    facet_wrap(~gene_name, scales = "free_y") + 
    scale_color_viridis_c() +
    theme(panel.background = element_blank()) +
    labs(color = "relative\nexpression")

ggplot(exp_dt[gene_name == "IKZF1"], aes(x = x, y = y, color = total_reads)) + 
    geom_point() +
    # facet_wrap(~gene_name, scales = "free_y") + 
    scale_color_viridis_c() +
    theme(panel.background = element_blank()) +
    labs(color = "total reads", title = "Library size")

pg_tsne = cowplot::plot_grid(p_tsne_groups, p_tsne_expression, rel_widths = c(1, 1.3))
ggsave("ikaros_status_goi_tsne.png", pg_tsne, width = 15.4, height = 6)
ggsave("ikaros_status_cell_and_phase.png", pg_cell_phase, width = 15.4, height = 6)

p_status_violin
p_status_box = p_status_box + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = 'RNA-seq CPM', title = "Expression of Ikaros relevant genes by Ikaros status")
ggsave("ikaros_status_goi_boxplots.png", p_status_box, width = 11, height = 8)

tsne_res = merge(tsne_res, sample_dt, by = "sample_id")
ggplot(tsne_res, aes(x = x, y = y, color = splicing)) +
    geom_point() +
    facet_wrap(~splicing)

setnames(tsne_res, "patient_id.x", "patient_id")
tsne_res.kary = merge(tsne_res, surv_dt[, .(patient_id, Karyotype)], by = "patient_id")

tsne_res.kary$Karyotype

tsne_res.kary[grepl("dic\\(9", Karyotype),]

ggplot(tsne_res.kary, aes(x = x, y = y, color = grepl("dic\\(9", Karyotype))) +
    geom_point() +
    facet_wrap(~splicing)
