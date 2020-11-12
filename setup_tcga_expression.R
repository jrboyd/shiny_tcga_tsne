library(genefu)

#TCGA seems to use gencode.v22
#expression files should have a first column of 
#gene_names named gene_name and expression after that
expression_files = list(
    BRCA_tiny = "data/BRCA_TCGA_expression.tiny.csv",
    BRCA = "data/BRCA_TCGA_expression.csv"
)

expression_loaded = lapply(expression_files, function(x)NULL)

# ref_gr = rtracklayer::import.gff("~/gencode.v22.annotation.gtf.gz", format = "gtf", feature.type = "gene")
# names(ref_gr) = sub("\\..+", "", ref_gr$gene_id)
# names(ref_gr) = ref_gr$gene_id