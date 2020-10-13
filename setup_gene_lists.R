library(genefu)
library(data.table)

gene_list_file = "data/BRCA_gene_signatures.Rds"
if(file.exists(gene_list_file)){
    message("loading cached gene lists")
    gene_lists = readRDS(gene_list_file)
}else{
    make_pam50 = function(){
        data("pam50.robust", envir = environment())
        pam50 = pam50.robust$centroids
        
        subGN = function(match, replace, df = pam50){
            rownames(df)[rownames(df) == match] = replace
            df
        }
        pam50 = subGN("CDCA1", "NUF2")
        pam50 = subGN("KNTC2", "NDC80")
        pam50 = subGN("ORC6L", "ORC6")
        rownames(pam50)
    }
    
    ref_gr = rtracklayer::import.gff("~/gencode.v22.annotation.gtf.gz")
    valid_gn = unique(ref_gr$gene_name)
    
    valid_gn[grepl("AA555029", valid_gn)]
    
    fix_gene_names = function(x){
        #these could be in a file
        x = sub("^CDCA1$", "NUF2", x)
        x = sub("^KNTC2$", "NDC80", x)
        x = sub("^ORC6L$", "ORC6", x)
        x = sub("^CTSL2$", "CTSV", x)
        x = sub("^GUS$", "GUSB", x)
        x = sub("^HER2$", "ERBB2", x)
        x = sub("^Ki67$", "MKI67", x)
        x = sub("^RPLPO$", "RPLP0", x)
        x = sub("^STK15$", "AURKA", x)
        x = sub("^TRFC$", "TFRC", x)
        x = sub("^AYTL2$", "LPCAT1", x)
        x = sub("^JHDM1D$", "KDM7A", x)
        x = sub("^LGP2$", "DHX58", x)
        x = sub("^PECI$", "ECI2", x)
        x = sub("^QSCN6L1$", "QSOX2", x)
        x = sub("^ZNF533$", "ZNF385B", x)
        x = sub("^IL17BR$", "IL17RB", x)
        x = sub("^HTF9C$", "TRMT2A", x)
        x = sub("^C20orf46$", "TMEM74B", x)
        x = sub("^C16orf61$", "CMC2", x)
        x = sub("^C9orf30$", "MSANTD3", x)
        x = sub("^AA555029_RC$", "LOC286052", x)
        x = sub("^LOC286052$", "TMEM65", x)
        x
    }
    
    gene_lists = list(
        PAM50 = make_pam50()    
    )
    
    gl_dt = fread("data/BRCA_gene_signatures.Table_1.PMC6131478.csv")
    setnames(gl_dt, c("Oncotype DX", "Prosigna/PAM50"), c("OncotypeDX", "PAM50"))
    gl_dt = lapply(as.list(gl_dt), function(x)x[x != ""])
    gl_dt = lapply(gl_dt, function(x){
        fix_gene_names(x)
    })
    seqsetvis::ssvFeatureVenn(list(genefu = gene_lists$PAM50, pmc = gl_dt$PAM50))
    
    lapply(gl_dt, function(x){
        setdiff(x, valid_gn)
    })
    
    gl_dt = lapply(gl_dt, function(x){
        intersect(x, valid_gn)
    })
    gl_dt
    saveRDS(gl_dt, gene_list_file)
    gene_lists = gl_dt
}

gene_lists$all_signatures = unique(unlist(gene_lists))
