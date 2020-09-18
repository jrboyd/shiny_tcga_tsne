library(genefu)

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

gene_lists = list(
    PAM50 = make_pam50()    
)



