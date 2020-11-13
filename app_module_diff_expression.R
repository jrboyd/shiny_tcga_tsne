library("DESeq2")
library("BiocParallel")
library(ggplot2)
NCORES = 10

run_group.fast = function(exp_mat, a, b, min_fraction = .5){
    browser()
    dt = as.data.table(reshape2::melt(exp_mat[, c(a, b)]))
    setnames(dt, c("Var1", "Var2", "value"), c("gene_name", "id", "expression"))
    dt[, group := ifelse(id %in% a, "A", "B")]
    dt
}

run_DE.fast = function(dt){
    # system.time({
    # dt[, .(mean_val = mean(expression), 
    #        sd_val = sd(expression), 
    #        q25_val = quantile(expression, .25),
    #        q50_val = quantile(expression, .50),
    #        q75_val = quantile(expression, .75)
    # ), .(gene_name, group)]
    # })
    
    # system.time({
    p_dt = dt[, .(mean_val = mean(expression)
    ), .(gene_name, group)]
    # })
    p_dt = dcast(p_dt, gene_name~group, value.var = "mean_val")
    
    ps = 100
    p_dt[, lg2_fc := log2((B + ps) / (A + ps))]
    p_dt[, lg2_min := log2(pmin(A, B)+ps)]
    
    p_dt
}

run_DE = function(exp_mat, a, b, min_fraction = .5){
    sample_info = data.table(id = c(a, b), group = factor(c(rep("A", length(A)), rep("B", length(b)))))
    
    
    register(MulticoreParam(NCORES))
    
    input_mat = round(exp_mat[, c(a, b)])
    
    #filter based on percent of samples where genes are detected
    rf_a = apply(input_mat[,a], 1, function(x)sum(x>0)/length(x))
    rf_b = apply(input_mat[,b], 1, function(x)sum(x>0)/length(x))
    k = rf_a> min_fraction | rf_b > min_fraction
    input_mat = input_mat[k,]
    dds <- DESeqDataSetFromMatrix(countData = input_mat,
                                  colData = sample_info,
                                  design = ~ group)
    dds = DESeq(dds)#, parallel = TRUE)
    res <- as.data.frame(results(dds))
    res.sig = subset(res, padj < .01 & baseMean > 10 & (log2FoldChange) > 3)
    dim(res.sig)
    
    if(FALSE){    
        gplots::heatmap.2(log10(exp_mat[rownames(res.sig), c(a, b)]+1), trace = "n", 
                          Colv = FALSE, dendrogram = "row",
                          ColSideColors = c(rep("purple", length(a)), rep("orange", length(b))), 
                          key.title = "", density.info = "none", col = viridis::viridis(20), margins = c(16, 10))
        
        tmp = log2(exp_mat[rownames(res.sig), c(a, b)]+1)
        tmp = melt(as.data.table(tmp, id.col = "id", keep.rownames = "id"), id.vars = "id")
        tmp[, group := ifelse(variable %in% a, "A", "B")]
        ggplot(tmp[id %in% sample(unique(id), min(9, length(unique(id))))], aes(x = id, y = value, fill = group)) +
            geom_boxplot() +
            facet_wrap(~id, scales = "free")
    }
    res.sig
}

.test_run_DE = function(){
    #testing vars
    exp_dt = fread("data/BRCA_TCGA_expression.tiny.csv", sep = ",", header = TRUE)
    exp_dt = melt(exp_dt, id.vars = "gene_name")
    exp_dt = exp_dt[, .(value = sum(value)), .(gene_name, variable)]
    # exp_dt = unique(exp_mat$gene_name[duplicated(exp_mat$gene_name)])
    exp_dt = dcast(exp_dt, gene_name~variable, value.var = "value")
    exp_mat = as.matrix(exp_dt[, -1])
    rownames(exp_mat) = exp_dt$gene_name
    
    a = sample(colnames(exp_mat), 20)
    b = setdiff(colnames(exp_mat), a)
    
    run_DE(exp_mat, a, b)
}
