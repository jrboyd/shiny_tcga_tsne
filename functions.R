library(data.table)
library(openxlsx)

parse_gl = function(txt){
    gl = strsplit(txt, "[, \n]")[[1]]
    gl = gl[gl != ""]
    toupper(gl)
}

load_expression = function(f){
    dat = fread(f)
    mat = as.matrix(dat[,-1])
    rownames(mat) = dat$gene_name
    mat
}

#parse any csv, txt, or xlsx file
decide_parse_FUN = function(f, name = f){
    stopifnot(is.character(f))
    stopifnot(length(f) == 1)
    ext = sub(".+\\.", "", name)
    parse_FUN = tryCatch(
        {
            # FUN = NULL
            FUN = switch (ext,
                          csv = parse_csv,
                          txt = parse_txt,
                          xlsx = parse_xlsx
            )
            if(is.null(FUN)){
                stop()
            }
            FUN
        },
        error = function(e){
            # message('Could not determine file type from extension for ', f)
            # NULL
            paste0('Could not determine file type from extension for ', f)
        })
    if(is.character(parse_FUN)) return(parse_FUN)
    #parse_FUN must take a single character argument that is the path to a file
    #it must return a named list of data.tables/data.frames
    
    tryCatch(
        {
            out = suppressWarnings(parse_FUN(f, name))
        },
        error = function(e){
            
        })
    stopifnot(is.list(out))
    stopifnot(!is.null(names(out)))
    stopifnot(sapply(out, is.data.frame))
    return(out)
}

parse_csv = function(f, name){
    out = list(fread(f))
    names(out) = name
    out
}
parse_txt = function(f, name){
    out = list(fread(f))
    names(out) = name
    out
}
parse_xlsx = function(f, name){
    sn = openxlsx::getSheetNames(f)
    if(length(sn) == 1){
        out = list(openxlsx::read.xlsx(xlsxFile = f))
        names(out) = name
        out
    }else{
        names(sn) = sn
        lapply(sn, openxlsx::read.xlsx, xlsxFile = f)    
    }
}

#' get class for every column in df
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
get_column_classes = function(df){
    sapply(df, class)
}

#' locate the first column in df containing valid_genes
#' enforces gene_name as column name
#'
#' @param df 
#' @param valid_genes 
#'
#' @return
#' @export
#'
#' @examples
locate_genes_in_df = function(df, valid_genes){
    to_check = colnames(df)[get_column_classes(df) %in% c("character", "factor")]    
    gene_col = NULL
    for(tc in to_check){
        if(any(toupper(df[[tc]]) %in% valid_genes)){
            gene_col = tc
            break
        }
    }
    message(gene_col)
    if(!is.null(gene_col)){
        setnames(df, gene_col, "gene_name")
    }else{
        if(valid_genes[1] == "PASTE")
        warning("could not locate valid genes in provided data.")
    }
    df
}

#' determine which genes in gene_name are valid
#'
#' @param df 
#' @param valid_genes 
#' @param valid_var 
#'
#' @return
#' @export
#'
#' @examples
validate_genes_in_df = function(df, valid_genes, valid_var = "valid"){
    # stopifnot("gene_name" %in% colnames(df))
    if(is.null(df[["gene_name"]])){
        df$gene_name = "missing"   
    }
    df[[valid_var]] = df$gene_name %in% valid_genes
    var = c("gene_name", "valid")
    
    df = df[, c(var, setdiff(colnames(df), var))]    
    df
}


sampleCap = function(x, n = 500){
    n = min(n, length(unique(x)))
    out = sample(unique(x), n)
    if(is.factor(out)) out = as.character(out)
    out
}

km_clust = function(tsne_res, k = 5, id_var = "patient_id", x_var = "x", y_var = "y", nsamp = Inf){
    if(k < 2){
        k = 2
        warning("increasing nn to 2")
    }
    if(k > nrow(tsne_res)/2){
        k = round(nrow(tsne_res)/2)
        warning("decreasing k to ", k)
    }
    mat = t(as.matrix(tsne_res[, c(x_var, y_var), with = FALSE]))
    colnames(mat) = tsne_res[[id_var]]
    mat = mat[, sampleCap(seq(ncol(mat)), nsamp)]
    
    km_res = kmeans(t(mat), centers = k)
    tsne_res$cluster_id = paste("cluster", km_res$cluster[tsne_res[[id_var]]])
    tsne_res
}

h_clust = function(tsne_res, n_clust = 5, id_var = "patient_id", x_var = "x", y_var = "y", nsamp = Inf){
    if(n_clust < 2){
        n_clust = 2
        warning("increasing n_clust to 2")
    }
    if(n_clust > nrow(tsne_res)/2){
        n_clust = round(nrow(tsne_res)/2)
        warning("decreasing n_clust to ", n_clust)
    }
    mat = t(as.matrix(tsne_res[, c(x_var, y_var), with = FALSE]))
    colnames(mat) = tsne_res[[id_var]]
    mat = mat[, sampleCap(seq(ncol(mat)), nsamp)]
    
    h_res = hclust(dist(t(mat)))
    
    
    tsne_res$cluster_id = paste("cluster", cutree(h_res, n_clust)[tsne_res[[id_var]]])
    tsne_res
}

#from seqtsne
nn_clust = function(tsne_res, nn = 100, auto_nn_fraction = 5, id_var = "patient_id", x_var = "x", y_var = "y", nsamp = Inf){
    if(nn < 2){
        nn = 2
        warning("increasing nn to 2")
    }
    if(nn > nrow(tsne_res)/2){
        nn = round(nrow(tsne_res)/2)
        warning("decreasing nn to ", nn)
    }
    if(auto_nn_fraction*nn > nrow(tsne_res)){
        warning("Automatically reducing nearest neighbors")
        nn = nrow(tsne_res)/auto_nn_fraction
        if(nn < 1){
            stop("not enough samples to cluster. ", nrow(tsne_res), " samples submitted.")
        }
    }
    mat = t(as.matrix(tsne_res[, c(x_var, y_var), with = FALSE]))
    colnames(mat) = tsne_res[[id_var]]
    mat = mat[, sampleCap(seq(ncol(mat)), nsamp)]
    knn.info <- RANN::nn2(t(mat), k=nn)
    knn <- knn.info$nn.idx
    colnames(knn) = c(id_var, paste0("V", seq(nn-1)))
    knn = as.data.table(knn)
    mknn = melt(knn, id.vars = id_var)
    
    # adj <- matrix(0, ncol(mat), ncol(mat))
    # rownames(adj) <- colnames(adj) <- colnames(mat)
    # for(i in seq_len(ncol(mat))) {
    #     adj[i,colnames(mat)[knn[i,]]] <- 1
    # }
    ADJ = Matrix::Matrix(0, ncol(mat), ncol(mat))
    ADJ[cbind(mknn[[id_var]], mknn$value)] = 1
    rownames(ADJ) = colnames(mat)
    colnames(ADJ) = colnames(mat)
    g <- igraph::graph.adjacency(ADJ, mode="undirected")
    g <- igraph::simplify(g) ## remove self loops
    # V(g)$color <- rainbow(G)[group[names(V(g))]] ## color nodes by group
    # plot(g, vertex.label=NA)
    km <- igraph::cluster_walktrap(g)
    ## community membership
    com <- km$membership
    names(com) <- km$names
    com_dt = data.table(V1 = names(com), cluster_id = paste("cluster", com))
    setnames(com_dt, "V1", id_var)
    
    p_dt = merge(tsne_res, com_dt, by = id_var)
    
    # p_dt[, coms := paste("cluster", com)]
    # ggplot(p_dt, aes(x = tx, y = ty, color = coms)) +
    #     annotate("point", x  = p_dt$tx, y = p_dt$ty, size = .2) +
    #     geom_point(size = .5) +
    #     facet_wrap("coms")
    
    # p = ggplot(p_dt, aes_string(x = x_var, y = y_var, color = "cluster_id")) +
    #     labs(color = "cluster_id") +
    #     # annotate("point", x  = p_dt$tx, y = p_dt$ty, size = .2) +
    #     geom_point(size = .5) #+
    # facet_wrap("coms")
    # return(list(data = p_dt, plot = p))
    
    p_dt
}

#this is not used anywhere
reinit_data = function(sel){
    browser()
    sel = input$sel_data
    #meta data
    if(!sel %in% names(clinical_loaded)){
        stop(sel, " not found in loaded metadata.")
    }
    if(is.null(clinical_loaded[[sel]])){
        clinical_loaded[[sel]] = fread(clinical_files[[sel]])
    }
    meta_data(clinical_loaded[[sel]])
    
    #expression data
    if(!sel %in% dataset_names){
        stop(sel, " not found in expression_files.")
    }
    if(is.null(expression_loaded[[sel]])){
        expression_loaded[[sel]] = load_expression(expression_files[[sel]])
    }
    dat = expression_loaded[[sel]]
    exp_dt = as.data.table(dat)
    exp_dt$gene_name = rownames(dat)
    exp_dt = melt(exp_dt, id.vars = "gene_name")
    exp_dt = exp_dt[, .(value = max(value)), .(gene_name, variable)]
    # exp_dt = unique(exp_mat$gene_name[duplicated(exp_mat$gene_name)])
    exp_dt = dcast(exp_dt, gene_name~variable, value.var = "value")
    exp_mat = as.matrix(exp_dt[, -1])
    rownames(exp_mat) = exp_dt$gene_name
    
    stopifnot(!any(duplicated(rownames(exp_mat))))
    tcga_data(exp_mat)
    
    #reset downstream
    
}

library("DESeq2")
library("BiocParallel")
library(ggplot2)
NCORES = 10

run_group.fast = function(exp_mat, a, b, min_fraction = .5){
    # browser()
    if("reactiveVal" %in% class(a)){
        a = a()
    }
    if("reactiveVal" %in% class(b)){
        b = b()
    }
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

run_survival_by_patients = function(surv_dt, impacted_patients, legend_title = "impacted"){
    surv_dt[, impacted := patient_id %in% impacted_patients]
    
    res = TCGAbiolinks::TCGAanalyze_survival(as.data.frame(surv_dt), 
                               "impacted", legend = legend_title,
                               height = 10, 
                               width=10, 
                               filename = NULL)
    pval = res$plot$layers[[5]]$geom_params$label
    pval = as.numeric(sub("p [=<>] ", "", pval))
    res$pval = pval
    res
}

run_DE = function(exp_mat, a, b, min_fraction = .5, on_windows = FALSE, max_padj = .01, min_baseMean = 10, min_log2FoldChange = 3){
    if("reactiveVal" %in% class(a)){
        a = a()
    }
    if("reactiveVal" %in% class(b)){
        b = b()
    }
    # saveRDS(list(exp_mat = exp_mat, a = a, b = b, min_fraction = min_fraction ), "diff_data.Rds")
    sample_info = data.table(id = c(a, b), group = factor(c(rep("A", length(a)), rep("B", length(b)))))
    
    if(on_windows){
        register(SnowParam(NCORES))    
    }else{
        register(MulticoreParam(NCORES))
    }
    
    
    
    input_mat = round(exp_mat[, c(a, b)])
    
    #filter based on percent of samples where genes are detected
    rf_a = apply(input_mat[,a], 1, function(x)sum(x>0)/length(x))
    rf_b = apply(input_mat[,b], 1, function(x)sum(x>0)/length(x))
    k = rf_a> min_fraction | rf_b > min_fraction
    input_mat = input_mat[k,]
    rm = apply(input_mat, 1, max)
    
    input_mat = input_mat[rm<5e7,]
    library(BiocFileCache)
    bfc = BiocFileCache::BiocFileCache()
    bfcif = ssvRecipes::bfcif
    res <- bfcif(bfc, paste("DESeq2_v2", digest::digest(list(input_mat, sample_info))), function(){
        dds = DESeqDataSetFromMatrix(countData = input_mat,
                               colData = sample_info,
                               design = ~ group)
        dds = DESeq(dds)#, parallel = TRUE)
        as.data.frame(results(dds))
    })
 

    res.sig = subset(res, padj < max_padj & baseMean > min_baseMean & abs(log2FoldChange) > min_log2FoldChange)
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

# run_tsne = function(exp_mat, genes, perplexity = 20){
#     library(BiocFileCache)
#     bfc = BiocFileCache::BiocFileCache()
#     bfcif = ssvRecipes::bfcif
#     expression_matrix = exp_mat[genes,]
#     tsne_patient = ssvRecipes::bfcif(bfc, digest::digest(list(expression_matrix, perplexity)), function(){
#         Rtsne::Rtsne(t(expression_matrix), 
#                      num_threads = 20, 
#                      check_duplicates = FALSE,
#                      perplexity = perplexity)    
#     })
#     
#     tsne_df = as.data.table(tsne_patient$Y)
#     colnames(tsne_df) = c("x", "y")
#     tsne_df$sample_id = colnames(expression_matrix)
#     tsne_df
# }

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