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

km_clust = function(tsne_res, k = 5, id_var = "submitter_id", x_var = "x", y_var = "y", nsamp = Inf){
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

h_clust = function(tsne_res, n_clust = 5, id_var = "submitter_id", x_var = "x", y_var = "y", nsamp = Inf){
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
nn_clust = function(tsne_res, nn = 100, auto_nn_fraction = 5, id_var = "submitter_id", x_var = "x", y_var = "y", nsamp = Inf){
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
# ex_files = dir("example_data", full.names = TRUE)
# # funs = lapply(ex_files, decide_parse_FUN)
# lapply(ex_files, function(f){
#     message(f)
#     decide_parse_FUN(f)
# })
# sapply(paste0(ex_files, "bad"), decide_parse_FUN)
