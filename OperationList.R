if(!exists("Operation")) source("Operation.R")

setClass("OperationList",
         representation(
             ops = "list"
         )
)

OperationList = function(){
    new("OperationList")
}

# setMethod("length", c("OperationList"), function(x){
#     length(x@ops)
# })

setMethod("length", "OperationList", function(x) length(x@ops))

setMethod("as.character", "OperationList", function(x){
    sapply(seq_along(x@ops), function(i){
        paste0(as.character(x@ops[[i]]))
    })
})

setMethod("show", "OperationList", function(object){
    txt = as.character(object)
    txt = sapply(seq_along(txt), function(i){
        paste0(i, ": ", txt[i])
    })
    message(paste(txt, collapse = "\n"))
})

setGeneric("addOperation", function(ops_list, new_op, i = NULL){standardGeneric("addOperation")})

setMethod("addOperation", c("OperationList", "Operation"), function(ops_list, new_op, i = NULL){
    stopifnot(i >= 1)
    stopifnot(i <= length(ops_list)+1)
    if(is.null(i)) i = length(ops_list) + 1
    if(i > length(ops_list)){
        ops_list@ops = c(ops_list@ops, list(new_op))    
    }else if(i <= 1){
        ops_list@ops = c(list(new_op), ops_list@ops)    
    }else{
        ops_list@ops = c(ops_list@ops[seq(1, i-1)], list(new_op), ops_list@ops[seq(i, length(ops_list))])    
    }
    ops_list
})

setGeneric("deleteOperation", function(ops_list, i){standardGeneric("deleteOperation")})

setMethod("deleteOperation", c("OperationList", "numeric"), function(ops_list, i){
    k = setdiff(seq_along(ops_list), i)
    ops_list@ops = ops_list@ops[k]
    ops_list
})

setGeneric("replaceOperation", function(ops_list, new_op, i){standardGeneric("replaceOperation")})

setMethod("replaceOperation", c("OperationList", "Operation", "numeric"), function(ops_list, new_op, i){
    stopifnot(length(i) == 1)
    stopifnot(i >= 1)
    stopifnot(i <= length(ops_list))
    if(is.null(i)) i = length(ops_list) + 1
    if(i == length(ops_list)){
        ops_list@ops = c(ops_list@ops[-i], list(new_op))    
    }else if(i <= 1){
        ops_list@ops = c(list(new_op), ops_list@ops[-i])    
    }else{
        ops_list@ops = c(ops_list@ops[seq(1, i-1)], list(new_op), ops_list@ops[seq(i+1, length(ops_list))])    
    }
    ops_list
})

setGeneric("executeOperationList", function(df, ops_list){standardGeneric("executeOperationList")})

setMethod("executeOperationList", c("data.frame", "OperationList"), function(df, ops_list){
    all_ops = ops_list@ops
    for(op in all_ops){
        df = applyOperation(df, op)
    }
    df
})


setMethod("[", c("OperationList", "ANY"), function(x, i){
    x@ops = x@ops[i]
    x
})

setMethod("[<-", c("OperationList", "ANY", "ANY", "Operation"), function(x, i, j, value){
    x@ops[i] = value
})

setMethod("[[", c("OperationList", "ANY"), function(x, i){
    x@ops[[i]]
})

setMethod("[[<-", c("OperationList", "ANY", "ANY", "Operation"), function(x, i, j, ..., value){
    x@ops[[i]] = value
})


if(FALSE){
    {
        library(magrittr)
        library(data.table)
        df = data.table::fread("installed_datasets/BRCA_tiny/clinical.csv")
        samp = data.table::fread("installed_datasets/BRCA_tiny/samples.csv")
        df = merge(df, samp[, .(sample_id, patient_id)], by = 'patient_id')
        
        
        exp = fread("installed_datasets/BRCA_tiny/expression.csv")
        
        vals = suppressWarnings({as.numeric(exp[gene_name == "TIMP1"])[-1]})
        names(vals) = colnames(exp)[-1]
        tst = ifelse(vals > mean(vals), "high", "low")
        
        k = df$sample_id %in% names(vals)
        setequal(df$sample_id, names(vals))
        all(k)
        
        append_exp = AppendNumeric("sample_id", "TIMP1_expression", vals)
        append_tst = AppendCharacter("sample_id", "TIMP1_highlow", tst)
        #operations
        num_filter = FilterNumeric("days_to_last_follow_up", min_val = 50, max_val = 1000)
        char_filter = FilterCharacter("ethnicity", one_of = c("not reported"), inverse = TRUE)
        char_filter2 = FilterCharacter("vital_status", one_of = c("alive"))
        
        
        x_copy = TransformCopy("year_of_birth", "year_of_birth_ranked")
        x_rank = TransformRank("year_of_birth_ranked", decreasing = TRUE)
        x_bin = TransformBin("year_of_birth_ranked")
        

        
        ops = OperationList()
        ops = addOperation(ops, x_copy)
        ops = addOperation(ops, x_rank)
        ops = addOperation(ops, x_bin)
        ops
        ops = addOperation(ops, x_bin, i = 2)
        ops
        deleteOperation(ops, c(2,4))
        deleteOperation(ops, c(-1))
        deleteOperation(ops, c(1, 3))
        replaceOperation(ops, num_filter, 2)
        
        ops = OperationList()
        ops = ops %>%
            addOperation(x_copy) %>%
            addOperation(x_rank) %>%
            addOperation(x_bin) %>%
            addOperation(num_filter) %>%
            addOperation(char_filter) %>%
            addOperation(char_filter2) %>%
            addOperation(append_exp) %>%
            addOperation(append_tst)
        
        
        ops  
        
        executeOperationList(df, ops[1:8])
    }
}

