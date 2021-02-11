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

setGeneric("addOperation", function(ops_list, new_op){standardGeneric("addOperation")})

setMethod("addOperation", c("OperationList", "Operation"), function(ops_list, new_op){
    ops_list@ops = c(ops_list@ops, list(new_op))
    ops_list
})


library(magrittr)
library(data.table)
df = data.table::fread("installed_datasets/BRCA_tiny/clinical.csv")


#operations
num_filter = FilterNumeric("days_to_last_follow_up", min_val = 50, max_val = 100)
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

ops = OperationList()
ops = ops %>%
    addOperation(x_copy) %>%
    addOperation(x_rank) %>%
    addOperation(x_bin) %>%
    addOperation(num_filter) %>%
    addOperation(char_filter) %>%
    addOperation(char_filter2)


ops
