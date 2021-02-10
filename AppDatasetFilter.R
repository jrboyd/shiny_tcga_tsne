setClass("AppDatasetFilter", 
         representation(
             var = "character",
             inverse = "logical"
         )
)

setClass("AppDatasetFilterNumeric", 
         contains = "AppDatasetFilter", 
         representation(
             min_val = "numeric", 
             max_val = "numeric",
             min_inclusive = "logical",
             max_inclusive = "logical"
         )
)

setClass("AppDatasetFilterCharacter", 
         contains = "AppDatasetFilter", 
         representation(
             one_of = "character"
         )
)

AppDatasetFilterNumeric = function(var, 
                                   inverse = FALSE, 
                                   min_val = -Inf, 
                                   max_val = Inf,
                                   min_inclusive = FALSE,
                                   max_inclusive = FALSE){
    new("AppDatasetFilterNumeric", 
        var = var, 
        inverse = inverse,
        min_val = min_val,
        max_val = max_val, 
        min_inclusive = min_inclusive, 
        max_inclusive = max_inclusive)
}

AppDatasetFilterCharacter = function(var, 
                                   inverse = FALSE, 
                                   one_of = character()){
    new("AppDatasetFilterCharacter", 
        var = var, 
        inverse = inverse,
        one_of = one_of)
}

setMethod("show", "AppDatasetFilter", function(object){
    msg = "NYI"
    return(msg)
})

setMethod("show", "AppDatasetFilterCharacter", function(object){
    msg = paste(object@var, "is not filtered.")
    if(length(object@one_of) > 0){
        msg = paste0("one of (", paste(object@one_of, collapse = ", "),")")
        if(object@inverse){
            msg = paste0(object@var, " must NOT be ", msg)
        }else{
            msg = paste0(object@var, " must be ", msg)
        }
    }
    message(msg)
    return(msg)
})

setMethod("show", "AppDatasetFilterNumeric", function(object){
    msg = paste(object@var, "is not filtered.")
    min_sym = ifelse(object@min_inclusive, ">=", ">")
    max_sym = ifelse(object@max_inclusive, "<=", "<")
    if(is.finite(object@min_val) & is.finite(object@max_val)){
        msg = paste(object@var, min_sym, object@min_val, "and", max_sym, object@max_val)
    }else if(is.finite(object@min_val)){
        msg = paste(object@var, min_sym, object@min_val)
    }else if(is.finite(object@max_val)){
        msg = paste(object@var, max_sym, object@max_val)
    }
    if(object@inverse){
        msg = paste0("not(", msg, ")")
    }
    message(msg)
    return(msg)
})

AppDatasetFilterNumeric("num")
AppDatasetFilterNumeric("num", min_val = 2)
AppDatasetFilterNumeric("num", min_val = 2, max_val = 4, min_inclusive = TRUE, max_inclusive = TRUE)
AppDatasetFilterNumeric("num", min_val = 2, max_val = 4, inverse = TRUE)

AppDatasetFilterCharacter("char")
AppDatasetFilterCharacter("char", one_of = c("a", "c"))
AppDatasetFilterCharacter("char", one_of = c("a", "c", "asdf"), inverse = TRUE)

