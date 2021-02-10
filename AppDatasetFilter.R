setClass("Filter", 
         representation(
             var = "character",
             inverse = "logical",
             drop_NA = "logical"
         )
)

setClass("FilterNumeric", 
         contains = "Filter", 
         representation(
             min_val = "numeric", 
             max_val = "numeric",
             min_inclusive = "logical",
             max_inclusive = "logical"
         )
)

setClass("FilterCharacter", 
         contains = "Filter", 
         representation(
             one_of = "character"
         )
)

FilterNumeric = function(var, 
                                   inverse = FALSE, 
                                   min_val = -Inf, 
                                   max_val = Inf,
                                   min_inclusive = FALSE,
                                   max_inclusive = FALSE, 
                                   drop_NA = TRUE){
    new("FilterNumeric", 
        var = var, 
        inverse = inverse,
        min_val = min_val,
        max_val = max_val, 
        min_inclusive = min_inclusive, 
        max_inclusive = max_inclusive,
        drop_NA = drop_NA)
}

FilterCharacter = function(var, 
                                     inverse = FALSE, 
                                     one_of = character(), 
                                     drop_NA = TRUE){
    new("FilterCharacter", 
        var = var, 
        inverse = inverse,
        one_of = one_of,
        drop_NA = drop_NA)
}

setMethod("show", "Filter", function(object){
    msg = "NYI"
    return(msg)
})

setMethod("as.character", "FilterCharacter", function(x){
    msg = paste(x@var, "is not filtered.")
    if(length(x@one_of) > 0){
        msg = paste0("in (", paste(x@one_of, collapse = ", "),")")
        if(x@inverse){
            msg = paste0(x@var, " NOT ", msg)
        }else{
            msg = paste0(x@var, " is ", msg)
        }
    }
    return(msg)
})

setMethod("as.character", "FilterNumeric", function(x){
    msg = paste(x@var, "not filtered.")
    min_sym = ifelse(x@min_inclusive, ">=", ">")
    max_sym = ifelse(x@max_inclusive, "<=", "<")
    if(is.finite(x@min_val) & is.finite(x@max_val)){
        if(x@min_val == x@max_val){
            msg = paste(x@var, "==", x@min_val)    
        }else{
            msg = paste(x@var, min_sym, x@min_val, "and", max_sym, x@max_val)    
        }
    }else if(is.finite(x@min_val)){
        msg = paste(x@var, min_sym, x@min_val)
    }else if(is.finite(x@max_val)){
        msg = paste(x@var, max_sym, x@max_val)
    }
    
    if(x@inverse){
        msg = paste0("not(", msg, ")")
    }
    return(msg)
})

setMethod("show", "FilterCharacter", function(object){
    message(as.character(object))
})

setMethod("show", "FilterNumeric", function(object){
    message(as.character(object))
})


.numeric_filter = function(df, var, min_val, max_val, min_inclusive, max_inclusive, inverse, drop_NA){
    vals = df[[var]]
    if(is.finite(min_val) & is.finite(max_val)){
        if(min_val == max_val){
            k = vals == min_val
        }else{
            if(min_inclusive){
                k1 = vals >= min_val
            }else{
                k1 = vals > min_val
            }
            if(max_inclusive){
                k2 = vals <= max_val
            }else{
                k2 = vals < max_val
            }
            k = k1 & k2
        }
    }else if(is.finite(min_val)){
        if(min_inclusive){
            k = vals >= min_val
        }else{
            k = vals > min_val
        }
    }else if(is.finite(max_val)){
        if(max_inclusive){
            k = vals <= max_val
        }else{
            k = vals < max_val
        }
    }
    if(inverse){
        k = !k
    }
    if(drop_NA){
        k[is.na(k)] = FALSE
    }
    df[k,]
}

.character_filter = function(df, var, one_of, inverse, drop_NA){
    vals = df[[var]]
    k = vals %in% one_of
    if(inverse){
        k = !k
    }
    if(drop_NA){
        k[is.na(k)] = FALSE
    }
    df[k,]
}


setGeneric("applyFilter", function(df, object){standardGeneric("applyFilter")})

setMethod("applyFilter", c("data.frame", "Filter"), function(df, object){
    stop("applyFilter: use FilterNumeric or FilterCharacter")
})

setMethod("applyFilter", c("data.frame", "FilterNumeric"), function(df, object){
    .numeric_filter(df,
                    var = object@var,
                    min_val = object@min_val,
                    max_val = object@max_val,
                    min_inclusive = object@min_inclusive,
                    max_inclusive = object@max_inclusive,
                    inverse = object@inverse,
                    drop_NA = object@drop_NA
    )
})
setMethod("applyFilter", c("data.frame", "FilterCharacter"), function(df, object){
    .character_filter(df,
                      var = object@var,
                      one_of = object@one_of,
                      inverse = object@inverse,
                      drop_NA = object@drop_NA)
})



FilterNumeric("num")


FilterNumeric("num", min_val = 2, max_val = 2)
FilterNumeric("num", min_val = 2, max_val = 2, inverse = TRUE)
FilterNumeric("num", min_val = 2, max_val = 4, min_inclusive = TRUE, max_inclusive = TRUE)
FilterNumeric("num", min_val = 2, max_val = 4, inverse = TRUE)

FilterCharacter("char")
FilterCharacter("char", one_of = c("a", "c"))
FilterCharacter("char", one_of = c("a", "c", "asdf"), inverse = TRUE)





object = FilterNumeric("days_to_last_follow_up", min_val = 50, max_val = 100)
var = object@var
min_val = object@min_val
max_val = object@max_val
min_inclusive = object@min_inclusive
max_inclusive = object@max_inclusive
inverse = object@inverse
drop_NA = object@drop_NA

num_filter = object
char_filter = FilterCharacter("ethnicity", one_of = c("not reported"), inverse = TRUE)
char_filter2 = FilterCharacter("vital_status", one_of = c("alive"))

as.character(num_filter)
as.character(char_filter)

library(data.table)
df = fread("installed_datasets/BRCA_tiny/clinical.csv")

#individually
applyFilter(df, num_filter)
applyFilter(df, char_filter)
applyFilter(df, char_filter2)

#nested
applyFilter(applyFilter(df, num_filter), char_filter)

#magrittr
library(magrittr)
df %>% 
    applyFilter(num_filter) %>% 
    applyFilter(char_filter) %>%
    applyFilter(char_filter2)

#programmatic loop
all_filters = list(num_filter, char_filter, char_filter2)
df_loop = df
for(fil in all_filters){
    df_loop = applyFilter(df_loop, fil)
}
df_loop
