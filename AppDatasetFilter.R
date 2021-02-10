setClass("AppDatasetFilter", 
         representation(
             var = "character",
             inverse = "logical",
             drop_NA = "logical"
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
                                   max_inclusive = FALSE, 
                                   drop_NA = TRUE){
    new("AppDatasetFilterNumeric", 
        var = var, 
        inverse = inverse,
        min_val = min_val,
        max_val = max_val, 
        min_inclusive = min_inclusive, 
        max_inclusive = max_inclusive,
        drop_NA = drop_NA)
}

AppDatasetFilterCharacter = function(var, 
                                     inverse = FALSE, 
                                     one_of = character(), 
                                     drop_NA = TRUE){
    new("AppDatasetFilterCharacter", 
        var = var, 
        inverse = inverse,
        one_of = one_of,
        drop_NA = drop_NA)
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
        if(object@min_val == object@max_val){
            msg = paste(object@var, "==", object@min_val)    
        }else{
            msg = paste(object@var, min_sym, object@min_val, "and", max_sym, object@max_val)    
        }
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

setMethod("applyFilter", c("data.frame", "AppDatasetFilter"), function(df, object){
    stop("applyFilter: use AppDatasetFilterNumeric or AppDatasetFilterCharacter")
})

setMethod("applyFilter", c("data.frame", "AppDatasetFilterNumeric"), function(df, object){
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
setMethod("applyFilter", c("data.frame", "AppDatasetFilterCharacter"), function(df, object){
    .character_filter(df,
                      var = object@var,
                      one_of = object@one_of,
                      inverse = object@inverse,
                      drop_NA = object@drop_NA)
})



AppDatasetFilterNumeric("num")


AppDatasetFilterNumeric("num", min_val = 2, max_val = 2)
AppDatasetFilterNumeric("num", min_val = 2, max_val = 2, inverse = TRUE)
AppDatasetFilterNumeric("num", min_val = 2, max_val = 4, min_inclusive = TRUE, max_inclusive = TRUE)
AppDatasetFilterNumeric("num", min_val = 2, max_val = 4, inverse = TRUE)

AppDatasetFilterCharacter("char")
AppDatasetFilterCharacter("char", one_of = c("a", "c"))
AppDatasetFilterCharacter("char", one_of = c("a", "c", "asdf"), inverse = TRUE)





object = AppDatasetFilterNumeric("days_to_last_follow_up", min_val = 50, max_val = 100)
var = object@var
min_val = object@min_val
max_val = object@max_val
min_inclusive = object@min_inclusive
max_inclusive = object@max_inclusive
inverse = object@inverse
drop_NA = object@drop_NA

num_filter = object


char_filter = AppDatasetFilterCharacter("ethnicity", one_of = c("not reported"), inverse = TRUE)
char_filter2 = AppDatasetFilterCharacter("vital_status", one_of = c("alive"))

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
