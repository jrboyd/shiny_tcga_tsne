setClass("Operation",
         representation(
             var = "character"
         )
)

setClass("Filter", 
         representation(
             var = "character",
             inverse = "logical",
             drop_NA = "logical"
         ), 
         contains = "Operation"
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

setGeneric("applyOperation", function(df, object){standardGeneric("applyOperation")})

setMethod("applyOperation", c("data.frame", "Filter"), function(df, object){
    stop("applyOperation: use FilterNumeric or FilterCharacter")
})

setMethod("applyOperation", c("data.frame", "FilterNumeric"), function(df, object){
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
setMethod("applyOperation", c("data.frame", "FilterCharacter"), function(df, object){
    .character_filter(df,
                      var = object@var,
                      one_of = object@one_of,
                      inverse = object@inverse,
                      drop_NA = object@drop_NA)
})

setClass("Transform", 
         representation(),
         contains = "Operation"
)

setMethod("applyOperation", c("data.frame", "Transform"), function(df, object){
    stop("applyOperation: use class that inherits Transform")
})
########
setClass("TransformZscore", 
         contains = "Transform")

TransformZscore = function(var){
    new("TransformZscore", var = var)
}

setMethod("as.character", c("TransformZscore"), function(x){
    paste(x@var, "z-score")
})
setMethod("show", c("TransformZscore"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformZscore"), function(df, object){
    vals = df[[object@var]]
    vals = (vals - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)
    df[[object@var]] = vals
    df
})
########
setClass("TransformNegate", 
         contains = "Transform")

TransformNegate = function(var){
    new("TransformNegate", var = var)
}

setMethod("as.character", c("TransformNegate"), function(x){
    paste(x@var, "z-score")
})
setMethod("show", c("TransformNegate"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformNegate"), function(df, object){
    vals = df[[object@var]]
    vals = (vals - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)
    df[[object@var]] = vals
    df
})

########
setClass("TransformLog", 
         contains = "Transform", 
         representation(
             power = "numeric"
         )
)

TransformLog = function(var, power = 10){
    new("TransformLog", var = var, power = power)
}

setMethod("as.character", c("TransformLog"), function(x){
    paste(x@var, paste0("log", x@power))
})
setMethod("show", c("TransformLog"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformLog"), function(df, object){
    vals = df[[object@var]]
    vals = log(vals, base = object@power)
    df[[object@var]] = vals
    df
})

########
setClass("TransformCrop", 
         contains = "Transform",
         representation(
             min_val = "numeric",
             max_val = "numeric"
         )
)

TransformCrop = function(var, min_val = -3, max_val = 3){
    new("TransformCrop", var = var, min_val = min_val, max_val = max_val)
}

setMethod("as.character", c("TransformCrop"), function(x){
    msg = paste(x@var, "crop")
    if(is.finite(x@min_val)){
        msg = paste(msg, "above", x@min_val)
    }
    if(is.finite(x@max_val)){
        msg = paste(msg, "below", x@max_val)
    }
    paste(msg)
})
setMethod("show", c("TransformCrop"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformCrop"), function(df, object){
    vals = df[[object@var]]
    vals[vals < object@min_val] = object@min_val
    vals[vals > object@max_val] = object@max_val
    df[[object@var]] = vals
    df
})

########
setClass("TransformRescale", 
         contains = "Transform",
         representation(
             min_val = "numeric",
             max_val = "numeric"
         )
)

TransformRescale = function(var, min_val = 0, max_val = 1){
    stopifnot(is.finite(min_val))
    stopifnot(is.finite(max_val))
    new("TransformRescale", var = var, min_val = min_val, max_val = max_val)
}

setMethod("as.character", c("TransformRescale"), function(x){
    paste(x@var, "rescale", x@min_val, "to", x@max_val )
})
setMethod("show", c("TransformRescale"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformRescale"), function(df, object){
    vals = df[[object@var]]
    vals = scales::rescale(vals, c(object@min_val, object@max_val))
    df[[object@var]] = vals
    df
})

########
setClass("TransformBin", 
         contains = "Transform",
         representation(
             n_bins = "numeric"
         )
)

TransformBin = function(var, n_bins = 5){
    new("TransformBin", var = var, n_bins = n_bins)
}

setMethod("as.character", c("TransformBin"), function(x){
    paste(x@var, "binned", x@n_bins)
})
setMethod("show", c("TransformBin"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformBin"), function(df, object){
    vals = df[[object@var]]
    vals = cut(vals, breaks = object@n_bins)
    df[[object@var]] = vals
    df
})

########
setClass("TransformRank", 
         contains = "Transform",
         representation(
             decreasing = "logical"
         )
)

TransformRank = function(var, decreasing = FALSE){
    new("TransformRank", var = var, decreasing = decreasing)
}

setMethod("as.character", c("TransformRank"), function(x){
    paste(x@var, "rank", ifelse(x@decreasing, "decreasing", "increasing"))
})
setMethod("show", c("TransformRank"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformRank"), function(df, object){
    vals = df[[object@var]]
    if(object@decreasing){
        vals = frank(-vals, ties.method = "first")
    }else{
        vals = frank(vals, ties.method = "first")
    }
    df[[object@var]] = vals
    df
})

########
setClass("TransformSort", 
         contains = "Transform",
         representation(
             decreasing = "logical"
         )
)

TransformSort = function(var, decreasing = FALSE){
    new("TransformSort", var = var, decreasing = decreasing)
}

setMethod("as.character", c("TransformSort"), function(x){
    paste(x@var, "rank", ifelse(x@decreasing, "decreasing", "increasing"))
})
setMethod("show", c("TransformSort"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformSort"), function(df, object){
    vals = df[[object@var]]
    o = order(vals, decreasing = object@decreasing)
    df = df[o,]
    df
})

########
setClass("TransformCharacterRemap", 
         contains = "Transform",
         representation(
             mapping = "character",
             NA_fill = "character"
         )
)

TransformCharacterRemap = function(var, mapping, NA_fill = "NA"){
    new("TransformCharacterRemap", var = var, mapping = mapping, NA_fill = NA_fill)
}

setMethod("as.character", c("TransformCharacterRemap"), function(x){
    paste(x@var, "remap", length(unique(names(x@mapping))), "to", length(unique(x@mapping)))
})
setMethod("show", c("TransformCharacterRemap"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformCharacterRemap"), function(df, object){
    
    vals = df[[object@var]]
    missed = setdiff(vals, names(object@mapping))
    names(missed) = missed
    dict = c(missed, object@mapping)
    vals = dict[vals]
    names(vals) = NULL
    df[[object@var]] = vals
    df
})

########
setClass("TransformCopy", 
         contains = "Transform",
         representation(
             new_var = "character"
         )
)

TransformCopy = function(var, new_var){
    new("TransformCopy", var = var, new_var = new_var)
}

setMethod("as.character", c("TransformCopy"), function(x){
    paste(x@var, "copy", x@var, "to", x@new_var)
})
setMethod("show", c("TransformCopy"), function(object){
    message(as.character(object))  
})
setMethod("applyOperation", c("data.frame", "TransformCopy"), function(df, object){
    vals = df[[object@var]]
    df[[object@new_var]] = vals
    df
})

#informal tests
if(FALSE){
    {
        library(data.table)
        df = data.table::fread("installed_datasets/BRCA_tiny/clinical.csv")
        
        
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
        
        
        
        #individually
        applyOperation(df, num_filter)
        applyOperation(df, char_filter)
        applyOperation(df, char_filter2)
        
        #nested
        applyOperation(applyOperation(df, num_filter), char_filter)
        
        #magrittr
        library(magrittr)
        df %>% 
            applyOperation(num_filter) %>% 
            applyOperation(char_filter) %>%
            applyOperation(char_filter2)
        
        #programmatic loop
        all_filters = list(num_filter, char_filter, char_filter2)
        df_loop = df
        for(fil in all_filters){
            df_loop = applyOperation(df_loop, fil)
        }
        df_loop
        
        ###
        tst = TransformZscore("year_of_birth")
        as.character(tst)
        x = tst
        object = tst
        applyOperation(df, tst)
        
        tst = TransformNegate("year_of_birth")
        as.character(tst)
        x = tst
        object = tst
        applyOperation(df, tst)
        
        tst = TransformLog("year_of_birth")
        as.character(tst)
        x = tst
        object = tst
        applyOperation(df, tst)
        
        tst = TransformCrop("year_of_birth")
        as.character(tst)
        x = tst
        object = tst
        applyOperation(df, tst)
        
        tst = TransformRescale("year_of_birth")
        as.character(tst)
        x = tst
        object = tst
        applyOperation(df, tst)
        
        tst = TransformBin("year_of_birth")
        as.character(tst)
        x = tst
        object = tst
        applyOperation(df, tst) %>% applyOperation(FilterCharacter("year_of_birth", one_of = c("(1.95e+03,1.97e+03]")))
        
        tst = TransformRank("year_of_birth", decreasing = TRUE)
        as.character(tst)
        x = tst
        object = tst
        x = df$year_of_birth
        y = applyOperation(df, tst)$year_of_birth
        y
        plot(x,y, xlab = "year", ylab = "rank")
        title(as.character(tst))
        
        tst = TransformSort("year_of_birth", decreasing = TRUE)
        as.character(tst)
        x = tst
        object = tst
        head(applyOperation(df, tst), n = 20)
        
        tst = TransformCharacterRemap("race", mapping = c("white" = "white or asian", "asian" = "white or asian"))
        as.character(tst)
        x = tst
        object = tst
        applyOperation(df, tst)
        table(applyOperation(df, tst)$race)
        
        tst = TransformCopy("year_of_birth", "year_of_birth_ranked")
        as.character(tst)
        x = tst
        object = tst
        df2 = applyOperation(df, tst) %>%
            applyOperation(TransformRank("year_of_birth_ranked"))
        plot(df2$year_of_birth, df2$year_of_birth_ranked)
    }   
}