df = fread("installed_datasets/BRCA_tiny/clinical.csv")


setClass("Transform", 
         representation(
             var = "character"
         )
)

setGeneric("applyTransform", function(df, object){standardGeneric("applyTransform")})
setMethod("applyTransform", c("data.frame", "Transform"), function(df, object){
    stop("applyTransform: use class that inherits Transform")
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
setMethod("applyTransform", c("data.frame", "TransformZscore"), function(df, object){
    vals = df[[object@var]]
    vals = (vals - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)
    df[[object@var]] = vals
    df
})

tst = TransformZscore("year_of_birth")
as.character(tst)
x = tst
object = tst
applyTransform(df, tst)
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
setMethod("applyTransform", c("data.frame", "TransformNegate"), function(df, object){
    vals = df[[object@var]]
    vals = (vals - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)
    df[[object@var]] = vals
    df
})

tst = TransformNegate("year_of_birth")
as.character(tst)
x = tst
object = tst
applyTransform(df, tst)
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
setMethod("applyTransform", c("data.frame", "TransformLog"), function(df, object){
    vals = df[[object@var]]
    vals = log(vals, base = object@power)
    df[[object@var]] = vals
    df
})

tst = TransformLog("year_of_birth")
as.character(tst)
x = tst
object = tst
applyTransform(df, tst)
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
setMethod("applyTransform", c("data.frame", "TransformCrop"), function(df, object){
    vals = df[[object@var]]
    vals[vals < object@min_val] = object@min_val
    vals[vals > object@max_val] = object@max_val
    df[[object@var]] = vals
    df
})

tst = TransformCrop("year_of_birth")
as.character(tst)
x = tst
object = tst
applyTransform(df, tst)
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
setMethod("applyTransform", c("data.frame", "TransformRescale"), function(df, object){
    vals = df[[object@var]]
    vals = scales::rescale(vals, c(object@min_val, object@max_val))
    df[[object@var]] = vals
    df
})

tst = TransformRescale("year_of_birth")
as.character(tst)
x = tst
object = tst
applyTransform(df, tst)
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
setMethod("applyTransform", c("data.frame", "TransformBin"), function(df, object){
    vals = df[[object@var]]
    vals = cut(vals, breaks = object@n_bins)
    df[[object@var]] = vals
    df
})

tst = TransformBin("year_of_birth")
as.character(tst)
x = tst
object = tst
applyTransform(df, tst) %>% applyFilter(FilterCharacter("year_of_birth", one_of = c("(1.95e+03,1.97e+03]")))
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
setMethod("applyTransform", c("data.frame", "TransformRank"), function(df, object){
    vals = df[[object@var]]
    if(object@decreasing){
        vals = frank(-vals, ties.method = "first")
    }else{
        vals = frank(vals, ties.method = "first")
    }
    df[[object@var]] = vals
    df
})

tst = TransformRank("year_of_birth", decreasing = TRUE)
as.character(tst)
x = tst
object = tst
x = df$year_of_birth
y = applyTransform(df, tst)$year_of_birth
y
plot(x,y, xlab = "year", ylab = "rank")
title(as.character(tst))
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
setMethod("applyTransform", c("data.frame", "TransformSort"), function(df, object){
    vals = df[[object@var]]
    o = order(vals, decreasing = object@decreasing)
    df = df[o,]
    df
})

tst = TransformSort("year_of_birth", decreasing = TRUE)
as.character(tst)
x = tst
object = tst
head(applyTransform(df, tst), n = 20)
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
setMethod("applyTransform", c("data.frame", "TransformCharacterRemap"), function(df, object){
    
    vals = df[[object@var]]
    missed = setdiff(vals, names(object@mapping))
    names(missed) = missed
    dict = c(missed, object@mapping)
    vals = dict[vals]
    names(vals) = NULL
    df[[object@var]] = vals
    df
})

tst = TransformCharacterRemap("race", mapping = c("white" = "white or asian", "asian" = "white or asian"))
as.character(tst)
x = tst
object = tst
applyTransform(df, tst)
table(applyTransform(df, tst)$race)
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
setMethod("applyTransform", c("data.frame", "TransformCopy"), function(df, object){
    vals = df[[object@var]]
    df[[object@new_var]] = vals
    df
})

tst = TransformCopy("year_of_birth", "year_of_birth_ranked")
as.character(tst)
x = tst
object = tst
df2 = applyTransform(df, tst) %>%
    applyTransform(TransformRank("year_of_birth_ranked"))
plot(df2$year_of_birth, df2$year_of_birth_ranked)
