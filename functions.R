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

# ex_files = dir("example_data", full.names = TRUE)
# # funs = lapply(ex_files, decide_parse_FUN)
# lapply(ex_files, function(f){
#     message(f)
#     decide_parse_FUN(f)
# })
# sapply(paste0(ex_files, "bad"), decide_parse_FUN)
