convert.barcode <- function(cb.vector){
    if(length(grep("\\.", x = cb.vector)) != 0){ #if periods, change to dashes
        new.cb <- gsub(x = cb.vector, pattern = "\\.", replacement = "\\-")
        return(new.cb)
    } else if(length(grep("\\-", x = cb.vector)) != 0){ #if dashes, change to periods
        new.cb <- gsub(x = cb.vector, pattern = "\\-", replacement = "\\.")
        return(new.cb)
    } else { #if neither, return original string
        return(cb.vector)   
    }
}

corner <- function(input.mat){
    if(nrow(input.mat) > 4){
        row.out = 5
    } else { row.out = nrow(input.mat)}
    
    if(ncol(input.mat) > 7){
        col.out = 8
    } else { col.out = ncol(input.mat)}
    
    return(input.mat[1:row.out, 1:col.out])
}
