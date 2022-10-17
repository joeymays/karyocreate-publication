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

#twoArmPalette adapted from Blake R Mills's MetBrewer package, Cross palette.
#https://github.com/BlakeRMills/MetBrewer
twoArmPalette <- c("#FFB4AD", 
                   "#CEA1CA",
                   "#53939b",
                   "#004d64",
                   "#081f51",
                   "#ffba00",
                   "#fd6f00",
                   "#e12f34",
                   "#7e9a65")

twoArmPaletteOrder <- c("loss_gain", "gain_loss", "neutral_loss", "loss_neutral", "loss_loss",
                        "neutral_gain", "gain_neutral", "gain_gain", "neutral_neutral")
