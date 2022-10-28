
# Modified version of copykat algorithm that outputs smoothed CNA levels across supplied breakpoints.
# Joey Mays 2021-11-03


### Normalize Function, Modified from the original, excludes some steps originally used for breakpoint discovery, includes more options for normalization.
#rawmat: raw gene count matrix
#LOW.DR: excludes genes not expressed in a minimum % of cells for normalization, smoothing, and baseline adjustment
#UP.DR: excludes genes not expressed in a minimum % of cells for output and subsequent CNA analysis
#norm.cell.names: cell barcodes for control cells; if not supplied, the whole dataset is used to find the baseline
#n.cores: specify number of cores to use; used for smoothing, and baseline generation without specified norm.cell.names
#report.cytoband: FALSE (default), reports chromosome arms in annotation; TRUE, reports cytoband
#use.smoothing: toggles DLM smoothing step
#prefetched.gene.annotations: annotations retrieved using getGeneMetadata().
#use.package.annotation: superceeded by setting prefetched.gene.annotations;
    #TRUE will use gene annotations shipped with original copykat package;
    #FALSE will grab annotations on-the-fly using biomaRt
#mirror: mirror for biomaRt if use.package.annotation = F; "uswest", "useast"
#normalization.method: specifies normalization method; "FT" (default), Freeman-Tukey; "shiftedLog" y = log(x + pseudocount), etc.
#pseudocount: pseudocount for normalization if applicable

## original copykat parameters, generally ignored or obsolete for modified analysis
#id.type:"S" for gene Symbol, uses gene symbol if using copykat package gene annotation; "E"for Ensembl ID
#cell.line: used for baseline generation if no baseline control supplied
#ngene.chr: sets minimum number of genes per chromosome; obsolete
#sam.name: used for output file names; obsolete
#distance: used for breakpoint finding; obsolete
#output.seg; obsolete

copykatSimpleNormalize <- function(rawmat=rawdata, id.type="S", cell.line="no",
                                      ngene.chr=5,LOW.DR=0.05, UP.DR=0.1, norm.cell.names="",
                                    sam.name="", distance="euclidean", output.seg="FALSE", n.cores=1,
                                   report.cytoband = F, use.package.annotation = T, prefetched.gene.annotations = NULL,
                                   use.smoothing = T, mirror = "uswest", normalization.method = "FT",
                                   pseudocount = 1){

    #load copykat library
    if(!exists("copykat")){
        library(copykat)
    }

    if(norm.cell.names[1] == ""){
        message("WARNING: No control cells specified.")
    }

    ### SETUP

    set.seed(123)
    #sample.name <- paste(sam.name,"_copykat_", sep="")

    print("running custom script, modified from copykat v1.0.5 updated 07/15/2021")
    print("step1: read and filter data ...")
    print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", sep=""))

    ### FILTERING 1

    # cells with less than 200 genes are filtered out, only prints if cells are filtered.
    # get number of genes with greater than 0 exp for each cells (i.e. no. of genes exp in each cell)
    genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))
    if(sum(genes.raw<100)>1){
        rawmat <- rawmat[, -which(genes.raw< 200)]
        print(paste("filtered out ", sum(genes.raw<=200), " cells with less than 200 genes; remaining ", ncol(rawmat), " cells", sep=""))
    }

    # LOW.DR: Keeps genes that are expressed in at least <LOW.DR> percent of cells (5% Default)
    #for each gene, find proportion of cells the gene is expressed in.
    der <- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)

    #if there are genes that are expressed in over 5% (LOW.DR default) of cells
    if(sum(der>LOW.DR)>=1){
        rawmat <- rawmat[which(der > LOW.DR), ] #keep the genes that are expressed in at least 5% (LOW.DR default) of cells
        print(paste(nrow(rawmat)," genes past LOW.DR filtering", sep=""))
    }

    WNS1 <- "data quality is ok"
    if(nrow(rawmat) < 7000){ #data is considered low quality if number of genes is less than 7000; LOW.DR is then set to UP.DR
        WNS1 <- "low data quality"
        UP.DR<- LOW.DR
        print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
    }

    ### GENE ANNOTATIONS
    #Supplied biomaRt gene annotations are used if provided (recommended). If not, user can set use.package.annotation = T to use the annotation supplied with copykat package.
    #If running copykat locally, the gene annotations can be imported from biomaRt on the fly by specifying use.package.annotation = F & prefetched.gene.annotations = NULL.
    #Importing from biomaRt on the HPC has previously thrown errors, probably due to some authentication issues.

    if(is.null(prefetched.gene.annotations)){
        if(!use.package.annotation == F){
            print("Using local annotation...")
            anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type)
            anno.mat <- anno.mat[order(anno.mat$abspos, decreasing = FALSE),]
        } else {
            if(!exists("getGeneMetadata")){
                stop("getGeneMetadata.R script not loaded for remote annotation.")
            }
            print("Using remote annotation...")
            gene.annotations <- getGeneMetadata(gene.list = rownames(rawmat), name.type = "symbol", sex.chr = c("X","Y"), mirror = mirror)
            gene.annotations <- gene.annotations[,c("chromosome_name","start_position","end_position",
                                                    "ensembl_gene_id","hgnc_symbol","band")]

           #change X,Y to 23,24 for easier processing
            gene.annotations[which(gene.annotations$chromosome_name == "X"), "chromosome_name"] <- 23
            gene.annotations[which(gene.annotations$chromosome_name == "Y"), "chromosome_name"] <- 24
            gene.annotations$chromosome_name <- as.numeric(gene.annotations$chromosome_name)

            gene.annotations <- gene.annotations[order(gene.annotations$chromosome_name, gene.annotations$start_position),]
            gene.annotations <- cbind(data.frame(rank = 1:nrow(gene.annotations)), gene.annotations)

            rawmat.premerge <- cbind(data.frame(hgnc_symbol = rownames(rawmat)), rawmat)

            anno.mat <- merge(gene.annotations, rawmat.premerge,
                              by = "hgnc_symbol", all.y = F, all.x = T)
            anno.mat <- anno.mat[order(anno.mat$rank),]
        }
    } else {
        print("Using provided annotation...")
        gene.annotations <- prefetched.gene.annotations[,c("chromosome_name","start_position","end_position",
                                                "ensembl_gene_id","hgnc_symbol","band")]

        #change X,Y to 23,24 for easier processing
        gene.annotations[which(gene.annotations$chromosome_name == "X"), "chromosome_name"] <- 23
        gene.annotations[which(gene.annotations$chromosome_name == "Y"), "chromosome_name"] <- 24
        gene.annotations$chromosome_name <- as.numeric(gene.annotations$chromosome_name)

        gene.annotations <- gene.annotations[order(gene.annotations$chromosome_name, gene.annotations$start_position),]
        gene.annotations <- cbind(data.frame(rank = 1:nrow(gene.annotations)), gene.annotations)

        rawmat.premerge <- cbind(data.frame(hgnc_symbol = rownames(rawmat)), rawmat)

        anno.mat <- merge(gene.annotations, rawmat.premerge,
                          by = "hgnc_symbol", all.y = F, all.x = F)
        anno.mat <- anno.mat[order(anno.mat$rank),]
    }


    # HLA and Cell cycle genes are removed
    HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
    toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))
    if(length(toRev)>0){
        anno.mat <- anno.mat[-toRev, ]
    }

    ### NORMALIZATION
    #Freeman-Tukey normalization is used as default, I've added a few more options to try.

    rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
    rownames(rawmat3) <- anno.mat$hgnc_symbol

    if(normalization.method == "FT"){
        print("Normalization: Freeman-Tukey")
        #freeman-tukey normalization
        norm.mat<- log(sqrt(rawmat3)+sqrt(rawmat3+1))
        norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x))) #centering
        colnames(norm.mat) <-  colnames(rawmat3)
    } else if(normalization.method == "shiftedLog"){
        print(paste("Normalization: Shifted-Log, Pseudocount = ", pseudocount))
        norm.mat <- log(rawmat3 + pseudocount)
        norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x))) #centering
        colnames(norm.mat) <-  colnames(rawmat3)
    } else if(normalization.method == "sctransform"){
        print(paste("Normalization: scTransform"))
        norm.mat <- sctransform::vst(rawmat3)$y
        norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x))) #centering
        colnames(norm.mat) <-  colnames(rawmat3)
    } else if(normalization.method == "cp10k"){
        print(paste("Normalization: CP10K"))
        norm.mat <- apply(rawmat3, 2, function(x) x/sum(x))
        norm.mat <- norm.mat * 10000
        norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x))) #centering
        colnames(norm.mat) <-  colnames(rawmat3)
    } else if(normalization.method == "raw"){
        print(paste("Normalization: Raw Counts"))
        norm.mat <- rawmat3
        colnames(norm.mat) <-  colnames(rawmat3)
    } else {
        print("No Normalization Selected.")
    }

    #smoothing data to remove outliers, can be turned on or off
    if(!use.smoothing==F){
        print("step 3: smoothing data with dlm ...")
        dlm.sm <- function(c){
            model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
            x <- dlm::dlmSmooth(norm.mat[, c], model)$s
            x<- x[2:length(x)]
            x <- x-mean(x)
        }

        test.mc <-parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = n.cores)
        norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)
        colnames(norm.mat.smooth) <- colnames(norm.mat)
        rownames(norm.mat.smooth) <- rownames(norm.mat)

        #norm.mat.smooth is outlier-smoothed normalized matrix centered around 0
    } else {
        print("Skipping Step 3: No smoothing used...")
        norm.mat.smooth <- norm.mat
    }


    ### SET CONTROL BASELINE
    #Uses control cell population to set baseline

    print("step 4: measuring baselines ...")
    if (cell.line=="yes"){
        print("running pure cell line mode")
        relt <- baseline.synthetic(norm.mat=norm.mat.smooth, min.cells=10, n.cores=n.cores)
        norm.mat.relat <- relt$expr.relat
        CL <- relt$cl
        WNS <- "run with cell line mode"
        preN <- NULL
    } else if(length(norm.cell.names)>1){ #if basesline (i.e. diploid control) sample is supplied

        NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% norm.cell.names)]) # reports No. of control cells
        print(paste(NNN, " known normal cells found in dataset", sep=""))

        if (NNN==0) stop("known normal cells provided; however none existing in testing dataset") #stops if normal cells aren't in actual data
        print("run with known normal...")

        #basel is median normalized exp level of each gene from control cells
        basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% norm.cell.names)],1,median); print("baseline is from known input")

        WNS <- "run with known normal"
        preN <- norm.cell.names

        #normalized matrix is reduced by the baseline level of expression, norm.mat.relat is normalized matrix relative to control cells
        norm.mat.relat <- norm.mat.smooth-basel

    }else { #if no baseline control supplied
        basa <- baseline.norm.cl(norm.mat.smooth=norm.mat.smooth, min.cells=5, n.cores=n.cores)
        basel <- basa$basel
        WNS <- basa$WNS
        preN <- basa$preN
        CL <- basa$cl
        if (WNS =="unclassified.prediction"){
            basa <- baseline.GMM(CNA.mat=norm.mat.smooth, max.normal=5, mu.cut=0.05, Nfraq.cut=0.99,RE.before=basa,n.cores=n.cores)
            basel <-basa$basel
            WNS <- basa$WNS
            preN <- basa$preN
        }
        norm.mat.relat <- norm.mat.smooth-basel
    }

    #norm.mat.relat is normalized matrix relative to control cells

    ### FILTERING 3
    #filters genes for use in CNA calling, by UP.DR

    #rawmat3 is filtered, non-normalized matrix
    #for each gene in filtered matrix, get proportion of cells with exp of a given gene
    DR2 <- apply(rawmat3,1,function(x)(sum(x>0)))/ncol(rawmat3)

    #norm.mat.relat is now relative expression to baseline control, with genes expressed in under 10% of cells filtered out
    norm.mat.relat <- norm.mat.relat[which(DR2>=UP.DR),]

    #keeps genes with expression in at least 10% (UP.DR default) of cells; anno.mat is same numbers of rawmat3
    anno.mat2 <- anno.mat[which(DR2>=UP.DR), ]
    print(paste(nrow(norm.mat.relat)," genes past UP.DR filtering", sep=""))

    #re-annotate
        if(!use.package.annotation == F & is.null(prefetched.gene.annotations)){
            norm.mat.relative.annotated <- annotateGenes.hg20(mat = norm.mat.relat, ID.type = id.type)
        } else {
            normmat.premerge <- cbind(data.frame(hgnc_symbol = rownames(norm.mat.relat)), norm.mat.relat)
            norm.mat.relative.annotated <- merge(gene.annotations, normmat.premerge,
                                                 by = "hgnc_symbol", all.y = F, all.x = F)
            norm.mat.relative.annotated <- norm.mat.relative.annotated[order(norm.mat.relative.annotated$rank),]
        }


    #factorize chromosome names
    norm.mat.relative.annotated$chromosome_name <- factor(norm.mat.relative.annotated$chromosome_name)

    #replace band column with arms if specified
    if(!report.cytoband == T){
        print("Reporting Chromosome Arms...")
        norm.mat.relative.annotated$band <- paste0(norm.mat.relative.annotated$chromosome_name, substr(norm.mat.relative.annotated$band, 1, 1))
        colnames(norm.mat.relative.annotated)[7] <- "arm"

        norm.mat.relative.annotated$arm <- factor(norm.mat.relative.annotated$arm, levels = gtools::mixedsort(unique(norm.mat.relative.annotated$arm)))
    } else {
        print("Reporting Cytobands...")
    }

    return(norm.mat.relative.annotated)
}

### Segmentation Function: generates mean CNA values for copykat normalized matrix, given specific breakpoint (chromosome/arm) ranges
#norm.mat.relat: normalized matrix from copykatSimpleNormalize()
#anno.mat: annotated matrix from copykatSimpleNormalize()
#breakpoint.ranges: ranges to group genes by chromosome or arm; easily generated by generateBreakpointRanges()
#return.results: return mean CNA matrix as object
#write.CNA.table: write mean CNA matrix to file
#sample.name: file name prefix if write.CNA.table = T
#n.cores: specify number of cores to use for Monte Carlo simulations

copykatWithSpecificBreakpoints <- function(norm.mat.relat, anno.mat, breakpoint.ranges = NULL, write.CNA.table = T, return.results = F,
                                          sample.name = "", n.cores = 1, actual.means = F){

    ### SEGMENTATION

    if(is.null(breakpoint.ranges)){
        stop("No Breakpoints Specified.")
    }

    print("step 5: segmentation...")

    #standard segmentation and results function
    #results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = KS.cut, n.cores=n.cores)

    #modified function, does not generate breakpoints, only uses custom breakpoints to get 1 smoothed CNA estimate (mean) for each region
    BR <- breakpoint.ranges
    fttmat <- norm.mat.relat

    n <- nrow(fttmat)
    norm.mat.sm <- exp(fttmat)

    seg <- function(z) {
        x <- numeric(n)
        for(i in 1:nrow(BR)){
            
            if(actual.means == T){
                x[BR[i, "start"]:BR[i, "stop"]] <- mean(norm.mat.sm[BR[i, "start"]:BR[i, "stop"], z])
            } else {
            
            a <- max(mean(norm.mat.sm[BR[i, "start"]:BR[i, "stop"], z]), 0.001)
            posterior1 <- MCMCpack::MCpoissongamma(norm.mat.sm[BR[i, "start"]:BR[i, "stop"], z], a, 1, mc = 1000)
            x[BR[i, "start"]:BR[i, "stop"]] <- mean(posterior1)
            }
        }
        x <- log(x)
    }

    seg.test <- parallel::mclapply(1:ncol(norm.mat.sm), seg, mc.cores = n.cores) #segment each cell based on breakpoints
    logCNA <- matrix(unlist(seg.test), ncol = ncol(norm.mat.sm), byrow = FALSE)
    colnames(logCNA) <- colnames(norm.mat.sm)

    logCNA <- apply(logCNA, 2, function(x)(x <- x-mean(x))) # for each cell in result matrix, subtract the mean of the cell exp. (shifts center to 0)

    annotated.results <- cbind(anno.mat[, 1:7], logCNA) #add annotations to results

    if(write.CNA.table){
        print("Writing CNA table to file...")
        write.table(annotated.results, paste(sample.name, "CNA_customBreaks_results_gene_by_cell.txt", sep=""), sep="\t", row.names = FALSE, quote = F)
    }

    if(return.results){
        return(annotated.results)
    }
}

### generates a list of breakpoint ranges based on the chromosome vector in the matrix annotation.
#i.e. takes the running length of each vector element and reports the start and stop index of each element run
generateBreakpointRanges <- function(chromosome.vector){
    chr.arm.summary <- summary(chromosome.vector)
    breakpoint.ranges <- data.frame(chr = names(chr.arm.summary), start = NA, stop = NA)
    breakpoint.ranges[1,"start"] <- 1
    breakpoint.ranges[1,"stop"] <- chr.arm.summary[1]

    for(i in 2:length(chr.arm.summary)){
        breakpoint.ranges[i,"start"] <- breakpoint.ranges[(i-1),"stop"] + 1
        breakpoint.ranges[i,"stop"] <- breakpoint.ranges[i, "start"] + (chr.arm.summary[i] - 1)
    }
    rownames(breakpoint.ranges) <- breakpoint.ranges$chr

    return(breakpoint.ranges)
}


###Ranks CNA values for the comprehensive heatmap by specific chromosome for each sample
#to.rank is list where names of items are the samples and items are lists of chromosome arms
#e.g.: to.rank <- list(`HTO-05` = c("7p", "7q"), `HTO-08` = c("7p", "7q"), `HTO-09` = c("18"))
#cell.lookup is data frame of sample IDs and cell barcodes
RankChromosomeHeatmap <- function(copykat.full.CNA.annotated, to.rank, cell.lookup){
    
    master.ordered.bcs <- c()
    
    for(i in 1:length(to.rank)){
        bcs <- cell.lookup %>% filter(hash.ID == names(to.rank)[i])
        bcs <- bcs$cell.bc
        
        if(length(to.rank[[i]]) > 1){
            current.block <- copykat.full.CNA.annotated[match(to.rank[[i]], copykat.full.CNA.annotated$arm),bcs]
            block.ranks <- rbind(rank(current.block[1,]), rank(current.block[2,]))
            block.ranks <- reshape2::melt(block.ranks)
            
            consensus.ranks.setup <- data.frame(Reviewer = as.character(block.ranks$Var1), Item = as.character(block.ranks$Var2), Ranking = block.ranks$value)
            block.consensus.ranks <- RankAggregator::consensusRanking(consensus.ranks.setup)
            block.consensus.ranks.sorted <- block.consensus.ranks[order(block.consensus.ranks[,2]),]
            ordered.bcs <- block.consensus.ranks.sorted$Item
            
        } else if(length(to.rank[[i]]) == 1) {
            current.block <- copykat.full.CNA.annotated %>% filter(arm == to.rank[[i]])
            current.block <- current.block[1,bcs]
            current.block.rank <- rank(current.block)
            current.block.rank <- current.block.rank[order(current.block.rank)]
            ordered.bcs <- names(current.block.rank)
            
        }
        master.ordered.bcs <- c(master.ordered.bcs, ordered.bcs)
    }
    return(master.ordered.bcs)
}
