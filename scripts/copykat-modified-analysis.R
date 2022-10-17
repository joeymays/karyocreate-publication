#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description = "Run Modified Copykat to evaluate CNVs")
parser$add_argument('--experiment.folder', type = "character", required = TRUE, help = "Run folder containing preprocessed Seurat object")
parser$add_argument('--output.prefix', type = "character", required = TRUE, help = "Prefix for output folder in experiment folder and output files")
parser$add_argument('--preprocessed.seurat.object', type = "character", required = TRUE, help = "QC'd and filtered Seurat Object, in experiment folder")
parser$add_argument('--gene.annotation.matrix', type = "character", required = TRUE, help = "Gene annotation matrix")
parser$add_argument('--control.samples', type = "character", required = TRUE, help = "Control sample identities (from Seurat Object)")
parser$add_argument('--evaluation.samples', type = "character", required = TRUE, help = "Sample identities to evaluate, comma delimited (Seurat object)")
parser$add_argument('--MAD.scalar', type = "double", required = FALSE, default = 2.5, help = "MAD scalar for aneuploidy call threshold. Default 2.5")

args <- parser$parse_args()

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
exp.folder <- args$experiment.folder
output.prefix <- args$output.prefix
figure.folder <- paste0(output.prefix, "-copykat-results")

sc.preprocessed.filename <- args$preprocessed.seurat.object
gene.annotation.filename <- args$gene.annotation.matrix

control.htos <- trimws(strsplit(args$control.samples, ',')[[1]])
evaluate.htos <- trimws(strsplit(args$evaluation.samples, ',')[[1]])

MAD.scalar <- args$MAD.scalar

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(copykat)
library(tidyr)
library(SeuratObject)
library(stringr)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#set filter to prefer dplyr::filter
filter <- dplyr::filter

#set seed
set.seed(123)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#experiment & figure folders
suppressWarnings(dir.create(file.path(exp.folder, figure.folder)))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#load custom copykat functions
source(file.path("copykat-scripts","CopykatWithSpecificBreaks.R"))

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
sc <- readRDS(file.path(exp.folder, sc.preprocessed.filename))
sc

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
DefaultAssay(sc) <- "RNA"
sc.raw <- as.matrix(sc@assays$RNA@counts)

print(sc@meta.data %>% group_by(hash.ID) %>% tally(), n = 12)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#generate cell lookup table
cell.lookup <- sc@meta.data %>% select(hash.ID)
cell.lookup$cell.bc <- rownames(cell.lookup)
saveRDS(file = file.path(exp.folder, paste0(output.prefix, "-cell-bc-lookup.RDS")), object = cell.lookup)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(sc) <- "hash.ID"
bcs.evaluate <- WhichCells(sc, idents = evaluate.htos)

cells.eval.raw <- sc.raw[,bcs.evaluate]


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
control.bcs <- WhichCells(sc, idents = control.htos)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
gene.annotations <- readRDS(file.path(gene.annotation.filename))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
copykatNormResult <- copykatSimpleNormalize(rawmat = cells.eval.raw,
                                            norm.cell.names = control.bcs,
                                            prefetched.gene.annotations = gene.annotations,
                                            n.cores = parallel::detectCores(), use.smoothing = T)

annotated.matrix <- copykatNormResult
relative.normal.exp.matrix <- copykatNormResult[,8:ncol(copykatNormResult)]


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
rm(cells.eval.raw)
rm(sc)
rm(sc.raw)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
chr.arms <- annotated.matrix$arm
summary(chr.arms)
breakpoint.ranges <- generateBreakpointRanges(chr.arms)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#combine 18 arms
breakpoint.ranges["18p", "stop"] <- breakpoint.ranges["18q", "stop"]
breakpoint.ranges <- breakpoint.ranges[-(which(breakpoint.ranges$chr == "18q")),]
rownames(breakpoint.ranges)[which(breakpoint.ranges$chr == "18p")] <- "18"
breakpoint.ranges$chr[which(breakpoint.ranges$chr == "18p")] <- "18"

levels(chr.arms)[levels(chr.arms)=="18p" | levels(chr.arms)=="18q"] <- "18"


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
copykat.full.CNA.annotated <- copykatWithSpecificBreakpoints(norm.mat.relat = relative.normal.exp.matrix, anno.mat = annotated.matrix,
                              breakpoint.ranges = breakpoint.ranges, write.CNA.table = F, return.results = T, n.cores = parallel::detectCores())

#combine 18
levels(copykat.full.CNA.annotated$arm)[which(levels(copykat.full.CNA.annotated$arm) %in% c("18p", "18q"))] <- 18

write.table(x = copykat.full.CNA.annotated, file = file.path(exp.folder, figure.folder, paste0(output.prefix, ".arms.CNA.matrix.txt")), sep = "\t", quote = F, row.names = F, col.names = T)

saveRDS(object = copykat.full.CNA.annotated, file = file.path(exp.folder, figure.folder, paste0(output.prefix, ".arms.CNA.matrix.RDS")))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
copykat.full.CNA <- copykat.full.CNA.annotated[,8:ncol(copykat.full.CNA.annotated)]


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
summary(annotated.matrix$chromosome_name)
wholechr.breakpoint.ranges <- generateBreakpointRanges(annotated.matrix$chromosome_name)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
copykat.full.CNA.wholechr.annotated <- copykatWithSpecificBreakpoints(norm.mat.relat = relative.normal.exp.matrix, anno.mat = annotated.matrix,
                              breakpoint.ranges = wholechr.breakpoint.ranges, write.CNA.table = F, return.results = T, n.cores = parallel::detectCores())

write.table(x = copykat.full.CNA.wholechr.annotated, file = file.path(exp.folder, figure.folder, paste0(output.prefix, ".wholechr.CNA.matrix.txt")), sep = "\t", quote = F, row.names = F, col.names = T)

saveRDS(object = copykat.full.CNA.wholechr.annotated, file = file.path(exp.folder, figure.folder, paste0(output.prefix, ".wholechr.CNA.matrix.RDS")))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
copykat.full.CNA.wholechr <- copykat.full.CNA.wholechr.annotated[,8:ncol(copykat.full.CNA.wholechr.annotated)]


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
CNA.means <- data.frame(cell.bc = colnames(copykat.full.CNA))
for(i in unique(copykat.full.CNA.annotated$arm)){
    CNA.means[,paste0("chr",i)] <- copykat.full.CNA.annotated %>% filter(arm == i) %>% dplyr::select(colnames(copykat.full.CNA)) %>% apply(MARGIN = 2, function(x) x[1])
}
CNA.means <- merge(CNA.means, cell.lookup, by = "cell.bc", all.x = T, all.y = F, sort = F) #add sample tag (hash.ID)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
CNA.means.tidy <- reshape2::melt(CNA.means, id.vars = c("cell.bc", "hash.ID"))
colnames(CNA.means.tidy)[3:4] <- c("chr", "cna.value")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
thresholdCalculation <- function(x){x.median <-  median(x); x.mad <- mad(x); return(c(x.median + (x.mad * MAD.scalar), x.median - (x.mad * MAD.scalar)))}


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#function for calling gains or losses based on the thresholds
thresholdCalling <- function(cna.tidy, CNA.thresholds){
  results.df <- data.frame()
  for(current.chr in unique(cna.tidy$chr)){
    cna.tidy.current <- cna.tidy %>% filter(chr == current.chr)
    upper.threshold <- as.numeric(CNA.thresholds %>% filter(chr == current.chr, position == "upper.threshold") %>% select(threshold))
    lower.threshold <- as.numeric(CNA.thresholds %>% filter(chr == current.chr, position == "lower.threshold") %>% select(threshold))
    cna.tidy.current$threshold.call <- sapply(cna.tidy.current$cna.value, FUN = function(x){
      if(x > upper.threshold){return("gain")
      } else if(x < lower.threshold){return("loss")
      } else {return("neutral")}})
    results.df <- rbind(results.df, cna.tidy.current)
  }
  return(results.df)
}


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
CNA.thresholds <- data.frame(CNA.means %>% filter(hash.ID == control.htos) %>% select(-matches(c("cell.bc", "hash.ID"))) %>% apply(MARGIN = 2, FUN = thresholdCalculation))
CNA.thresholds$position <- c("upper.threshold", "lower.threshold")
CNA.thresholds <- reshape2::melt(CNA.thresholds, id.vars = "position")
colnames(CNA.thresholds) <- c("position", "chr", "threshold")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#call gains or losses
CNA.means.tidy <- thresholdCalling(CNA.means.tidy, CNA.thresholds)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
# all chromosomes quantified
CNA.output <- CNA.means.tidy %>% filter(hash.ID %in% evaluate.htos) %>% group_by(hash.ID, chr) %>% count(threshold.call) %>% mutate(percent = round(100 * n/sum(n),1)) %>% select(-n) %>% pivot_wider(names_from = threshold.call, values_from = percent, values_fill = 0) %>% arrange(chr, hash.ID)

write.table(x = CNA.output, file = file.path(exp.folder, figure.folder, paste0(output.prefix, "-CNA-calls-pct-arms.txt")), quote = F, row.names = F, col.names = T, sep = '\t')

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#transpose output
CNA.output.tx.1 <- CNA.output %>% pivot_longer(cols = c(gain, loss, neutral),
                                                     names_to = "class", values_to = "value")
CNA.output.tx.2 <- CNA.output.tx.1 %>%
    pivot_wider(names_from = c(class, hash.ID), values_from = value, names_sort = T) %>%
    tibble::column_to_rownames(var = "chr") %>% t() %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "hash.ID")

write.table(x = CNA.output.tx.2, file = file.path(exp.folder, figure.folder, paste0(output.prefix, "-CNA-calls-pct-arms-transpose.txt")), quote = F, row.names = F, col.names = T, sep = '\t')

## ------------------------------
# determine frequency of both arms being gained or lost from arm data
both.arms.output <- CNA.means.tidy
both.arms.output$full.chr <- paste0("chr", both.arms.output$chr %>% str_extract(pattern = "\\d{1,2}"))

both.arms.output <- both.arms.output %>% group_by(cell.bc, full.chr, hash.ID) %>% summarize(both.arm.gain = all(threshold.call == "gain"), both.arm.loss = all(threshold.call == "loss")) %>% group_by(hash.ID, full.chr) %>% pivot_longer(cols = starts_with("both.arm."), names_to = "class", values_to = "class.call") %>% group_by(hash.ID, full.chr, class) %>% count(class.call) %>% mutate(percent = round(100 * n/sum(n),1)) %>% select(-n) %>% filter(class.call == TRUE) %>% select(-class.call) %>% pivot_wider(names_from = c(class, hash.ID), values_from = percent, names_sort = T, values_fill = 0) %>% tibble::column_to_rownames(var = "full.chr") %>% t() %>% as.data.frame()

both.arms.output <- both.arms.output[,gtools::mixedorder(colnames(both.arms.output))]

both.arms.output <- both.arms.output %>% tibble::rownames_to_column(var = "hash.ID")

write.table(x = both.arms.output, file = file.path(exp.folder, figure.folder, paste0(output.prefix, "-CNA-calls-pct-both-arms.txt")), quote = F, row.names = F, col.names = T, sep = '\t')

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
CNA.means <- data.frame(cell.bc = colnames(copykat.full.CNA.wholechr))
for(i in unique(copykat.full.CNA.wholechr.annotated$chromosome_name)){
    CNA.means[,paste0("chr",i)] <- copykat.full.CNA.wholechr.annotated %>% filter(chromosome_name == i) %>% dplyr::select(colnames(copykat.full.CNA.wholechr)) %>% apply(MARGIN = 2, function(x) x[1])
}
CNA.means <- merge(CNA.means, cell.lookup, by = "cell.bc", all.x = T, all.y = F, sort = F) #add sample tag (hash.ID)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
CNA.means.tidy <- reshape2::melt(CNA.means, id.vars = c("cell.bc", "hash.ID"))
colnames(CNA.means.tidy)[3:4] <- c("chr", "cna.value")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
CNA.thresholds <- data.frame(CNA.means %>% filter(hash.ID == control.htos) %>% select(-matches(c("cell.bc", "hash.ID"))) %>% apply(MARGIN = 2, FUN = thresholdCalculation))
CNA.thresholds$position <- c("upper.threshold", "lower.threshold")
CNA.thresholds <- reshape2::melt(CNA.thresholds, id.vars = "position")
colnames(CNA.thresholds) <- c("position", "chr", "threshold")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#call gains or losses
CNA.means.tidy <- thresholdCalling(CNA.means.tidy, CNA.thresholds)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------
# all chromosomes quantified
CNA.output <- CNA.means.tidy %>% filter(hash.ID %in% evaluate.htos) %>% group_by(hash.ID, chr) %>% count(threshold.call) %>% mutate(percent = round(100 * n/sum(n),1)) %>% select(-n) %>% pivot_wider(names_from = threshold.call, values_from = percent, values_fill = 0) %>% arrange(chr, hash.ID)

write.table(x = CNA.output, file = file.path(exp.folder, figure.folder, paste0(output.prefix, "-CNA-calls-pct-wholechr.txt")), quote = F, row.names = F, col.names = T, sep = '\t')

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
#transpose output
CNA.output.tx.1 <- CNA.output %>% pivot_longer(cols = c(gain, loss, neutral),
                                               names_to = "class", values_to = "value")
CNA.output.tx.2 <- CNA.output.tx.1 %>%
    pivot_wider(names_from = c(class, hash.ID), values_from = value, names_sort = T) %>%
    tibble::column_to_rownames(var = "chr") %>% t() %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "hash.ID")

write.table(x = CNA.output.tx.2, file = file.path(exp.folder, figure.folder, paste0(output.prefix, "-CNA-calls-pct-wholechr-transpose.txt")), quote = F, row.names = F, col.names = T, sep = '\t')

## ----------------------------------------------------------------------------------------------------------------------------------------------------------
sessioninfo::session_info()
