#!/usr/bin/env Rscript

library(getopt)

`%nin%` <- Negate(`%in%`)

spec <- matrix(c(
    "help"      , "h", 0 , "logical",
    "chr_length", "c", "", "character",
    "bw_dir"    , "b", "", "character",
    "genome"    , "g", "", "character",
    "organism"  , "o", "", "character",
    "annot"     , "a", "", "character"
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

# CLI tests
# opt$chr_length <- "chrNameLength.txt"
# opt$bw_dir <- "."
# opt$genome <- "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
# opt$annot <- "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf"
# opt$organism <- "Gallus gallus"

if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}
if (is.null(opt$chr_length)) {
    cat("MISSING OPTION: --chr_length \n")
    q(status = 1)
}
if (is.null(opt$bw_dir)) {
    cat("MISSING OPTION: --bw_dir \n")
    q(status = 1)
}
if (is.null(opt$genome)) {
    cat("MISSING OPTION: --genome \n")
    q(status = 1)
}
if (is.null(opt$organism)) {
    cat("MISSING OPTION: --organism \n")
    q(status = 1)
}
if (is.null(opt$annot)) {
    cat("MISSING OPTION: --annot \n")
    q(status = 1)
}

# Loading libraries
library(tidyverse)
library(CAGEfightR)
library(GenomeInfoDb)    # SeqInfo
library(rtracklayer)     # Export
library(GenomicFeatures) # makeTxDbFromGFF makeTxDbFromEnsembl

# Load chromosomes lengths and prepare columns for Seqinfo
cat(paste0("==> Load ", opt$chr_length, "\n"))
chrNameLength <-
    read_tsv(
        opt$chr_length,
        col_names = c("seqnames", "seqlengths")) %>%
    mutate(
        isCicular = if_else(str_detect(seqnames, "^M"), TRUE, FALSE),
        genome = opt$genome)

cat(paste0("==> Prepare Seqinfo object", "\n"))
seqinfo <- Seqinfo(
    seqnames   = chrNameLength$seqnames,
    seqlengths = chrNameLength$seqlengths,
    isCircular = chrNameLength$isCicular,
    genome     = chrNameLength$genome)

# List BigWig files and attach ids
cat(paste0("==> List bw files in ", opt$bw_dir, "\n"))
bw <-
    tibble(
        bw_plus  = list.files(path = opt$bw_dir, pattern = "*Unique.str1.out.bw"),
        bw_minus = list.files(path = opt$bw_dir, pattern = "*Unique.str2.out.bw"),
        run_id   = str_extract(bw_plus, "CAGE_(.+)_R[12]"))


processBigWig <- function(idx) {
    cat(paste0("####### Working on ", bw[idx,]$run_id, " #######\n"))
    # Prepare design object
    cat(paste0("==> Prepare design object", "\n"))
    design <-
        bw[idx,] %>%
        dplyr::select(
            Name        = run_id,
            BigWigPlus  = bw_plus,
            BigWigMinus = bw_minus) %>%
        as.data.frame()

    rownames(design) <- design$Name

    # Loading BigWig
    cat(paste0("==> Load bw files", "\n"))
    bw_plus  <- BigWigFileList(bw$bw_plus[idx])
    bw_minus <- BigWigFileList(bw$bw_minus[idx])

    names(bw_plus) <- names(bw_minus) <- bw$run_id[idx]

    # Quantify the number of CAGE-tags for each CTSSs from BigWig-files
    cat(paste0("==> Quantify CTSSs", "\n"))
    CTSSs <- quantifyCTSSs(
        plusStrand  = bw_plus,
        minusStrand = bw_minus,
        genome      = seqinfo,
        design      = design)

    # Unidirectionnal clustering (TSSs)
    cat(paste0("==> Unidirectional clustering", "\n"))
    TCs <-
        CTSSs %>%
        calcTPM(inputAssay = "counts", outputAssay = "TPM") %>% # Compute Tag Per Million
        calcPooled(inputAssay = "TPM") %>%                      # Compute TPM over all samples
        subsetBySupport(                                        # Filter CTSS shared by min 1 sample with at least 10 reads support
            inputAssay   = "counts",
            outputColumn = "support",
            unexpressed  = 10,
            minSamples   = 0) %>%
        clusterUnidirectionally(
            pooledCutoff = 2,
            mergeDist = 20)
    
    cat(paste0("===> Export ", bw$run_id[idx], "_TCs_min1sample_min10reads.bed", "\n"))
    export.bed(TCs, con = paste0(bw$run_id[idx], "_TCs_min1sample_min10reads.bed"))


    # Birectionnal clustering (enhancers)
    cat(paste0("==> Bidirectionnal clustering", "\n"))
    BCs <-
        CTSSs %>%
        calcTPM(inputAssay = "counts", outputAssay = "TPM") %>% # Compute Tag Per Million
        calcPooled(inputAssay = "TPM") %>%                      # Compute TPM over all samples
        subsetBySupport(                                        # Filter CTSS shared by min 1 sample with at least 10 reads support
            inputAssay   = "counts",
            outputColumn = "support",
            unexpressed  = 10,
            minSamples   = 0) %>%
        clusterBidirectionally(balanceThreshold = 0.95)

    cat(paste0("====> Export ", bw$run_id[idx], "_BCs_min1sample_min10reads.bed", "\n"))
    export.bed(object = BCs, con = paste0(bw$run_id[idx], "_BCs_min1sample_min10reads.bed"))

    # Assign transcripts (TSS and enhancers)
    cat(paste0("==> Assign transcripts", "\n"))
    txdb <- makeTxDbFromGFF(
        file     = opt$annot,
        format   = "gtf",
        organism = opt$organism)

    TCs <- assignTxID  (TCs, txModels = txdb)
    TCs <- assignTxType(TCs, txModels = txdb)
    cat(paste0("====> Export ", bw$run_id[idx], "_TCs_min1sample_min10reads.rds", "\n"))
    saveRDS(object = TCs, file = paste0(bw$run_id[idx], "_TCs_min1sample_min10reads.rds"))

    BCs_all <- assignTxID  (BCs,     txModels = txdb)
    BCs_all <- assignTxType(BCs_all, txModels = txdb)
    cat(paste0("====> Export ", bw$run_id[idx], "_BCs_min1sample_min10reads.rds", "\n"))
    saveRDS(object = BCs_all, file = paste0(bw$run_id[idx], "_BCs_min1sample_min10reads.rds"))
}

map(seq_len(length.out = nrow(bw)), processBigWig)
