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
chrNameLength <-
    read_tsv(
        opt$chr_length,
        col_names = c("seqnames", "seqlenghts")) %>%
    mutate(
        isCicular = if_else(str_detect(seqnames, "^M"), TRUE, FALSE),
        genome = opt$genome)

ov_ars <- Seqinfo(
    seqnames   = chrNameLength$seqnames,
    seqlengths = chrNameLength$seqlenghts,
    isCircular = chrNameLength$isCicular,
    genome     = chrNameLength$genome)

# List BigWig files and attach ids
bw <-
    tibble(
        bw_plus  = list.files(path = opt$bw_dir, pattern = "*Unique.str1.out.bw"),
        bw_minus = list.files(path = opt$bw_dir, pattern = "*Unique.str2.out.bw"),
        run_id   = str_extract(bw_plus, "CAGE_(.+)_R[12]"))

# Prepare design object
design <-
    bw %>%
    dplyr::select(
        Name        = run_id,
        BigWigPlus  = bw_plus,
        BigWigMinus = bw_minus) %>%
    as.data.frame()

rownames(design) <- design$Name

# Loading BigWig
bw_plus  <- BigWigFileList(bw$bw_plus)
bw_minus <- BigWigFileList(bw$bw_minus)

names(bw_plus) <- names(bw_minus) <- bw$sample_id

# Quantify the number of CAGE-tags for each CTSSs from BigWig-files
CTSSs <- quantifyCTSSs(
    plusStrand  = bw_plus,
    minusStrand = bw_minus,
    genome      = ov_ars,
    design      = design)

# Unidirectionnal clustering (TSSs)
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

export.bed(TCs, con = "TCs_min1sample_min10reads.bed")


# Birectionnal clustering (enhancers)
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

export.bed(object = BCs, con = "BCs_min1sample_min10reads.bed")

# Assign transcripts (TSS and enhancers)
txdb <- makeTxDbFromGFF(
    file     = opt$gtf,
    format   = "gtf",
    organism = opt$organism)

TCs <- assignTxID(TCs, txModels = txdb)
TCs <- assignTxType(TCs, txModels = txdb)

BCs_all <- assignTxID(BCs, txModels = txdb)
BCs_all <- assignTxType(BCs_all, txModels = txdb)
