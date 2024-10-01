# Set up work env 
setwd("/home/sguizard/Work/Projects/Intronless/CAGEseq")
.libPaths(c(.libPaths(), "/home/sguizard/miniforge3/envs/cagefightr/lib/R/library"))
`%nin%` <- Negate(`%in%`)

# Loading libraries
library(tidyverse)
library(CAGEfightR)
library(GenomeInfoDb)    # SeqInfo
library(rtracklayer)     # Export
library(GenomicFeatures) # makeTxDbFromGFF makeTxDbFromEnsembl
source("https://raw.githubusercontent.com/sguizard/bioinfo-scripts/refs/heads/master/utils/readGxf.R")

# Load rds files
BCs <- readRDS("BCs_min1sample_min10reads.rds")
TCs <- readRDS("TCs_min1sample_min10reads.rds")

# Extract SEG and write gtf file
extractSEG(
  file     = "chicken_mat_ensembl_gs.gtf", 
  out_file = "chicken_mat_ensembl_gs_SEG.gtf"
)

seg_gtf <- read_gtf("chicken_mat_ensembl_gs_SEG.gtf", separate_attributes = TRUE)

txdb <- makeTxDbFromGFF(
  file     = "chicken_mat_ensembl_gs_SEG.gtf",
  format   = "gtf",
  organism = "Gallus gallus")

TCs$txID   <- NULL
TCs$txType <- NULL

TCs <- assignTxID  (TCs, txModels = txdb)
TCs <- assignTxType(TCs, txModels = txdb)

# ### Overlap summary: ###
# txType count percentage
# 1   promoter   896        2.7
# 2   proximal   856        2.6
# 3    fiveUTR     0        0.0
# 4   threeUTR     0        0.0
# 5        CDS     0        0.0
# 6       exon   681        2.1
# 7     intron     0        0.0
# 8  antisense   627        1.9
# 9 intergenic 29787       90.7

TCs_seg <- TCs[TCs$txType %in% c("promoter", "proximal", "exon")]
export.bed(TCs_seg, "TCs_min1sample_min10reads_SEG.bed")

TCs_seg_ids <- 
  as_tibble(TCs_seg) %>%
  distinct(txID) %>%
  separate_longer_delim(cols = txID, delim = ";") %>%
  mutate(txID = 
    str_replace_all(txID, '"', '') %>%
    str_replace_all("\\.\\d+$", "")) %>% 
  pull(txID)

seg_gtf %>% 
  filter(gene_id %in% TCs_seg_ids) %>% 
  dplyr::select(sequence:gene_id, transcript_id) %>% 
  mutate(attributes = paste0('gene_id="', gene_id, '";transcript_id="', transcript_id, '"')) %>% 
  dplyr::select(-gene_id, -transcript_id) %>% 
  write_tsv("chicken_mat_ensembl_gs_SEG_with_CAGE_peak.gtf")

BCs$txID   <- NULL
BCs$txType <- NULL

BCs <- assignTxID  (BCs, txModels = txdb)
BCs <- assignTxType(BCs, txModels = txdb)

# ### Overlap summary: ###
# txType count percentage
# 1   promoter   279       12.0
# 2   proximal   145        6.2
# 3    fiveUTR     0        0.0
# 4   threeUTR     0        0.0
# 5        CDS     0        0.0
# 6       exon    57        2.4
# 7     intron     0        0.0
# 8  antisense     0        0.0
# 9 intergenic  1849       79.4

