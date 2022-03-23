library(tidyverse)

read_bed12 <- function(file) {require(readr);read_tsv(file, col_names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"), col_types = cols(chrom = col_character(), chromStart = col_integer(), chromEnd = col_integer(), name = col_character(), score = col_double(), strand = col_character(), thickStart = col_integer(), thickEnd = col_integer(), itemRgb = col_character(), blockCount = col_integer(), blockSizes = col_character(), blockStarts = col_character()))}

read_gff <- function(file) {require(readr); read_tsv(file, col_names = c("sequence", "source", "feature","start","end","score","strand","phase","attributes"), col_types = cols(sequence = "c", source = "c", feature = "c", start = "i", end = "i", score = "c", strand = "c", phase = "c", attributes = "c"), comment = "#")}

ens <- 
    read_bed12('chicken_ensembl.bed') %>%  
    mutate( 
        tama_gene_id = str_extract(name, 'G\\d+'),  
        tama_transcript_id = str_extract(name, 'G\\d+\\.\\d+'), 
        ens_gene = str_extract(name, 'ENSGALG\\d+'),  
        ens_transcript = str_extract(name, 'ENSGALT\\d+'))

id2source <- 
    read_tsv('chicken_ensembl_trans_report.txt') %>% 
    select(tama_transcript_id = transcript_id, sources)

new2old <- 
    read_tsv('chicken_ensembl_trans_report.txt') %>%  
    select(tama_transcript_id = transcript_id, old_id = all_source_trans) %>%  
    separate_rows(old_id, sep = ",") %>%   
    mutate(old_id = str_replace_all(old_id, '(ensembl_|isoseq_)', '')) %>% 
    distinct()

ens <- left_join(ens, id2source)

new2old <- 
    read_tsv('chicken_read_support.txt') %>%  
    select(old_id = merge_trans_id, trans_read_count, source_line) %>%  
    right_join(new2old)

new2old <-  
    read_gff("novel.gtf") %>%  
    filter(feature == "transcript") %>%  
    mutate( 
        old_id = 
            str_extract(attributes, 'G\\d+\\.\\d+'),  
        feelnc_biotype = 
            str_extract(attributes, 'feelnc_biotype \".+\"') %>% 
            str_replace_all('feelnc_biotype ', '') %>% 
            str_replace_all('\"', '')) %>%  
    distinct(old_id, feelnc_biotype) %>%  
    right_join(new2old)

ens <- 
    ens %>% 
    left_join(new2old) %>% 
    write_tsv('chicken_ensembl_all_info.bed')
