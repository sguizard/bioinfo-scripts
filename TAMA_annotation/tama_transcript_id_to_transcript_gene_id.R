library(tidyverse)

read_bed12 <- function(file) { 
        require(readr) 
        read_tsv( 
            file,  
            col_names = c( 
                "chrom",  
                "chromStart",  
                "chromEnd",  
                "name",  
                "score",  
                "strand",  
                "thickStart",  
                "thickEnd",  
                "itemRgb",  
                "blockCount",  
                "blockSizes",  
                "blockStarts"),  
            col_types = cols( 
                chrom = col_character(), 
                chromStart = col_integer(), 
                chromEnd = col_integer(), 
                name = col_character(), 
                score = col_double(), 
                strand = col_character(), 
                thickStart = col_integer(), 
                thickEnd = col_integer(), 
                itemRgb = col_character(), 
                blockCount = col_integer(), 
                blockSizes = col_character(), 
                blockStarts = col_character())) 
}

read_bed12('isoseq_ensembl.bed') %>%  
    select(name) %>%  
    mutate( 
        tama_name = str_replace(name, ';ENS.+', ''),  
        ens_transcript_name  = str_replace(name, '^G\\d+;G\\d+\\.\\d+;?', '')) %>% 
    separate(col = tama_name, into = c("tama_gene_id", "tama_transcript_id"), sep = ";") %>%  
    select(-name) %>%  
    separate_rows(ens_transcript_name, sep = ';') %>%  
    select(tama_transcript_id, ens_transcript_name) %>%  
    filter(str_detect(ens_name, 'ENSGALT')) %>%  
    distinct() %>%  
    write_tsv("tama_transcript_id_to_ens_transcript_id.tsv")
