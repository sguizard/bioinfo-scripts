library(tidyverse)

read_gff <- function(file) { 
    require(readr) 
    read_tsv( 
        file,  
        col_names = c( 
            "sequence",  
            "source",  
            "feature",  
            "start",  
            "end",  
            "score",  
            "strand",  
            "phase",  
            "attributes"),  
        col_types = cols( 
            sequence   = col_character(),  
            source     = col_character(),  
            feature    = col_character(),  
            start      = col_integer(),  
            end        = col_integer(),  
            score      = col_double(),  
            strand     = col_character(),  
            phase      = col_character(),  
            attributes = col_character()),  
    comment = "#", skip = 1) 
}

gff     <- read_gff('31SamplesFilesConfCountCorrected.gtf')
tgi_egi <- read_tsv('tama_gene_id_to_ens_gene_id.tsv')
tti_eti <- read_tsv("tama_transcript_id_to_ens_transcript_id.tsv")

gff %>%  
    mutate( 
        tama_gene_id       = str_extract(attributes, 'gene_id "[^"]+";')       %>% str_extract('G[^"]+'),  
        tama_transcript_id = str_extract(attributes, 'transcript_id "[^"]+";') %>% str_extract('G[^"]+')) %>%  
    left_join(tgi_egi) %>%  
    left_join(tti_eti) %>%  
    select(-tama_gene_id, -tama_transcript_id) %>%  
    mutate( 
        ens_gene_name       = paste0('ens_gene_name "', ens_gene_name, '"'),  
        ens_transcript_name = paste0('ens_transcript_name "', ens_transcript_name, "'"),  
        attributes          = paste0(attributes, "; ", ens_gene_name, "; ", ens_transcript_name)) %>%  
    select(-ens_gene_name, -ens_transcript_name) %>%  
    write_tsv("31SamplesFilesConfCountCorrectedEnsNames.gtf")
