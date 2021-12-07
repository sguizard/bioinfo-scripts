library(tidyverse)

d <- read_tsv("31SamplesFiles_read_support.txt")

high_conf <-  
    d %>%  
    filter( 
        trans_read_count >= 2 | 
        ( trans_read_count == 1 & str_detect(support_line, "transcript") ) 
    ) %>% 
    select("merge_gene_id", "merge_trans_id", "support_line") %>% 
    write_tsv("high_conf.bed")

low_conf <-  
    d %>%  
    filter( 
        !( 
            trans_read_count >= 2 | 
            ( trans_read_count == 1 & str_detect(support_line, "transcript") ) 
        ) 
    ) %>% 
    select("merge_gene_id", "merge_trans_id", "support_line") %>% 
    write_tsv("low_conf.bed")
