library(tidyverse)
source("read_gff.R")

d <- 
    read_gff("chicken-21samples.gtf") %>% 
    mutate(
        merge_gene_id = str_extract(attributes, '\"G\\d+\"') %>% str_replace_all(., '"', ''),
        merge_trans_id = str_extract(attributes, '\"G\\d+\\.\\d+\"') %>% str_replace_all(., '"', ''))

lc <- 
    read_tsv("low_conf.bed") %>% 
    mutate(confidence = "low")

hc <- 
    read_tsv("high_conf.bed") %>% 
    mutate(confidence = "high")

confidences <- bind_rows(lc, hc)

dc <-  
    d %>%  
    left_join(confidences) %>% 
    select(-gene_read_count, -trans_read_count) %>% 
    write_tsv("chicken-21samples_conf.gtf")
