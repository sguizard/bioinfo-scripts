library(tidyverse)

h <- 
    read_tsv("high_conf.bed") %>% 
    mutate(conf = "high")

l <- 
    read_tsv("low_conf.bed") %>% 
    mutate(conf = "low") 

hl <- 
    bind_rows(h,l) %>% 
    separate_rows(support_line, sep = ";") %>%  
    separate(support_line, into = c("sample", "reads"), sep = ":") %>%  
    separate_rows(reads, sep = ",")

count <-  
    read_tsv("31SamplesFiles_cluster_read_support.txt") %>%  
    select("sample", "reads" = "merge_gene_id", "trans_read_count")

trc <-  
    hl %>%  
    left_join(count) %>% 
    mutate(trans_read_count = ifelse(str_detect(reads, "ccs"), 1, trans_read_count)) %>%  
    group_by(merge_gene_id, merge_trans_id) %>%  
    summarise(trans_read_count = sum(trans_read_count))

grc <-  
    trc %>%  
    group_by(merge_gene_id) %>% 
    summarise(gene_read_count = sum(trans_read_count))

corrected_count <- trc %>% left_join(grc)

cnc <-  
    bind_rows(h,l) %>%  
    select(-support_line) %>%  
    left_join(corrected_count) %>%  
    arrange("merge_gene_id", "merge_trans_id") %>%  
    write_tsv("annotation_conf_and_count.tsv", col_names = TRUE)

cncg <-  
    cnc %>%  
    select(merge_gene_id, gene_read_count) %>%  
    distinct()

cncg <-  
    cnc %>%  
    group_by(merge_gene_id) %>%  
    count(conf) %>%  
    pivot_wider(names_from = conf, values_from = n) %>%  
    replace_na(replace = list(high = 0, low = 0)) %>% 
    right_join(cncg)

gtf <- read_tsv("31SamplesFilesConf.gtf")

gtf_corrected <- 
    gtf %>%     
    left_join(cncg) %>%     
    left_join(cnc) %>%     
    select(-conf, -merge_gene_id, -merge_trans_id) %>%     
    mutate(    
        attributes       = str_replace_all(attributes,       '""', '"'                     ) %>% str_replace_all(';$', '' ),   
        support_line     = str_replace_all(support_line,     '^',  'support_line "'        ) %>% str_replace_all('$',  '"'),   
        confidence       = str_replace_all(confidence,       '^',  'confidence "'          ) %>% str_replace_all('$',  '"'),   
        gene_read_count  = str_replace_all(gene_read_count,  '^',  'gene_read_count "'     ) %>% str_replace_all('$',  '"'),   
        trans_read_count = str_replace_all(trans_read_count, '^',  'trans_read_count "'    ) %>% str_replace_all('$',  '"'), 
        high             = str_replace_all(high,             '^',  'transcript_low_count "' ) %>% str_replace_all('$',  '"'), 
        low              = str_replace_all(low,              '^',  'transcript_high_count "') %>% str_replace_all('$',  '"')) %>%   
    unite(attributes, high, low, support_line, confidence, gene_read_count, trans_read_count, col = attributes, sep = "; ") %>%   
    mutate(attributes = str_replace_all(attributes, '; NA', '')) %>%  
    write_tsv("31SamplesFilesConfCountCorrected.gtf", escape = "none")
