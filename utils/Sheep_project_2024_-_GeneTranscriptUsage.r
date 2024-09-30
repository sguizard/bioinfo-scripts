setwd("/home/sguizard/Work/Projects/Emily_Clark_-_Sheep_data/BlackFace_x_Texel_-_Rambouillet_v3/results/08_flair_collapse")

library(tidyverse)
source("~/Work/Dev/github/sguizard/bioinfo-scripts/utils/readGxf.R")

countGeneAndTranscripts <- function(file) {
    cat("==> Processing ", file, "\n")
    d <-
        read_gtf(file, separate_attributes = TRUE) %>%
        select(gene_id, transcript_id)
    ng <- d %>% distinct(gene_id) %>% nrow()
    nt <- d %>% distinct(transcript_id) %>% nrow()

    tibble(
        file = file,
        run_id = str_replace(file, ".isoforms.gtf", ""),
        nGene = ng,
        nTranscript = nt)
}

plotCounts <- function(stg) {
    cat("==> Processing ", stg, "\n")
    
    p <-
        rs %>%
        filter(stage == stg) %>%
        pivot_longer(nGene:nTranscript, names_to = "feature", values_to = "counts") %>%
        ggplot(aes(time, counts)) +
        geom_col() +
        facet_grid(feature ~ tissue, scales = "free") + 
        ggtitle(stg)

    print(p)
    ggsave(paste0("Gene_Transcript_count_", stg, ".png"), plot = p)
}


s <- 
    read_tsv("sample_ids_hifi.tsv") %>%
    rename(run_id = run_accession) %>%
    mutate(
        time =
            ifelse(time == "1 week", "1w",
            ifelse(time == "8 week", "8w",
            ifelse(time == "100 days", "100d",
            ifelse(time == "23 days", "23d",
            ifelse(time == "35 days", "35d",
            ifelse(time == "0 days", "0d", "sloubi")))))),
        time  = factor(time, levels = c("23d", "35d", "100d", "0d" ,"1w", "8w")),
        stage = 
            ifelse(stage == "lamb", "Lamb", 
            ifelse(stage == "mother", "Mother", 
            ifelse(stage == "embryo", "Embryo", "sloubi"))),
        stage = factor(stage, levels = c("Lamb", "Mother", "Embryo")))

r <- map_df(Sys.glob(paths = "*.gtf"), countGeneAndTranscripts)

rs <- 
    left_join(r, s) %>% 
    write_tsv("Gene_Transcript_counts.tsv")

map(
    rs %>% distinct(stage) %>% mutate(stage = as.character(stage)) %>% pull(stage), 
    plotCounts)
