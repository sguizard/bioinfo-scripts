#!/usr/bin/env Rscript

library(getopt)
library(tidyverse)
library(viridis)
library(this.path)
source(paste0(this.dir(), "/../utils/readGxf.R"))

`%nin%` <- Negate(`%in%`)

spec <- matrix(c(
    "help", "h", 0,  "logical",
    "data", "d", "", "character",
    "ensembl", "e", "", "character",
    "genome", "g", "", "character"
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}
if (is.null(opt$data)) {
    cat("MISSING OPTION: --data \n")
    q(status = 1)
}
if (is.null(opt$ensembl)) {
    cat("MISSING OPTION: --ensembl \n")
    q(status = 1)
}
if (is.null(opt$genome)) {
    cat("MISSING OPTION: --genome \n")
    q(status = 1)
}
if (opt$genome %nin% c("pig", "gg6", "gg7b", "gg7w")) {
    cat("OPTION ERROR: --genome accepted values pig, gg6, gg7b, gg7w")
    q(status = 1)
}

specie        <- ""
stages        <- c("")
timepoints    <- c("")
tissues       <- rev(c(
    "Brain", "Ileum", "Kidney", "Liver",
    "Lung", "Muscle", "Skin"))
match_tissues <- "_(Brain|Ileum|Kidney|Liver|Lung|Muscle|Skin)"
feelnc_cat    <- rev(c("mRNA", "TUCp", "lncRNA", "noORF"))

if (opt$genome == "pig") {
    specie  <- "Sscrofa11.1"
    stages  <- rev(c("30dpf", "70dpf", "NB"))
    timepoints <- rev(c(
        "30dpf_Brain",  "70dpf_Brain",  "NB_Brain",
        "30dpf_Ileum",  "70dpf_Ileum",  "NB_Ileum",
        "30dpf_Kidney", "70dpf_Kidney", "NB_Kidney",
        "30dpf_Liver",  "70dpf_Liver",  "NB_Liver",
        "30dpf_Lung",   "70dpf_Lung",   "NB_Lung",
        "30dpf_Muscle", "70dpf_Muscle", "NB_Muscle",
        "30dpf_Skin",   "70dpf_Skin",   "NB_Skin"))
    match_stages <- "(30dpf|70dpf|NB)_"

} else if (startsWith(opt$genome, "gg")) {
    stages <- rev(c("E8", "E15", "HC"))
    timepoints <- rev(c(
        "E8_Brain",  "E15_Brain",  "HC_Brain",
        "E8_Ileum",  "E15_Ileum",  "HC_Ileum",
        "E8_Kidney", "E15_Kidney", "HC_Kidney",
        "E8_Liver",  "E15_Liver",  "HC_Liver",
        "E8_Lung",   "E15_Lung",   "HC_Lung",
        "E8_Muscle", "E15_Muscle", "HC_Muscle",
        "E8_Skin",   "E15_Skin",   "HC_Skin"))
    match_stages <- "(E8|E15|HC)_"

    if (opt$genome == "gg6") {
        specie  <- "GRCg6a"
    } else if (opt$genome == "gg7b") {
        specie  <- "GRCg7b"
    } else if (opt$genome == "gg7w") {
        specie  <- "GRCg7w"
    }
}



# Load data
cat("00 - Loading data...\n")
d <- read_tsv(file = opt$data)

cat("00.1 - Counting transcripts...\n")
ensembl <-
    read_gtf(opt$ensembl) %>%
    filter(feature == "transcript") %>%
    nrow()

cat("00.2 - Extract exons...\n")
ensembl_exon <-
    read_gtf(file = opt$ensembl, separate_attributes = TRUE) %>%
    group_by(transcript_id) %>% 
    count(name = "Exon") %>% 
    ungroup() %>% 
    count(Exon) %>% 
    mutate(Project = "Ensembl", Source = "Ensembl 108") %>% 
    select(Project, Source, Exon, n)

ensembl_exon2 <-
    read_gtf(opt$ensembl, separate_attributes = TRUE) %>%
    filter(feature == "exon") %>%
    count(transcript_id, name = "Exon") %>% 
    rename(merge_trans_id = transcript_id) %>% 
    mutate(Project = "Ensembl 108")

cat("00.3 - Counting transcripts per gene...\n")
ensembl_transcripts_by_gene <-
    read_gtf(file = opt$ensembl, separate_attributes = TRUE) %>%
    filter(feature == "transcript") %>%
    count(gene_id) %>%
    mutate(Transcript = case_when(
        n == 1 ~ "=1",
        n == 2 ~ "=2",
        between(n, 3, 10) ~ "[3,10]",
        between(n, 11, 100) ~ "[11,100]",
        .default = ">100")) %>%
    mutate(p = n / sum(n) * 100) %>%
    group_by(Transcript) %>%
    summarise(
        n = sum(n),
        p = sum(p)) %>%
    mutate(
        Source = "Ensembl 108",
        Project = "Ensembl")




### Total Number of transcripts
cat("01 - Ploting Total Number of Transcripts...\n")
d %>%
    select(merge_trans_id, sources_id_gene) %>%
    filter(sources_id_gene %in% c("ISOseq", "Ensembl,ISOseq")) %>%
    distinct(merge_trans_id) %>%
    count() %>%
    mutate(src = "GENE-SWitCH") %>%
    add_row(src = "Ensembl 108", n = ensembl) %>%
    mutate(src = factor(src, levels = c("Ensembl 108", "GENE-SWitCH"))) %>%
    ggplot(aes(n, src)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -.1) +
        scale_x_continuous(expand = expansion(mult = .20)) +
        labs(
            title = paste0(specie, " - Total number of transcripts"),
            x = "Number of transcripts",
            y = "")

ggsave(
    paste0(
        "TRANSCRIPTS_01-geom_col_-_", 
        specie, 
        "_-_Total_number_of_transcripts.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


### Total Number of transcripts Ensembl
cat("02 - Ploting Total Number of Transcripts Ensembl...\n")
d %>% 
    select(Source = source_summary) %>%
    mutate(Source = str_replace_all(Source, ":\\d+", "")) %>%
    count(Source, name = "Count") %>%
    add_row(Source = "Ensembl 108", Count = ensembl) %>% 
    mutate(Project =
        ifelse(Source == "Ensembl 108", "Ensembl", "GENE-SWitCH") %>% 
        factor(levels = c("GENE-SWitCH", "Ensembl"))) %>%
    ggplot(aes(Count, Source)) +
        geom_col() +
        facet_grid(Project ~ ., scales = "free") +
        geom_text(aes(label = Count), hjust = -.1) +
        scale_x_continuous(expand = expansion(mult = .20)) +
        labs(
            title = paste0(specie, " - Total number of transcripts"),
            x = "Number of transcripts",
            y = "")

ggsave(
    paste0(
        "TRANSCRIPTS_02-geom_col_-_", 
        specie, 
        "_-_Total_number_of_transcript_ensembl.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)



# Transcript Number Ensembl + Exon number
cat("03 - Ploting Transcript Number Ensembl + Exon number...\n")
cat("03.1 - Ploting Transcript Number Ensembl + Exon number (Count)...\n")
d %>%
    select(Source = source_summary, Exon = blockCount) %>% 
    mutate(Source = str_replace_all(Source, ":\\d+", "")) %>%  
    count(Source, Exon) %>%
    mutate(Project = "GENE-SWitCH") %>% 
    rbind(ensembl_exon) %>% 
    mutate(Project = factor(Project, levels = c("GENE-SWitCH", "Ensembl"))) %>%
    mutate(Exon = case_when(
        Exon == 1 ~ "1",
        Exon == 2 ~ "2",
        Exon == 3 ~ "3",
        Exon == 4 ~ "4",
        Exon == 5 ~ "5",
        Exon == 6 ~ "6",
        Exon == 7 ~ "7",
        Exon == 8 ~ "8",
        Exon == 9 ~ "9",
        Exon >= 10 ~ ">= 10")) %>% 
    mutate(Exon = factor(
        Exon,
        levels = rev(c(
            "1", "2", "3", "4", "5",
            "6", "7", "8", "9", ">= 10")))) %>%
    ggplot(aes(n, Source, fill = Exon)) +
        geom_col() +
        facet_grid(Project ~ ., scales = "free") +
        scale_x_continuous(expand = expansion(mult = .15)) +
        labs(
            title = paste0(
                specie,
                " - Total number of transcript + Number of Exons"),
            x = "Number of transcripts",
            y = "") +
        scale_fill_viridis(discrete = TRUE) +
        theme(legend.position = "bottom") +
        guides(fill = guide_legend(
            nrow = 1,
            reverse = TRUE,
            title.position = "top"))

ggsave(
    paste0(
        "TRANSCRIPTS_03-geom_col_-_",
        specie,
        "_-_Total_number_of_transcript_ensembl_plus_exons.png"),
    units = "px",
    width = 755,
    height = 450,
    dpi = 125)



# Transcript Number Ensembl + Exon number (Percentages)
cat("03.2 - Ploting Transcript Number Ensembl + Exon number (Percentage by Project)...\n")
d %>%
    select(Source = source_summary, Exon = blockCount) %>% 
    mutate(Source = str_replace_all(Source, ":\\d+", "")) %>%  
    count(Source, Exon) %>%
    mutate(Project = "GENE-SWitCH") %>% 
    rbind(ensembl_exon) %>% 
    mutate(Project = factor(Project, levels = c("GENE-SWitCH", "Ensembl"))) %>%
    mutate(Exon = case_when(
        Exon == 1 ~ "1",
        Exon == 2 ~ "2",
        Exon == 3 ~ "3",
        Exon == 4 ~ "4",
        Exon == 5 ~ "5",
        Exon == 6 ~ "6",
        Exon == 7 ~ "7",
        Exon == 8 ~ "8",
        Exon == 9 ~ "9",
        Exon >= 10 ~ ">= 10")) %>% 
    mutate(Exon = factor(
        Exon,
        levels = rev(c(
            "1", "2", "3", "4", "5",
            "6", "7", "8", "9", ">= 10")))) %>%
    group_by(Project) %>%
        mutate(p = n / sum(n) * 100) %>%
    ggplot(aes(p, Source, fill = Exon)) +
        geom_col() +
        facet_grid(Project ~ ., scales = "free") +
        scale_x_continuous(expand = expansion(mult = .15)) +
        labs(
            title = paste0(specie, " - Total number of transcripts"),
            x = "Number of transcripts + Exon Number",
            y = "") +
        scale_fill_viridis(discrete = TRUE) +
    theme(legend.position = "bottom") +
        guides(fill = guide_legend(
            nrow = 1,
            reverse = TRUE,
            title.position = "top"))

ggsave(
    paste0(
        "TRANSCRIPTS_03-geom_col_-_",
        specie,
        "_-_Total_number_of_transcripts_ensembl_plus_exons_percentages_by_project.png"),
    units = "px",
    width = 755,
    height = 450,
    dpi = 125)



cat("03.3 - Ploting Transcript Number Ensembl + Exon number (Percentage by source)...\n")
d %>%
    select(Source = source_summary, Exon = blockCount) %>% 
    mutate(Source = str_replace_all(Source, ":\\d+", "")) %>%  
    count(Source, Exon) %>%
    mutate(Project = "GENE-SWitCH") %>% 
    rbind(ensembl_exon) %>% 
    mutate(Project = factor(Project, levels = c("GENE-SWitCH", "Ensembl"))) %>%
    mutate(Exon = case_when(
        Exon == 1 ~ "1",
        Exon == 2 ~ "2",
        Exon == 3 ~ "3",
        Exon == 4 ~ "4",
        Exon == 5 ~ "5",
        Exon == 6 ~ "6",
        Exon == 7 ~ "7",
        Exon == 8 ~ "8",
        Exon == 9 ~ "9",
        Exon >= 10 ~ ">= 10")) %>% 
    mutate(Exon = factor(
        Exon,
        levels = rev(c(
            "1", "2", "3", "4", "5",
            "6", "7", "8", "9", ">= 10")))) %>%
    group_by(Source) %>%
        mutate(p = n / sum(n) * 100) %>%
    ggplot(aes(p, Source, fill = Exon)) +
        geom_col() +
        facet_grid(Project ~ ., scales = "free") +
        scale_x_continuous(expand = expansion(mult = .15)) +
        labs(
            title = paste0(specie, " - Total number of transcripts"),
            x = "Number of transcripts + Exon Number",
            y = "") +
        scale_fill_viridis(discrete = TRUE) +
    theme(legend.position = "bottom") +
        guides(fill = guide_legend(
            nrow = 1,
            reverse = TRUE,
            title.position = "top"))

ggsave(
    paste0(
        "TRANSCRIPTS_03-geom_col_-_",
        specie,
        "_-_Total_number_of_transcripts_ensembl_plus_exons_percentages_by_source.png"),
    units = "px",
    width = 755,
    height = 450,
    dpi = 125)



### Read support
cat("04.1 - Ploting Read Support...\n")
d %>% 
    select(merge_trans_id, trans_read_count) %>% 
    mutate(cat =
        ifelse(trans_read_count == 1, "=1",
        ifelse(trans_read_count == 2, "=2",
        ifelse(between(trans_read_count, 3, 10), "[3,10]",
        ifelse(between(trans_read_count, 11, 100), "[11,100]", ">100"))))) %>% 
    count(cat) %>%
    mutate(cat = factor(
        cat,
        levels = rev(c("=1", "=2", "[3,10]", "[11,100]", ">100")))) %>%
    ggplot(aes(n, cat)) +
        geom_col() +
        labs(
            title = paste0(specie, " - Number of reads supporting transcripts"),
            x = "Number of transcripts",
            y = "Number of reads") +
        geom_text(aes(label = n), hjust = -.1) +
        scale_x_continuous(expand = expansion(mult = .16))

ggsave(
    paste0(
        "TRANSCRIPTS_04-geom_col_-_", 
        specie, 
        "_-_Number_of_reads_supporting_transcripts.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)



### Read support Exons
cat("04.2 - Ploting Read Support Exons...\n")
d %>% 
    select(merge_trans_id, trans_read_count, Exon = blockCount) %>% 
    mutate(cat =
        ifelse(trans_read_count == 1, "=1",
        ifelse(trans_read_count == 2, "=2",
        ifelse(between(trans_read_count, 3, 10), "[3,10]",
        ifelse(between(trans_read_count, 11, 100), "[11,100]", ">100"))))) %>% 
    count(Exon, cat) %>%
    mutate(
        cat = factor(
            cat,
            levels = rev(c("=1", "=2", "[3,10]", "[11,100]", ">100"))),
        Exon =
            ifelse(
                Exon >= 10,
                ">= 10",
                as.character(Exon)) %>%
            factor(levels = rev(c(
                "1", "2", "3", "4", "5",
                "6", "7", "8", "9", ">= 10")))) %>%
    ggplot(aes(n, cat, fill = Exon)) +
        geom_col() +
        labs(
            title = paste0(specie, " - Number of reads supporting transcripts"),
            x = "Number of transcripts",
            y = "Number of reads") +
        scale_x_continuous(expand = expansion(mult = .05)) +
        scale_fill_viridis(discrete = TRUE) +
        theme(legend.position = "bottom") +
        guides(fill = guide_legend(reverse = TRUE, nrow = 1))

ggsave(
    paste0(
        "TRANSCRIPTS_04-geom_col_-_", 
        specie, 
        "_-_Number_of_reads_supporting_transcripts_exons.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)



### Transcripts per Dev. stage
cat("05 - Ploting Transcripts per Dev. stage...\n")
d %>%
    select(merge_trans_id, source_line) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    mutate(
        source_line = factor(source_line, levels = timepoints),
        Tissue = factor(
            str_replace(source_line, match_stages, ""),
            levels = tissues),
        Stage = factor(
            str_replace(source_line, match_tissues, ""),
            levels = stages)) %>%
    distinct(merge_trans_id, Stage) %>%
    count(Stage) %>%
    ggplot(aes(n, Stage)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -.1) +
        labs(
            title = paste0(specie, " - Number of transcripts per Dev. Stage"),
            x = "Number of transcripts",
            y = "Development stage") +
        scale_x_continuous(expand = expansion(mult = .15))

ggsave(
    paste0(
        "TRANSCRIPTS_05-geom_col_-_", 
        specie, 
        "_-_Number_of_transcripts_per_dev_stage.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


### Transcripts per Tissue
cat("06 - Ploting Transcripts per Tissue...\n")
d %>%
    select(merge_trans_id, source_line) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    mutate(
        source_line = factor(source_line, levels = timepoints),
        Tissue = factor(
            str_replace(source_line, match_stages, ""),
            levels = tissues),
        Stage = factor(
            str_replace(source_line, match_tissues, ""),
            levels = stages)) %>%
    distinct(merge_trans_id, Tissue) %>%
    count(Tissue) %>%
    mutate(Tissue = fct_reorder(Tissue, n)) %>%
    ggplot(aes(n, Tissue)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -.1) +
        labs(
            title = paste0(specie, " - Number of transcripts per Dev. Stage"),
            x = "Number of transcripts",
            y = "Development stage") +
        scale_x_continuous(expand = expansion(mult = .15))

ggsave(
    paste0(
        "TRANSCRIPTS_06-geom_col_-_", 
        specie, 
        "_-_Number_of_transcripts_per_tissue.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)



### Transcripts per Timepoint
cat("07 - Ploting Transcripts per Timepoint...\n")
d %>%
    select(merge_trans_id, source_line) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    mutate(
        source_line = factor(source_line, levels = timepoints),
        Tissue = factor(
            str_replace(source_line, match_stages, ""),
            levels = tissues)) %>% 
    count(source_line, Tissue) %>% 
    ggplot(aes(n, source_line, fill = Tissue)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -.1) +
        labs(
            title = paste0(specie, " - Number of transcripts per Timepoint"),
            x = "Number of transcripts",
            y = "Timepoint") +
        scale_x_continuous(expand = expansion(mult = .20)) +
        guides(fill = guide_legend(reverse = TRUE))

    ggsave(
        paste0(
            "TRANSCRIPTS_07-geom_col_-_",
            specie,
            "_-_Number_of_transcripts_per_timepoint.png"),
        units = "px",
        width = 755,
        height = 610,
        dpi = 125)



### Feelnc annotation
cat("08 - Ploting Feelnc annotation...\n")
d %>% 
    select(merge_trans_id, feelnc_biotype) %>% 
    count(feelnc_biotype) %>% 
    filter(!is.na(feelnc_biotype)) %>% 
    mutate(feelnc_biotype = factor(feelnc_biotype, levels = feelnc_cat)) %>% 
    ggplot(aes(n, feelnc_biotype)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -.1) +
    labs(
        title = paste0(specie, " - Feelnc lncRNA annotation"),
        x = "Number of transcripts",
        y = element_blank()) +
    scale_x_continuous(expand = expansion(mult = .20)) +
        guides(fill = guide_legend(reverse = TRUE))

ggsave(
    paste0(
        "TRANSCRIPTS_08-geom_col_-_", 
        specie, 
        "_-_Feelnc_annotation.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)



### Number of exons per transcript
cat("09 - Ploting Number of exons per transcript...\n")
d %>% 
    select(merge_trans_id, Exon = blockCount) %>% 
    mutate(Project = "GENE-SWitCH") %>% 
    rbind(ensembl_exon2) %>% 
    mutate(
        Exon = 
            ifelse(Exon >= 20, ">=20", as.character(Exon)) %>% 
            factor(c(seq(1,19,1), ">=20")), 
        Project = factor(
            Project, 
            levels = c("GENE-SWitCH", "Ensembl 108"))) %>% 
    count(Project, Exon) %>% 
    ggplot(aes(Exon, n)) + 
    geom_col() + 
    facet_wrap(vars(Project), nrow = 2) + 
    # geom_text(aes(label = n), hjust = -.1, angle = 90) +
        labs(
            title = paste0(specie, " - Number of Exons per Transcript"),
            x = "Number of Exons",
            y = "Number of transcripts") # +
    # scale_y_continuous(expand = expansion(mult = .50))

ggsave(
    paste0(
        "TRANSCRIPTS_09-geom_col_-_", 
        specie, 
        "_-_Exons_per_transcripts.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


### Number of specific transcripts
cat("10 - Ploting Number of specific transcripts...\n")
d %>% 
    select(merge_trans_id, source_line) %>% 
    separate_rows(source_line, sep = ",") %>% 
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>% 
    group_by(merge_trans_id) %>% 
    summarise(source_line = paste0(source_line, collapse = ",")) %>% 
    mutate(nSamples = str_count(source_line, ",")) %>%
    filter(nSamples == 0) %>%
    count(source_line) %>% 
    mutate(
        source_line = factor(
            source_line,
            levels = timepoints %>% fct_reorder(n)),
        Tissue = factor(
            str_replace(source_line, match_stages, ""),
            levels = tissues)) %>%
    ggplot(aes(n, source_line, fill = Tissue)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -0.1) +
        labs(
            title = paste0(specie, " - Number of specific transcripts"),
            x = "Number of specific transcripts",
            y = "Timepoint") +
        scale_x_continuous(expand = expansion(mult = .15)) +
        guides(fill = guide_legend(reverse = TRUE))

ggsave(
    paste0(
        "TRANSCRIPTS_10-geom_col_-_", 
        specie, 
        "_-_Number_of_specific_transcripts.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)



cat("11 - Ploting Number of Dev. Stage specific transcripts...\n")
d %>%
    select(merge_trans_id, source_line) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    mutate(
        source_line = factor(source_line, levels = timepoints),
        Tissue = factor(
            str_replace(source_line, match_stages, ""),
            levels = tissues),
        Stage = factor(
            str_replace(source_line, match_tissues, ""),
            levels = stages)) %>%
    distinct(merge_trans_id, Stage) %>%
    group_by(merge_trans_id) %>%
    summarise(Stage = paste0(Stage, collapse = ",")) %>%
    mutate(n = str_count(Stage, ",")) %>%
    filter(n == 0) %>%
    count(Stage) %>%
    mutate(Stage = factor(Stage, levels = stages)) %>%
    ggplot(aes(n, Stage)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -.1) +
        labs(
            title = paste0(
                specie,
                " - Number of specific transcipts per Dev. Stage"),
            x = "Number of transcipts",
            y = "Dev. Stage") +
        scale_x_continuous(expand = expansion(mult = .15))

ggsave(
    paste0(
        "TRANSCRIPTS_11-geom_col_-_",
        specie,
        "_-_Number_of_dev_stage_specific_transcipts.png"),
    units = "px",
    width = 755,
    height = 400,
    dpi = 125)



cat("12 - Ploting Number of Tissue specific transcripts...\n")
d %>%
    select(merge_trans_id, source_line) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    mutate(
        source_line = factor(source_line, levels = timepoints),
        Tissue = factor(
            str_replace(source_line, match_stages, ""),
            levels = tissues),
        Stage = factor(
            str_replace(source_line, match_tissues, ""),
            levels = stages)) %>%
    distinct(merge_trans_id, Tissue) %>%
    group_by(merge_trans_id) %>%
    summarise(Tissue = paste0(Tissue, collapse = ",")) %>%
    mutate(n = str_count(Tissue, ",")) %>%
    filter(n == 0) %>%
    count(Tissue) %>%
    mutate(Tissue =
        factor(Tissue, levels = tissues) %>%
        fct_reorder(n)) %>%
    ggplot(aes(n, Tissue)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -.1) +
        labs(
            title = paste0(
                specie,
                " - Number of specific transcripts per Tissue"),
            x = "Number of transcripts",
            y = "Tissue") +
        scale_x_continuous(expand = expansion(mult = .15))

ggsave(
    paste0(
        "TRANSCRIPTS_12-geom_col_-_",
        specie,
        "_-_Number_of_tissue_specific_transcripts.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)
