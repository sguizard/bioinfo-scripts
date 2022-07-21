library(getopt)
library(tidyverse)

`%nin%` <- Negate(`%in%`)

spec <- matrix(c(
    "help"  , "h", 0,  "logical",
    "data"  , "d", "", "character",
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
if (is.null(opt$genome)) {
    cat("MISSING OPTION: --genome \n")
    q(status = 1)
}
if (opt$genome %nin% c("pig", "gg6", "gg7b", "gg7w")) {
    cat("OPTION ERROR: --genome accepted values pig, gg6, gg7b, gg7w")
    q(status = 1)
}

specie        <- ""
ensembl       <- 0
kuo           <- 0
stages        <- c("")
timepoints    <- c("")
tissues       <- rev(c("Brain", "Ileum", "Kidney", "Liver", "Lung", "Muscle", "Skin"))
match_tissues <- "_(Brain|Ileum|Kidney|Liver|Lung|Muscle|Skin)"
feelnc_cat    <- rev(c("mRNA", "lncRNA", "TUCp", "noORF"))

if (opt$genome == "pig") {
    specie  <- "Sscrofa11.1"
    ensembl <- 31908
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
    xlim_02 <- c(0, 540000)
    xlim_03 <- c(0, 225000)
    xlim_07 <- c(0, 190000)
    xlim_08 <- c(0, 190000)
    xlim_09 <- c(0, 30000)
    xlim_10 <- c(0, 31000)
    xlim_11 <- c(0, 50000)
    xlim_12 <- c(0, 125000)
    xlim_13 <- c(0, 21000)

} else if (startsWith(opt$genome, "gg")) {
    ensembl <- 24356
    kuo     <- 29013
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
    xlim_02 <- c(0, 400000)
    xlim_03 <- c(0, 160000)
    xlim_07 <- c(0, 100000)
    xlim_08 <- c(0, 100000)
    xlim_09 <- c(0, 33000)
    xlim_10 <- c(0, 33000)
    xlim_11 <- c(0, 63000)
    xlim_12 <- c(0, 40000)
    xlim_13 <- c(0, 18000)

    if (opt$genome == "gg6") {
        specie  <- "GRCg6a"
    } else if (opt$genome == "gg7b") {
        specie  <- "GRCg7b"
    } else if (opt$genome == "gg7w") {
        specie  <- "GRCg7w"
    }
}

d <- read_tsv(file = opt$data)

### Total Number of transcripts
cat("01 - Ploting Total Number of Transcripts...\n")
if (specie == "Sscrofa11.1") {
d %>%
    filter(trans_read_count >= 2) %>%
    select(
        tama_gene_id,   tama_transcript_id, ens_gene,
        ens_transcript, sources,    old_id, trans_read_count) %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    distinct(tama_transcript_id) %>%
    count() %>%
    mutate(src = "GENE-SWitCH") %>%
    add_row(src = "Ensembl 105", n = ensembl) %>%
    mutate(src = factor(src, levels = c("Ensembl 105", "GENE-SWitCH"))) %>%
    ggplot(aes(n, src)) +
    geom_col() +
    geom_text(aes(label = n), hjust = -0.1) +
    labs(
        title = paste0(specie, " - Total number of transcripts"),
        x = "Number of transcripts",
        y = "source") +
    xlab(label = "Number of transcripts") +
    xlim(xlim_02)
} else {
d %>%
    filter(trans_read_count >= 2) %>%
    select(
        tama_gene_id,   tama_transcript_id, ens_gene,
        ens_transcript, sources,    old_id, trans_read_count) %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    distinct(tama_transcript_id) %>%
    count() %>%
    mutate(src = "GENE-SWitCH") %>%
    add_row(src = "Ensembl 105", n = ensembl) %>%
    add_row(src = "Kuo et al 2017", n = kuo) %>%
    mutate(src = factor(src, levels = c("Kuo et al 2017", "Ensembl 105", "GENE-SWitCH"))) %>%
    ggplot(aes(n, src)) +
    geom_col() +
    geom_text(aes(label = n), hjust = -0.1) +
    labs(
        title = paste0(specie, " - Total number of transcripts"),
        x = "Number of transcripts",
        y = "source") +
    xlab(label = "Number of transcripts") +
    xlim(xlim_02)
}

ggsave(
    paste0("TRANSCRIPTS_01-geom_col_-_", specie, "_-_FILTERED_-_Total_number_of_transcripts.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


### Total Number of transcripts Ensembl
cat("02 - Ploting Total Number of Transcripts Ensembl...\n")
if (specie == "Sscrofa11.1") {
    d %>%
        filter(trans_read_count >= 2|is.na(trans_read_count)) %>%
        select(
            tama_gene_id,   tama_transcript_id, ens_gene,
            ens_transcript, sources,    old_id, trans_read_count) %>%
        filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
        distinct(tama_transcript_id, sources) %>%
        group_by(tama_transcript_id) %>%
        mutate(count = str_count(sources, ",")) %>%
        top_n(n = 1, wt = count) %>%
        ungroup() %>%
        count(sources) %>%
        mutate(sources = str_replace(sources, "ensembl,isoseq", "ensembl")) %>%
        rename(known = sources) %>%
        mutate(sources = "GENE-SWitCH") %>%
        add_row(sources = "Ensembl 105", known = "ensembl", n = ensembl) %>%
        mutate(sources = factor(
            sources,
            levels = c("Ensembl 105", "GENE-SWitCH"))) %>%
        ggplot(aes(n, sources, fill = known)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = n), hjust = -0.1, position = position_dodge(0.9)) +
        xlim(xlim_02) +
        theme(legend.title = element_blank()) +
        labs(
            title = paste0(specie, " - Total number of transcripts"),
            x = "Number of transcripts",
            y = "source")
} else {
    d %>%
        filter(trans_read_count >= 2|is.na(trans_read_count)) %>%
        select(
            tama_gene_id,   tama_transcript_id, ens_gene,
            ens_transcript, sources,    old_id, trans_read_count) %>%
        filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
        distinct(tama_transcript_id, sources) %>%
        group_by(tama_transcript_id) %>%
        mutate(count = str_count(sources, ",")) %>%
        top_n(n = 1, wt = count) %>%
        ungroup() %>%
        count(sources) %>%
        mutate(sources = str_replace(sources, "ensembl,isoseq", "ensembl")) %>%
        rename(known = sources) %>%
        mutate(sources = "GENE-SWitCH") %>%
        add_row(sources = "Ensembl 105",    known = "ensembl",   n = ensembl) %>%
        add_row(sources = "Kuo et al 2017", known = "Kuo et al", n = kuo) %>%
        mutate(sources = factor(
            sources,
            levels = c("Kuo et al 2017", "Ensembl 105", "GENE-SWitCH"))) %>%
        ggplot(aes(n, sources, fill = known)) +
        geom_col(position = "dodge") +
        geom_text(aes(label = n), hjust = -0.1, position = position_dodge(0.9)) +
        xlim(xlim_02) +
        theme(legend.title = element_blank()) +
        labs(
            title = paste0(specie, " - Total number of transcripts"),
            x = "Number of transcripts",
            y = "source")
}

ggsave(
    paste0("TRANSCRIPTS_02-geom_col_-_", specie, "_-_FILTERED_-_Total_number_of_transcript_ensembl.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


### Read support
cat("03 - Ploting Read Support...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_transcript_id, sources, trans_read_count) %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    na.omit() %>%
    mutate(cat =
        ifelse(trans_read_count == 1, "=1",
        ifelse(trans_read_count == 2, "=2",
        ifelse(between(trans_read_count, 3, 10), "[3,10]",
        ifelse(between(trans_read_count, 11, 100), "[11,100]", ">100"))))) %>%
    count(cat) %>%
    group_by(cat) %>%
    summarise(n = sum(n)) %>%
    mutate(cat = factor(cat, levels = rev(c("=1", "=2", "[3,10]", "[11,100]", ">100")))) %>%
    ggplot(aes(n, cat)) +
    geom_col() +
    scale_fill_viridis_d() +
    labs(
        title = paste0(specie, " - Number of reads supporting transcripts"),
            x = "Number of genes",
            y = "Number of reads") +
    geom_text(aes(label = n), hjust = -0.1) +
    xlim(xlim_03)

ggsave(
    paste0("TRANSCRIPTS_03-geom_col_-_", specie, "_-_FILTERED_-_Number_of_reads_supporting_transcripts.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


### Read support Exons
cat("04 - Ploting Read Support Exons...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_transcript_id, sources, trans_read_count, blockCount) %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    na.omit() %>%
    mutate(cat =
        ifelse(trans_read_count == 1, "=1",
        ifelse(trans_read_count == 2, "=2",
        ifelse(between(trans_read_count, 3, 10), "[3,10]",
        ifelse(between(trans_read_count, 11, 100), "[11,100]", ">100"))))) %>%
    count(cat, blockCount) %>%
    mutate(blockCount = ifelse(blockCount > 15, ">15", blockCount)) %>%
    group_by(cat, blockCount) %>%
    summarise(n = sum(n)) %>%
    mutate(
        cat = factor(cat, levels = rev(c("=1", "=2", "[3,10]", "[11,100]", ">100"))),
        Exons = factor(
            blockCount,
            levels = rev(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", ">15")))) %>%
    ggplot(aes(n, cat, fill = Exons)) +
    geom_col() +
    scale_fill_viridis_d() +
    labs(
        title = paste0(specie, " - Number of reads supporting transcripts"),
            x = "Number of genes",
            y = "Number of reads") +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 10))

ggsave(
    paste0("TRANSCRIPTS_04-geom_col_-_", specie, "_-_FILTERED_-_Number_of_reads_supporting_transcripts_exons.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


### Transcripts per Dev. stage
cat("05 - Ploting Transcripts per Dev. stage...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_transcript_id, source_line) %>%
    filter(!is.na(source_line)) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    mutate(Stage = factor(str_replace(source_line, match_tissues, ""), levels = stages)) %>%
    distinct(tama_transcript_id, Stage) %>%
    count(Stage) %>%
    ggplot(aes(n, Stage)) +
    geom_col() +
    geom_text(aes(label = n), hjust = 1.1) +
    labs(
        title = paste0(specie, " - Number of transcripts per Dev. stage"),
        x = "Number of transcripts",
        y = "Development stage")

ggsave(
    paste0("TRANSCRIPTS_05-geom_col_-_", specie, "_-_FILTERED_-_Number_of_transcripts_per_dev_stage.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


### Transcripts per Tissue
cat("06 - Ploting Transcripts per Tissue...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_transcript_id, source_line) %>%
    filter(!is.na(source_line)) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    mutate(Tissue = factor(str_replace(source_line, match_stages, ""), levels = tissues)) %>%
    distinct(tama_transcript_id, Tissue) %>%
    count(Tissue) %>%
    ggplot(aes(n, Tissue)) +
    geom_col() +
    geom_text(aes(label = n), hjust = 1.1) +
    labs(
        title = paste0(specie, " - Number of transcripts per Tissue"),
        x = "Number of transcripts",
        y = "Tissue")

ggsave(
    paste0("TRANSCRIPTS_06-geom_col_-_", specie, "_-_FILTERED_-_Number_of_transcripts_per_tissue.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


### Transcripts per Timepoint
cat("07 - Ploting Transcripts per Timepoint...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_transcript_id, source_line) %>%
    filter(!is.na(source_line)) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    count(source_line) %>%
    mutate(
        source_line = factor(source_line,                                levels = timepoints),
        Tissue      = factor(str_replace(source_line, match_stages, ""), levels = tissues)) %>%
    ggplot(aes(n, source_line, fill = Tissue)) +
    geom_col() +
    geom_text(aes(label = n), hjust = 1.1) +
    labs(
        title = paste0(specie, " - Number of transcripts per Timepoint"),
        x = "Number of transcripts",
        y = "Tissue") +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 10)) +
    xlim(xlim_07)

    ggsave(
        paste0("TRANSCRIPTS_07-geom_col_-_", specie, "_-_FILTERED_-_Number_of_transcripts_per_timepoint.png"),
        units = "px",
        width = 755,
        height = 610,
        dpi = 125)


# Gene per Timepoint (Ensembl)
cat("08 - Ploting Gene per Timepoint (Ensembl)...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_transcript_id, source_line, sources) %>%
    filter(!is.na(source_line)) %>%
    separate_rows(source_line, sep = ',') %>%
    mutate(source_line = str_replace(source_line, '_P[12]', '')) %>%
    distinct() %>%
    count(source_line, sources) %>%
    mutate(source_line = factor(source_line, levels = timepoints)) %>%
    mutate(Tissue      = factor(str_replace(source_line, match_stages, ''), levels = tissues)) %>%
    mutate(sources = ifelse(sources == 'isoseq', sources, 'ensembl')) %>%
    ggplot(aes(n, source_line, fill = sources)) +
        geom_col(position = "dodge") +
        geom_text(
            aes(label = n),
            hjust = -0.1,
            position = position_dodge(1),
            size = 3) +
        labs(
            title = paste0(specie, " - Number of transcripts per Timepoint (Ensembl)"),
            x = "Number of transcripts",
            y = "Timepoints") +
        xlim(xlim_08)

ggsave(
    paste0("TRANSCRIPTS_08-geom_col_-_", specie, "_-_FILTERED_-_Number_of_transcripts_per_timepoint_ensembl.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


### Feelnc annotation
cat("09 - Ploting Feelnc annotation...\n")
d %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    filter(trans_read_count >= 2) %>%
    count(feelnc_biotype) %>%
    filter(!is.na(feelnc_biotype)) %>%
    mutate(feelnc_biotype = factor(feelnc_biotype, levels = feelnc_cat)) %>%
    ggplot(aes(n, feelnc_biotype)) +
    geom_col() +
    geom_text(aes(label = n), hjust = -0.1) +
    xlim(xlim_09) +
    labs(
        title = paste0(specie, " - Feelnc lncRNA annotation"),
        x = "Number of transcripts",
        y = element_blank())

ggsave(
    paste0("TRANSCRIPTS_09-geom_col_-_", specie, "_-_FILTERED_-_Feelnc_annotation.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


### Feelnc annotation - Ensembl
cat("10 - Ploting Feelnc annotation - Ensembl...\n")
d %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    filter(trans_read_count >= 2) %>%
    count(sources, feelnc_biotype) %>%
    filter(!is.na(feelnc_biotype)) %>%
    add_row(sources = "ensembl,isoseq", feelnc_biotype = "noORF", n = 0) %>%
    mutate(feelnc_biotype = factor(feelnc_biotype, levels = feelnc_cat)) %>%
    mutate(sources = ifelse(sources == "isoseq", sources, "ensembl")) %>%
    ggplot(aes(n, feelnc_biotype, fill = sources)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = n), position = position_dodge(0.9), hjust = -0.1) +
    xlim(xlim_10) +
    labs(
        title = paste0(specie, " - Feelnc lncRNA annotation (Known genes)"),
        x = "Number of transcripts",
        y = element_blank()) +
    theme(legend.position = "bottom", legend.title = element_blank())

ggsave(
    paste0("TRANSCRIPTS_10-geom_col_-_", specie, "_-_FILTERED_-_Feelnc_annotation_ensembl.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


## Feelnc annotation - read_support
cat("11 - Ploting Feelnc annotation - read_support...\n")
d %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    mutate(read_count = ifelse(trans_read_count == 1, "=1", ">=2")) %>%
    count(read_count, feelnc_biotype) %>%
    filter(!is.na(feelnc_biotype)) %>%
    add_row(read_count = ">=2", feelnc_biotype = "noORF", n = 0) %>%
    mutate(feelnc_biotype = factor(feelnc_biotype, levels = feelnc_cat)) %>%
    ggplot(aes(n, feelnc_biotype, fill = read_count)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = n), position = position_dodge(0.9), hjust = -0.1) +
    labs(
        title = paste0(specie, " - Feelnc lncRNA annotation (Known genes)"),
        x = "Number of transcripts",
        y = element_blank()) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlim(xlim_11)

ggsave(
    paste0("TRANSCRIPTS_11-geom_col_-_", specie, "_-_FILTERED_-_Feelnc_annotation_read_support.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


### Number of exons per transcript
cat("12 - Ploting Number of exons per transcript...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    mutate(blockCount = ifelse(blockCount >= 30, ">=30", blockCount)) %>%
    count(blockCount) %>%
    mutate(blockCount = factor(
		blockCount,
		levels = rev(c(seq(from = 1, to = 29, by = 1), ">=30")))) %>%
    ggplot(aes(n, blockCount)) +
    geom_col() +
    labs(
        title = paste0(specie, " - Number of exons per transcript"),
        x     = "Number of transcripts",
        y     = "Number of exons") +
    geom_text(aes(label = n), hjust = -0.1) +
    xlim(xlim_12)

ggsave(
    paste0("TRANSCRIPTS_12-geom_col_-_", specie, "_-_FILTERED_-_Exons_per_transcripts.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


### Number of specific transcripts
cat("13 - Ploting Number of specific transcripts...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_transcript_id, source_line) %>%
    filter(!is.na(source_line)) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    group_by(tama_transcript_id) %>%
    summarise(source_line = paste0(source_line, collapse = ",")) %>%
    mutate(nSamples = str_count(source_line, ",")) %>%
    filter(nSamples == 0) %>%
    count(source_line) %>%
    mutate(
        source_line = factor(source_line,                                levels = timepoints),
        Tissue      = factor(str_replace(source_line, match_stages, ""), levels = tissues)) %>%
    ggplot(aes(n, source_line, fill = Tissue)) +
    geom_col() +
    geom_text(aes(label = n), hjust = -0.1) +
    labs(
        title = paste0(specie, " - Number of specific transcripts"),
        x = "Number of specific transcripts",
        y = "Timepoint") +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 10)) +
    xlim(xlim_13)

ggsave(
    paste0("TRANSCRIPTS_13-geom_col_-_", specie, "_-_FILTERED_-_Number_of_specific_transcripts.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)