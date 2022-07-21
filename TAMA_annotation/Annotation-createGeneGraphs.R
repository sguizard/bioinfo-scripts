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
    match_stages <- '(30dpf|70dpf|NB)_'
    xlim_02 <- c(0, 85000)
    xlim_08 <- c(0, 40000)
    xlim_09 <- c(0, 65000)
    xlim_10 <- c(0, 8000)

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
    match_stages <- '(E8|E15|HC)_'
    xlim_08 <- c(0, 17000)
    xlim_09 <- c(0, 12500)
    xlim_10 <- c(0, 1500)

    if (opt$genome == "gg6") {
        specie  <- "GRCg6a"
    } else if (opt$genome == "gg7b") {
        specie  <- "GRCg7b"
    } else if (opt$genome == "gg7w") {
        specie  <- "GRCg7w"
    }
}

# Load data
d <- read_tsv(file = opt$data)

# Gene Number
cat("01 - Ploting Gene Number...\n")
if (specie == "Sscrofa11.1") {
    d %>%
        select(
            tama_gene_id,   tama_transcript_id, ens_gene,
            ens_transcript, sources,    old_id, trans_read_count) %>%
        filter(trans_read_count >= 2) %>%
        filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
        distinct(tama_gene_id) %>%
        count() %>%
        mutate(src = "GENE-SWitCH") %>%
        add_row(src = "Ensembl 105", n = ensembl) %>%
        mutate(src = factor(src, levels = c("Ensembl 105", "GENE-SWitCH"))) %>%
        ggplot(aes(n, src)) +
            geom_col() +
            geom_text(aes(label = n), hjust = 1.1) +
            labs(
                title = paste0(specie, " - Total number of genes"),
                x = "Number of genes",
                y = "source")
} else {
    d %>%
        select(
            tama_gene_id,   tama_transcript_id, ens_gene,
            ens_transcript, sources,    old_id, trans_read_count) %>%
        filter(trans_read_count >= 2) %>%
        filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
        distinct(tama_gene_id) %>%
        count() %>%
        mutate(src = "GENE-SWitCH") %>%
        add_row(src = "Ensembl 105", n = ensembl) %>%
        add_row(src = "Kuo et al 2017", n = kuo) %>%
        mutate(src = factor(src, levels = c("Kuo et al 2017", "Ensembl 105", "GENE-SWitCH"))) %>%
        ggplot(aes(n, src)) +
            geom_col() +
            geom_text(aes(label = n), hjust = 1.1) +
            labs(
                title = paste0(specie, " - Total number of genes"),
                x = "Number of genes",
                y = "source")
}

ggsave(
    paste0("GENE_01-geom_col_-_", specie, "_-_FILTERED_-_Total_number_of_genes.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


# Gene Number Ensembl
cat("02 - Ploting Gene Number Ensembl...\n")
if (specie == "Sscrofa11.1") {
    d %>%
        select(
            tama_gene_id,   tama_transcript_id, ens_gene,
            ens_transcript, sources,    old_id, trans_read_count) %>%
        filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
        filter(trans_read_count >= 2) %>%
        distinct(tama_gene_id, sources) %>%
        group_by(tama_gene_id) %>%
        mutate(count = str_count(sources, ",")) %>%
        top_n(n = 1, wt = count) %>%
        ungroup() %>%
        count(sources) %>%
        mutate(sources = str_replace(sources, "ensembl,isoseq", "ensembl")) %>%
        rename(known = sources) %>%
        mutate(sources = "GENE-SWitCH") %>%
        add_row(sources = "Ensembl 105", known = "ensembl", n = ensembl) %>%
        mutate(sources = factor(sources, levels = c("Ensembl 105", "GENE-SWitCH"))) %>%
        ggplot(aes(n, sources, fill = known)) +
            geom_col(position = "dodge") +
            geom_text(
                aes(label = n),
                hjust = -.1,
                position = position_dodge(0.9)) +
            labs(
                title = paste0(specie, " - Total number of genes"),
                x = "Number of genes",
                y = "source") +
            xlim(xlim_02)
} else {
    d %>%
        select(
            tama_gene_id,   tama_transcript_id, ens_gene,
            ens_transcript, sources,    old_id, trans_read_count) %>%
        filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
        filter(trans_read_count >= 2) %>%
        distinct(tama_gene_id, sources) %>%
        group_by(tama_gene_id) %>%
        mutate(count = str_count(sources, ",")) %>%
        top_n(n = 1, wt = count) %>%
        ungroup() %>%
        count(sources) %>%
        mutate(sources = str_replace(sources, "ensembl,isoseq", "ensembl")) %>%
        rename(known = sources) %>%
        mutate(sources = "GENE-SWitCH") %>%
        add_row(sources = "Ensembl 105",    known = "ensembl",   n = ensembl) %>%
        add_row(sources = "Kuo et al 2017", known = "Kuo et al", n = kuo) %>%
        mutate(sources = factor(sources, levels = c("Kuo et al 2017", "Ensembl 105", "GENE-SWitCH"))) %>%
        ggplot(aes(n, sources, fill = known)) +
            geom_col(position = "dodge") +
            geom_text(
                aes(label = n),
                hjust = 1.1,
                position = position_dodge(0.9)) +
            labs(
                title = paste0(specie, " - Total number of genes"),
                x = "Number of genes",
                y = "source")
}

ggsave(
    paste0("GENE_02-geom_col_-_", specie, "_-_FILTERED_-_Total_number_of_genes_ensembl.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


# Read support
cat("03 - Ploting Read support...\n")
d %>%
    select(
        tama_gene_id,   tama_transcript_id, ens_gene,
        ens_transcript, sources,    old_id, trans_read_count) %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    filter(trans_read_count >= 2) %>%
    group_by(tama_gene_id) %>%
    summarise(gene_read_count = sum(trans_read_count, na.rm = TRUE)) %>%
    mutate(cat =
        ifelse(gene_read_count == 1, "=1",
        ifelse(gene_read_count == 2, "=2",
        ifelse(between(gene_read_count, 3, 10), "[3,10]",
        ifelse(between(gene_read_count, 11, 100), "[11,100]", ">100"))))) %>%
    count(cat) %>%
    mutate(cat = factor(cat, levels = rev(c("=1", "=2", "[3,10]", "[11,100]", ">100")))) %>%
    ggplot(aes(n, cat)) +
        geom_col() +
        labs(
            title = paste0(specie, " - Number of reads supporting genes"),
            x = "Number of genes",
            y = "Number of reads") +
        geom_text(aes(label = n), hjust = 1.1)

ggsave(
    paste0("GENE_03-geom_col_-_", specie, "_-_FILTERED_-_Number_of_reads_supporting_genes.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


# Transcripts per gene
cat("04 - Ploting Transcripts per gene...\n")
d %>%
    filter(trans_read_count >= 2|is.na(trans_read_count)) %>%
    select(tama_gene_id, tama_transcript_id, sources) %>%
    mutate(sources = str_replace_all(sources, "ensembl,", "")) %>%
    group_by(tama_gene_id, sources) %>%
    count() %>%
    mutate(cat =
        ifelse(n == 1, "=1",
        ifelse(n == 2, "=2",
        ifelse(between(n, 3, 10), "[3,10]",
        ifelse(between(n, 11, 100), "[11,100]", ">100"))))) %>%
    ungroup() %>%
    count(sources, cat) %>%
    mutate(cat = factor(cat, levels = rev(c("=1", "=2", "[3,10]", "[11,100]", ">100")))) %>%
    group_by(sources) %>%
    mutate(n = n / sum(n) * 100) %>%
    ggplot(aes(n, sources, fill = cat)) +
        geom_col() +
        labs(
            title = paste0(specie, " - Number of transcripts per gene"),
            x = element_blank(),
            y = element_blank())

ggsave(
    paste0("GENE_04-geom_col_-_", specie, "_-_FILTERED_-_Number_of_transcripts_per_genes.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


# Gene per Dev. stage
cat("05 - Ploting Gene per Dev. stage...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_gene_id, source_line) %>%
    filter(!is.na(source_line)) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    mutate(
        source_line = factor(source_line,                                 levels = timepoints),
        Tissue      = factor(str_replace(source_line, match_stages, ""),  levels = tissues),
        Stage       = factor(str_replace(source_line, match_tissues, ""), levels = stages)) %>%
    select(tama_gene_id, Stage) %>%
    distinct() %>%
    count(Stage) %>%
    ggplot(aes(n, Stage)) +
        geom_col() +
        geom_text(aes(label = n), hjust = 1.1) +
        labs(
            title = paste0(specie, " - Number of genes per Dev. Stage"),
            x = "Number of genes",
            y = "Development stage")

ggsave(
    paste0("GENE_05-geom_col_-_", specie, "_-_FILTERED_-_Number_of_genes_per_dev_stage.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


# Gene per Tissue
cat("06 - Ploting Gene per Tissue...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_gene_id, source_line) %>%
    filter(!is.na(source_line)) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    mutate(
        source_line = factor(source_line,                                 levels = timepoints),
        Tissue      = factor(str_replace(source_line, match_stages, ""),  levels = tissues),
        Stage       = factor(str_replace(source_line, match_tissues, ""), levels = stages)) %>%
    select(tama_gene_id, Tissue) %>%
    distinct() %>%
    count(Tissue) %>%
    ggplot(aes(n, Tissue)) +
        geom_col() +
        geom_text(aes(label = n), hjust = 1.1) +
        labs(
            title = paste0(specie, " - Number of genes per Tissue"),
            x = "Number of genes",
            y = "Tissue")

ggsave(
    paste0("GENE_06-geom_col_-_", specie, "_-_FILTERED_-_Number_of_genes_per_tissue.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)


# Gene per Timepoint
cat("07 - Ploting Gene per Timepoint...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_gene_id, source_line) %>%
    filter(!is.na(source_line)) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    count(source_line) %>%
    mutate(
        source_line = factor(source_line, levels = timepoints),
        Tissue      = factor(str_replace(source_line, match_stages, ""), levels = tissues)) %>%
    ggplot(aes(n, source_line, fill = Tissue)) +
        geom_col() +
        geom_text(aes(label = n), hjust = 1.1) +
        labs(
            title = paste0(specie, " - Number of genes per Timepoint"),
            x = "Number of genes",
            y = "Timepoint") +
        theme(legend.position = "bottom") +
        guides(fill = guide_legend(ncol = 10))

ggsave(
    paste0("GENE_07-geom_col_-_", specie, "_-_FILTERED_-_Number_of_genes_per_timepoint.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


# Gene per Timepoint (Ensembl)
cat("08 - Ploting Gene per Timepoint (Ensembl)...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_gene_id, source_line, sources) %>%
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
            title = paste0(specie, " - Number of genes per Timepoint (Ensembl)"),
            x = "Number of genes",
            y = "Timepoints") +
        xlim(xlim_08)

ggsave(
    paste0("GENE_08-geom_col_-_", specie, "_-_FILTERED_-_Number_of_genes_per_timepoint_ensembl.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


# Number of exons per genes
# NB: Number of exons == Maximum number of exons in gene's transcripts
cat("09 - Ploting Number of exons per genes...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    filter(sources %in% c("isoseq", "ensembl,isoseq")) %>%
    select(tama_gene_id, blockCount) %>%
    distinct() %>%
    group_by(tama_gene_id) %>%
    top_n(n = 1, wt = blockCount) %>%
    ungroup() %>%
    mutate(blockCount = ifelse(blockCount >= 20, ">=20", blockCount)) %>%
    mutate(blockCount = factor(blockCount, levels = rev(c(seq(1, 100, 1), ">=20")))) %>%
    count(blockCount) %>%
    ggplot(aes(n, blockCount)) +
        geom_col() +
        labs(
            title = paste0(specie, " - Number of exons per gene"),
            x = "Number of genes",
            y = "Number of exons") +
        geom_text(aes(label = n), hjust = -0.1) +
        xlim(xlim_09)


ggsave(
    paste0("GENE_09-geom_col_-_", specie, "_-_FILTERED_-_Number_of_exons_per_genes.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


# Number of specific genes
cat("10 - Ploting Number of specific genes...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(tama_gene_id, source_line) %>%
    filter(!is.na(source_line)) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    group_by(tama_gene_id) %>%
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
            title = paste0(specie, " - Number of specific genes"),
            x = "Number of specific genes",
            y = "Timepoint") +
        theme(legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 10)) +
    xlim(xlim_10)



ggsave(
    paste0("GENE_10-geom_col_-_", specie, "_-_FILTERED_-_Number_of_specific_genes.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)

# q(status=0)
