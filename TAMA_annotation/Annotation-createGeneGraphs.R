#!/usr/bin/env Rscript

library(getopt)
library(tidyverse)
library(viridis)
library(this.path)
source(paste0(this.dir(), "/../utils/readGxf.R"))

`%nin%` <- Negate(`%in%`)

spec <- matrix(c(
    "help", "h", 0, "logical",
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
    xlim_02 <- c(0, 85000)
    xlim_08 <- c(0, 40000)
    xlim_09 <- c(0, 65000)
    xlim_10 <- c(0, 8000)

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
cat("00 - Loading data...\n")
d <- read_tsv(file = opt$data)

cat("00.1 - Counting genes...\n")
ensembl <-
    read_gtf(opt$ensembl) %>%
    filter(feature == "gene") %>%
    nrow()

cat("00.2 - Extract exons...\n")
ensembl_exon <-
    read_gtf(file = opt$ensembl, separate_attributes = TRUE) %>%
    filter(feature == "exon") %>%
    count(gene_id, transcript_id, name = "Exon") %>%
    group_by(gene_id) %>%
    slice_max(order_by = Exon, with_ties = FALSE) %>%
    ungroup() %>%
    count(Exon) %>%
    mutate(Project = "Ensembl", Source = "Ensembl 108") %>%
    select(Project, Source, Exon, n)

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



# Gene Number
cat("01 - Ploting Gene Number...\n")
d %>%
    select(
        merge_gene_id, merge_trans_id,
        sources_id_gene, trans_read_count) %>%
    filter(trans_read_count >= 2) %>%
    filter(sources_id_gene %in% c("ISOseq", "Ensembl,ISOseq")) %>%
    distinct(merge_gene_id) %>%
    count() %>%
    mutate(src = "GENE-SWitCH") %>%
    add_row(src = "Ensembl 108", n = ensembl) %>%
    mutate(src = factor(src, levels = c("Ensembl 108", "GENE-SWitCH"))) %>%
    ggplot(aes(n, src)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -.1) +
        scale_x_continuous(expand = expansion(mult = .15)) +
        labs(
            title = paste0(specie, " - Total number of genes"),
            x = "Number of genes",
            y = "")

ggsave(
    paste0(
        "GENE_01-geom_col_-_",
        specie,
        "_-_Total_number_of_genes.png"),
    units = "px",
    width = 755,
    height = 350,
    dpi = 125)



# Gene Number Ensembl
cat("02 - Ploting Gene Number Ensembl...\n")
p <- d %>%
    select(
        merge_gene_id,   merge_trans_id,
        Source = sources_id_gene, trans_read_count) %>%
    filter(Source %in% c("ISOseq", "Ensembl,ISOseq")) %>%
    distinct(merge_gene_id, Source) %>%
    count(Source, name = "Count") %>%
    add_row(Source = "Ensembl 108", Count = ensembl) %>%
    mutate(Project = ifelse(
        str_detect(Source, "ISOseq"),
        "GENE-SWitCH", "Ensembl")) %>%
    mutate(Project = factor(Project, levels = c("GENE-SWitCH", "Ensembl"))) %>%
    ggplot(aes(Count, Source)) +
        geom_col() +
        facet_grid(Project ~ ., scales = "free") +
        geom_text(aes(label = Count), hjust = -.1) +
        scale_x_continuous(expand = expansion(mult = .15)) +
    labs(
        title = paste0(specie, " - Total number of genes"),
        x = "Number of genes",
        y = "") +
    theme(legend.position = "bottom")

# if (specie == "Sscrofa11.1"){
#     p  + xlim(xlim_02)
# }

ggsave(
    paste0(
        "GENE_02-geom_col_-_",
        specie,
        "_-_Total_number_of_genes_ensembl.png"),
    plot = p,
    units = "px",
    width = 755,
    height = 350,
    dpi = 125)



# Gene Number Ensembl + Exon number
cat("03 - Ploting Gene Number Ensembl + Exon number...\n")
cat("03.1 - Ploting Gene Number Ensembl + Exon number (Count)...\n")
d %>%
    select(
        merge_gene_id,
        merge_trans_id,
        Source = sources_id_gene,
        trans_read_count,
        Exon = blockCount) %>%
    filter(Source %in% c("ISOseq", "Ensembl,ISOseq")) %>%
    distinct(merge_gene_id, Source, Exon) %>%
    group_by(merge_gene_id, Source) %>%
    slice_max(order_by = Exon) %>%
    ungroup() %>%
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
                " - Total number of genes + Number of Exons"),
            x = "Number of genes",
            y = "") +
        scale_fill_viridis(discrete = TRUE) +
        theme(legend.position = "bottom") +
        guides(fill = guide_legend(
            nrow = 1,
            reverse = TRUE,
            title.position = "top"))

ggsave(
    paste0(
        "GENE_03-geom_col_-_",
        specie,
        "_-_Total_number_of_genes_ensembl_plus_exons.png"),
    units = "px",
    width = 755,
    height = 450,
    dpi = 125)



# Gene Number Ensembl + Exon number (Percentages)
cat("03.2 - Ploting Gene Number Ensembl + Exon number (Percentage by Project)...\n")
d %>%
    select(
        merge_gene_id,
        merge_trans_id,
        Source = sources_id_gene,
        trans_read_count,
        Exon = blockCount) %>%
    filter(Source %in% c("ISOseq", "Ensembl,ISOseq")) %>%
    distinct(merge_gene_id, Source, Exon) %>%
    group_by(merge_gene_id, Source) %>%
    slice_max(order_by = Exon) %>%
    ungroup() %>%
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
            title = paste0(specie, " - Total number of genes"),
            x = "Number of genes + Exon Number",
            y = "") +
        scale_fill_viridis(discrete = TRUE) +
    theme(legend.position = "bottom") +
        guides(fill = guide_legend(
            nrow = 1,
            reverse = TRUE,
            title.position = "top"))


ggsave(
    paste0(
        "GENE_03-geom_col_-_",
        specie,
        "_-_Total_number_of_genes_ensembl_plus_exons_percentages_by_project.png"),
    units = "px",
    width = 755,
    height = 450,
    dpi = 125)


cat("03.3 - Ploting Gene Number Ensembl + Exon number (Percentage by source)...\n")
d %>%
    select(
        merge_gene_id,
        merge_trans_id,
        Source = sources_id_gene,
        trans_read_count,
        Exon = blockCount) %>%
    filter(Source %in% c("ISOseq", "Ensembl,ISOseq")) %>%
    distinct(merge_gene_id, Source, Exon)  %>%
    group_by(merge_gene_id, Source) %>%
    slice_max(order_by = Exon) %>%
    ungroup() %>%
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
            title = paste0(specie, " - Total number of genes"),
            x = "Number of genes + Exon Number",
            y = "") +
        scale_fill_viridis(discrete = TRUE) +
        theme(legend.position = "bottom") +
        guides(fill = guide_legend(
            nrow = 1,
            reverse = TRUE,
            title.position = "top"))


ggsave(
    paste0(
        "GENE_03-geom_col_-_",
        specie,
        "_-_Total_number_of_genes_ensembl_plus_exons_percentages_by_source.png"),
    units = "px",
    width = 755,
    height = 450,
    dpi = 125)



# Read support
cat("04 - Ploting Read support...\n")
d %>%
    select(
        merge_gene_id,   merge_trans_id,
        sources_id_gene, trans_read_count) %>%
    filter(sources_id_gene %in% c("ISOseq", "Ensembl,ISOseq")) %>%
    filter(trans_read_count >= 2) %>%
    group_by(merge_gene_id) %>%
    summarise(gene_read_count = sum(trans_read_count, na.rm = TRUE)) %>%
    mutate(cat =
        ifelse(gene_read_count == 1, "=1",
        ifelse(gene_read_count == 2, "=2",
        ifelse(between(gene_read_count, 3, 10), "[3,10]",
        ifelse(between(gene_read_count, 11, 100), "[11,100]", ">100"))))) %>%
    count(cat) %>%
    mutate(cat = factor(
        cat,
        levels = rev(c("=1", "=2", "[3,10]", "[11,100]", ">100")))) %>%
    ggplot(aes(n, cat)) +
        geom_col() +
        labs(
            title = paste0(specie, " - Number of reads supporting genes"),
            x = "Number of genes",
            y = "Number of reads") +
        geom_text(aes(label = n), hjust = -.1) +
        scale_x_continuous(expand = expansion(mult = .15))

ggsave(
    paste0(
        "GENE_04-geom_col_-_",
        specie,
        "_-_Number_of_reads_supporting_genes.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)



# Transcripts per gene
cat("05 - Ploting Transcripts per gene...\n")
d %>%
    select(
        merge_gene_id,
        merge_trans_id,
        Source = sources_id_gene,
        trans_read_count,
        Exon = blockCount) %>%
    filter(Source %in% c("ISOseq", "Ensembl,ISOseq")) %>%
    count(Source, merge_gene_id) %>%
    mutate(Transcript = case_when(
        n == 1 ~ "=1",
        n == 2 ~ "=2",
        between(n, 3, 10) ~ "[3,10]",
        between(n, 11, 100) ~ "[11,100]",
        .default = ">100")) %>%
    count(Source, Transcript) %>%
    group_by(Source) %>%
    mutate(p = n / sum(n) * 100) %>%
    mutate(Project = "GENE-SWitCH") %>%
    rbind(ensembl_transcripts_by_gene) %>%
    mutate(
        Project = factor(Project, levels = c("GENE-SWitCH","Ensembl")),
        Source = factor(
            Source,
            levels = c("Ensembl,ISOseq", "ISOseq", "Ensembl 108")),
        Transcript = factor(
            Transcript,
            levels = rev(c("=1", "=2", "[3,10]", "[11,100]", ">100")))) %>%
    ggplot(aes(p, Source, fill = Transcript)) +
        geom_col() +
        facet_grid(Project ~ ., scales = "free") +
        scale_x_continuous(expand = expansion(mult = .15)) +
        labs(
            title = paste0(specie, " - Number of transcripts per gene"),
            x = element_blank(),
            y = element_blank()) +
        scale_fill_viridis(discrete = TRUE)

ggsave(
    paste0(
        "GENE_05-geom_col_-_",
        specie,
        "_-_Number_of_transcripts_per_genes.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)



# Gene per Dev. stage
cat("06 - Ploting Gene per Dev. stage...\n")
d %>%
    select(merge_gene_id, source_line) %>%
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
    distinct(merge_gene_id, Stage) %>%
    count(Stage) %>%
    ggplot(aes(n, Stage)) +
        geom_col() +
        geom_text(aes(label = n), hjust = -.1) +
        labs(
            title = paste0(specie, " - Number of genes per Dev. Stage"),
            x = "Number of genes",
            y = "Development stage") +
        scale_x_continuous(expand = expansion(mult = .15))

ggsave(
    paste0(
        "GENE_06-geom_col_-_",
        specie,
        "_-_Number_of_genes_per_dev_stage.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)



# Gene per Tissue
cat("07 - Ploting Gene per Tissue...\n")
d %>%
    filter(trans_read_count >= 2) %>%
    select(merge_gene_id, source_line) %>%
    filter(!is.na(source_line)) %>%
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
    select(merge_gene_id, Tissue) %>%
    distinct() %>%
    count(Tissue) %>%
    mutate(Tissue = fct_reorder(Tissue, n)) %>%
    ggplot(aes(n, Tissue)) +
        geom_col() +
        geom_text(aes(label = n), hjust = 1.1) +
        labs(
            title = paste0(specie, " - Number of genes per Tissue"),
            x = "Number of genes",
            y = "Tissue")

ggsave(
    paste0(
        "GENE_07-geom_col_-_",
        specie,
        "_-_Number_of_genes_per_tissue.png"),
    units = "px",
    width = 755,
    height = 305,
    dpi = 125)



# Gene per Timepoint
cat("08 - Ploting Gene per Timepoint...\n")
d %>%
    select(merge_gene_id, source_line) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
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
        geom_text(aes(label = n), hjust = -.1) +
        labs(
            title = paste0(specie, " - Number of genes per Timepoint"),
            x = "Number of genes",
            y = "Timepoint") +
        scale_x_continuous(expand = expansion(mult = .20)) +
        guides(fill = guide_legend(reverse = TRUE))

ggsave(
    paste0(
        "GENE_08-geom_col_-_",
        specie,
        "_-_Number_of_genes_per_timepoint.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)


# Not relevant for genes, it's the transcripts that hold the exons
# # Number of exons per genes
# # NB: Number of exons == Maximum number of exons in gene's transcripts
# cat("09 - Ploting Number of exons per genes...\n")
# d %>%
#     filter(trans_read_count >= 2) %>%
#     filter(sources_id_gene %in% c("ISOseq", "Ensembl,ISOseq")) %>%
#     select(merge_gene_id, blockCount) %>%
#     distinct() %>%
#     group_by(merge_gene_id) %>%
#     top_n(n = 1, wt = blockCount) %>%
#     ungroup() %>%
#     mutate(blockCount = ifelse(blockCount >= 20, ">=20", blockCount)) %>%
#     mutate(blockCount = factor(blockCount, levels = rev(c(seq(1, 100, 1), ">=20")))) %>%
#     count(blockCount) %>%
#     ggplot(aes(n, blockCount)) +
#         geom_col() +
#         labs(
#             title = paste0(specie, " - Number of exons per gene"),
#             x = "Number of genes",
#             y = "Number of exons") +
#         geom_text(aes(label = n), hjust = -0.1) +
#         xlim(xlim_09)

# ggsave(
#     paste0("GENE_09-geom_col_-_", specie, "_-_Number_of_exons_per_genes.png"),
#     units = "px",
#     width = 755,
#     height = 610,
#     dpi = 125)



# Number of specific genes
cat("09 - Ploting Number of Tissue specific genes...\n")
d %>%
    select(merge_gene_id, source_line) %>%
    separate_rows(source_line, sep = ",") %>%
    mutate(source_line = str_replace(source_line, "_P[12]", "")) %>%
    distinct() %>%
    group_by(merge_gene_id) %>%
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
            title = paste0(specie, " - Number of specific genes"),
            x = "Number of specific genes",
            y = "Timepoint") +
        scale_x_continuous(expand = expansion(mult = .15)) +
        guides(fill = guide_legend(reverse = TRUE))

ggsave(
    paste0(
        "GENE_10-geom_col_-_",
        specie,
        "_-_Number_of_specific_genes.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)



# Number of Dev. Stage specific genes
cat("10 - Ploting Number of Dev. Stage specific genes...\n")
d %>%
    select(merge_gene_id, source_line) %>%
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
    distinct(merge_gene_id, Stage) %>%
    group_by(merge_gene_id) %>%
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
                " - Number of specific genes per Dev. Stage"),
            x = "Number of genes",
            y = "Dev. Stage") +
        scale_x_continuous(expand = expansion(mult = .15))

ggsave(
    paste0(
        "GENE_11-geom_col_-_",
        specie,
        "_-_Number_of_dev_stage_specific_genes.png"),
    units = "px",
    width = 755,
    height = 400,
    dpi = 125)



# Number of Tissue specific genes
cat("11 - Ploting Number of Tissue specific genes...\n")
d %>%
    select(merge_gene_id, source_line) %>%
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
    distinct(merge_gene_id, Tissue) %>%
    group_by(merge_gene_id) %>%
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
                " - Number of specific genes per Tissue"),
            x = "Number of genes",
            y = "Tissue") +
        scale_x_continuous(expand = expansion(mult = .15))

ggsave(
    paste0(
        "GENE_12-geom_col_-_",
        specie,
        "_-_Number_of_tissue_specific_genes.png"),
    units = "px",
    width = 755,
    height = 610,
    dpi = 125)