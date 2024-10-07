require(readr)
require(dplyr)
require(tidyr)
require(purrr)
require(stringr)

# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff.md

#### read_gff ###########################################
read_gxf <- function(file) {
    readr::read_tsv(
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
            score      = col_character(),
            strand     = col_character(),
            phase      = col_character(),
            attributes = col_character()),
    comment = "#")
}



#### separate_attributes ########################################
separate_attributes <- function(obj) {
    add_att <- function(obj, name) {
        cat(paste0("==> Adding ", name, "\n"))
        pattern <- paste0(name, '[ =]"?([^"=;]+)"?')
        mutate(obj, "{name}" := str_match(attributes, pattern)[, 2])
    }

    cat("==> Listing attributes\n")
    to_add <-
        obj %>%
        separate_longer_delim(
            cols  = attributes,
            delim = stringr::regex(" ?; ?")) %>%
        filter(attributes != "") %>%
        separate_wider_delim(
            cols  = attributes,
            delim = stringr::regex("[ =]"),
            names = c("attributes", "value")) %>%
        distinct(attributes) %>%
        pull(attributes)

    for (i in to_add) {
        obj <- add_att(obj, i)
    }

    obj <-
        obj %>%
        dplyr::select(-attributes) %>%
        purrr::discard(~all(is.na(.)))

    return(obj)

}



#### read_gxf_and_separate ###########################################
read_gxf_and_separate <- function(file) {
    cat(paste0("==> Reading ", file, "\n"))
    gxf <- read_gxf(file = file)
    gxf <- separate_attributes(gxf)

    return(gxf)
}



#### extract_seg ###########################################
extract_seg <- function(file, out_file = NULL) {
    gtf <-
        read_gtf(file) %>%
        dplyr::mutate(
            transcript_id = str_match(attributes, 'transcript_id "([^"]+)";')[, 2],
            gene_id       = str_match(attributes, 'gene_id "([^"]+)";')[, 2])

    # List Single Exon Transcript ids
    set <-
        gtf %>%
        dplyr::filter(feature == "exon") %>%
        dplyr::count(transcript_id) %>%
        dplyr::filter(n == 1) %>%
        dplyr::pull(transcript_id)

    # List Single Exon Gene ids
    seg <-
        gtf %>%
        dplyr::filter(transcript_id %in% set) %>%
        dplyr::distinct(gene_id) %>%
        dplyr::pull(gene_id)

    # Extract rows Single Exon Transcripts and their associated gene line
    seg_gtf <-
        gtf %>%
        dplyr::filter(
            transcript_id %in% set | 
            (gene_id      %in% seg & is.na(transcript_id))) %>%
        dplyr::select(-transcript_id, -gene_id)

    if (!is.null(out_file)) {
        readr::write_tsv(seg_gtf, out_file, col_names = FALSE, escape = "none")
    }

    return(seg_gtf)
}