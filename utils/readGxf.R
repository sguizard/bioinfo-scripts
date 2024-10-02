read_gff <- function(file) {
    require(readr)

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

read_gtf <- function(file, separate_attributes = FALSE) {
    require(readr)
    require(dplyr)

    gtf <- readr::read_tsv(
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
            sequence   = "c",
            source     = "c",
            feature    = "c",
            start      = "i",
            end        = "i",
            score      = "c",
            strand     = "c",
            phase      = "c",
            attributes = "c"),
            comment    = "#")

    if (isTRUE(separate_attributes)) {
        gtf <-
            gtf %>%
            dplyr::mutate(
                gene_id                      = str_match(attributes, 'gene_id "([^"]+)";')[,2],
                gene_version                 = str_match(attributes, 'gene_version "([^"]+)";')[,2],
                gene_name                    = str_match(attributes, 'gene_name "([^"]+)";')[,2],
                gene_source                  = str_match(attributes, 'gene_source "([^"]+)";')[,2],
                gene_biotype                 = str_match(attributes, 'gene_biotype "([^"]+)";')[,2],
                transcript_id                = str_match(attributes, 'transcript_id "([^"]+)";')[,2],
                transcript_version           = str_match(attributes, 'transcript_version "([^"]+)";')[,2],
                transcript_name              = str_match(attributes, 'transcript_name "([^"]+)";')[,2],
                transcript_source            = str_match(attributes, 'transcript_source "([^"]+)";')[,2],
                transcript_biotype           = str_match(attributes, 'transcript_biotype "([^"]+)";')[,2],
                exon_number                  = str_match(attributes, 'exon_number "([^"]+)";')[,2],
                exon_id                      = str_match(attributes, 'exon_id "([^"]+)";')[,2],
                exon_version                 = str_match(attributes, 'exon_version "([^"]+)";')[,2],
                protein_id                   = str_match(attributes, 'protein_id "([^"]+)";')[,2],
                protein_version              = str_match(attributes, 'protein_version "([^"]+)";')[,2],
                tag                          = str_match(attributes, 'tag "([^"]+)";')[,2],
                projection_parent_transcript = str_match(attributes, 'projection_parent_transcript "([^"]+)";')[,2]) %>%
            dplyr::select(-attributes)
    }
    return(gtf)
}


extractSEG <- function(file, out_file = NULL) {
    require(dplyr)
    require(readr)

    gtf <- 
        read_gtf(file) %>% 
        dplyr::mutate(
            transcript_id = str_match(attributes, 'transcript_id "([^"]+)";')[,2],
            gene_id       = str_match(attributes, 'gene_id "([^"]+)";')[,2])

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