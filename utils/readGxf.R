read_gff <- function(file) {
    require(readr)
    read_tsv(
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
    gtf <- read_tsv(
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
            mutate(
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
            select(-attributes)
    }
    return(gtf)
}