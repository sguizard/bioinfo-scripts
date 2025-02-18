require(readr)
require(dplyr)
require(tidyr)
require(purrr)
require(stringr)
require(furrr)

# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff.md
# GFF1, GFF2 http://www.sanger.ac.uk/resources/software/gff/spec.html
# GFF3 http://www.sequenceontology.org/gff3.shtml
# GVF http://www.sequenceontology.org/resources/gvf.html
# GTF http://mblab.wustl.edu/GTF22.html


#### read_gff ###########################################
read_gxf <- function(file, threads = 1) {
    cat(paste0("==> Reading ", file, "\n"))

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
    comment = "#",
    num_threads = threads)
}



#### separate_attributes ########################################
separate_attributes <- function(obj, threads = 1, nlines = 50000) {
    extract_att <- function(.data) {
        .data %>%
            slice_head(n = nlines) %>%
            separate_longer_delim(
                cols  = attributes,
                delim = stringr::regex(" ?; ?")) %>%
            filter(attributes != "") %>%
            separate_wider_delim(
                cols  = attributes,
                delim = stringr::regex("[ =]"),
                names = c("attributes", "value")) %>%
            distinct(attributes) %>%
            pull(attributes) %>%
            return()
    }

    add_att <- function(obj, name) {
        cat(paste0("==> Adding ", name, "\n"))
        pattern <- paste0(name, '[ =]"?([^"=;]+)"?')
        mutate(obj, "{name}" := str_match(attributes, pattern)[, 2])
    }

    apply_add_att <- function(x) {
        for (i in to_add) {
            x <- add_att(x, i)
        }

        x <-
            x %>%
            dplyr::select(-attributes) %>%
            purrr::discard(~all(is.na(.)))

        return(x)
    }

    if (threads == 1) {
        cat("==> Listing attributes\n")
        to_add <-
            obj %>%
            extract_att()

        obj <- apply_add_att(obj)
    }
    else {
        cat("==> Setup workers\n")
        plan(multisession, workers = threads)
        options(future.globals.maxSize = 10240 * 1024^2)

        cat(paste0("==> Spliting data in ", threads, " parts\n"))
        obj <-
            obj %>%
            group_by((row_number() - 1) %/% (n() / threads)) %>%
            nest %>%
            pull(data)

        cat("==> Listing attributes\n")
        to_add <-
            obj[[1]] %>%
            slice_head(n = nlines) %>%
            extract_att()

        cat("==> Adding attributes\n")
        obj <-
            future_map(
                obj,
                apply_add_att,
                .options = furrr_options(scheduling = 1)) %>%
            bind_rows()
    }

    return(obj)
}



#### read_gxf_and_separate ###########################################
read_gxf_and_separate <- function(file, threads = 1) {
    if (threads == 1) {
        gxf <- read_gxf(file = file)
        gxf <- separate_attributes(gxf)
    }
    else {
        gxf <- read_gxf           (file = file, threads = threads)
        gxf <- separate_attributes(obj  = gxf , threads = threads)
    }
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