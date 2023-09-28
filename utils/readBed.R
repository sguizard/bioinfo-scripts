read_bed6 <- function(file) {
    require(readr)
    read_tsv(
        file,
        col_names = c(
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand"),
        col_types = cols(
            chrom = col_character(),
            chromStart = col_integer(),
            chromEnd = col_integer(),
            name = col_character(),
            score = col_double(),
            strand = col_character()))
}

read_bed12 <- function(file, skip = 0) {
    require(readr)
    read_tsv(
        file,
        col_names = c(
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "blockCount",
            "blockSizes",
            "blockStarts"),
        col_types = cols(
            chrom = col_character(),
            chromStart = col_integer(),
            chromEnd = col_integer(),
            name = col_character(),
            score = col_double(),
            strand = col_character(),
            thickStart = col_integer(),
            thickEnd = col_integer(),
            itemRgb = col_character(),
            blockCount = col_integer(),
            blockSizes = col_character(),
            blockStarts = col_character()),
        skip = skip
    )
}