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
      "attributes"
    ), 
    col_types = cols(
      sequence   = col_character(), 
      source     = col_character(), 
      feature    = col_character(), 
      start      = col_integer(), 
      end        = col_integer(), 
      score      = col_character(), 
      strand     = col_character(), 
      phase      = col_character(), 
      attributes = col_character()
    ), 
    comment = "#"
  )
}
