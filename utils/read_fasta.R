read_fasta <- function(file) {
    require(tibble)

    # Opening file in reading mode
    con  <- file(file, open = "r")

    # prepare vectors
    id <- c()
    seq <- c()
    seq_temp <- c()

    # parsing line by line
    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
        if (grepl('^>', line)) {
            id <- c(id, line)
            if (length(seq_temp) != 0) {
                seq = c(seq, seq_temp)
                seq_temp = c()
            }
        }
        else {
            seq_temp <- paste0(seq_temp, line)
        }
    }

    seq = c(seq, seq_temp)

    close(con)

    return(tibble(id = id, seq = seq))
}