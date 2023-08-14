read_blast7 <- function(file) {
    require(tidyverse)
    cmd    <- paste0("grep -e 'Fields:' ", file, "|head -n 1|sed '1s/# Fields: //'")
    header <- (system(cmd, intern = TRUE) %>% str_split(pattern = ", "))[[1]]

    header <- header %>% str_replace("query acc.", "qid")
    header <- header %>% str_replace("subject acc.", "sid")
    header <- header %>% str_replace("% identity", "pident")
    header <- header %>% str_replace("alignment length", "alength")
    header <- header %>% str_replace("mismatches", "mm")
    header <- header %>% str_replace("gap opens", "gaps")
    header <- header %>% str_replace("q. start", "qstart")
    header <- header %>% str_replace("q. end", "qend")
    header <- header %>% str_replace("s. start", "sstart")
    header <- header %>% str_replace("s. end", "send")
    header <- header %>% str_replace("bit score", "bitscore")
    header <- header %>% str_replace("query length", "qlen")
    header <- header %>% str_replace("% query coverage per subject", "qcovs")
    header <- header %>% str_replace("% query coverage per hsp", "qcovhsp")
    header <- header %>% str_replace("% query coverage per uniq subject", "qcovus")
    header <- header %>% str_replace("subject title", "stitle")

    read_tsv(
        file,
        col_names = header,
        comment   = "#",
        col_types = c(
            alength = "i",
            mm      = "i",
            gaps    = "i",
            qstart  = "i",
            qend    = "i",
            sstart  = "i",
            send    = "i")) %>%
        return
}
