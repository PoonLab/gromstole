library(here)
library(jsonlite)

pango_rep <- read_json(here("data", "representative_pangodiffs.json"))

# A horribly inefficient way to do this
pangodf <- data.frame(
    qname = sapply(pango_rep, function(x) x$qname),
    lineage = sapply(pango_rep, function(x) x$lineage),
    diffs = sapply(pango_rep, function(x) paste(x$diffs, collapse = "|")),
    missing = sapply(pango_rep, function(x) paste(unlist(x$missing), collapse = "|"))
    )
rm(pango_rep)

head(pangodf)

write.csv(pangodf, file = here("data", "representative_pangodiffs.csv"), row.names = FALSE)
