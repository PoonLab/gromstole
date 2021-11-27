library(here)
library(jsonlite)

ohmy <- read.csv(here("reports",
    "omicron-mutations.csv"))

head(ohmy)

pangodiffs <- readLines(here("data", "pangodiffs.json")) 

pd <- jsonlite::fromJSON(pangodiffs)

head(pd)
class(pd)

pd2 <- data.frame(
    qname = pd$qname, 
    lineage = pd$lineage, 
    diffs = sapply(pd$diffs, paste, collapse = "-"), 
    missing = sapply(pd$missing, paste, collapse = "-"),
    stringsAsFactors = FALSE)
str(pd2)

write.csv(pd2, file = "pangodiffs.csv", row.names = FALSE)
