library(here)

# Omicron mutations
ohmy <- read.csv(here("reports",
    "omicron-mutations.csv"))[, -1]
oh_diffs <- sapply(1:nrow(ohmy), function(x) paste(ohmy[x, 1:3], collapse = ""))

# Mutations from ALL sequences that have pangolin lineages
# Note: These come from old code that misses some, haven't done the new batch yet
pd <- read.csv(here("data", "pangodiffs.csv"))
pango_diffs <- strsplit(pd$diffs, split = "-")

# Initialize list
# oh_distinct[["0"]]: any omicron mutation that appears in NONE of the pangolin sequences
# oh_distinct[["1"]]: omicron mutations that only appear in ONE of the other sequences
oh_distinct <- list()

for(i in oh_diffs){
    pango_match <- sapply(pango_diffs, function(x) i %in% x)
    count_of_pango_matches <- as.character(sum(pango_match))
    oh_distinct[[count_of_pango_matches]] <- c(oh_distinct[[count_of_pango_matches]], i)
}

oh_distinct[order(as.numeric(names(oh_distinct)))]
