library(here)

# Omicron mutations
ohmy <- read.csv(here("reports",
    "omicron-mutations.csv"))[, -1]
oh_diffs <- sapply(1:nrow(ohmy), function(x) paste(ohmy[x, 1:3], collapse = ""))

# Mutations from ALL sequences that have pangolin lineages
# Note: These come from old code that misses some, haven't done the new batch yet
pd <- read.csv(here("data", "pangodiffs.csv"))
pango_diffs <- strsplit(pd$diffs, split = "-")

# Initialize dataframe
# pango_matches is the number of different sequences from pangodiffs
# that particular mutation appears in.
oh_distinct <- data.frame(mutation = oh_diffs, freq = ohmy$freq,
    pango_matches = 0)

for(i in oh_diffs){
    pango_match <- sapply(pango_diffs, function(x) i %in% x)
    count_of_pango_matches <- as.character(sum(pango_match))
    oh_distinct$pango_matches[oh_distinct$mutation == i] <- count_of_pango_matches
}

oh_distinct

write.csv(oh_distinct, file = here("data", "omicron_v_pango.csv"))
