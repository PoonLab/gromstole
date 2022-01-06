library(jsonlite)
library(here)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if("-q" %in% args) {
    q <- args[which(args == "-q") + 1]
} else {
    q <- 0.95
}

# All variants listed in the cov-lineages/constellations repo,
    # but we're going to use our own mutations list
VoC <- c('A.23.1', 'AV1', 'AY.4.2', 'AY.4', 'B.1.1.318', 'B.1.1.529', 
    'B.1.1.7', 'B.1.351', 'B.1.427', 'B.1.429', 'B.1.525', 'B.1.526', 
    'B.1.617.1', 'B.1.617.2', 'B.1.617.3', 'B.1.621', 'BA.1', 'BA.2', 
    'BA.3', 'C.37', 'P.1', 'P.2', 'P.3')

recent_sequences <- read.csv(here("data", "sequence_recency.csv"), 
    stringsAsFactors = FALSE)
VoC <- recent_sequences$lineage[order(-recent_sequences$last6)][1:66] # TODO: Sensitivity analysis
# Ensure that omicron is in the data.
VoC <- c("B.1.1.529", "BA.1", "BA.2", VoC)
VoC <- unique(VoC)

mutations <- read_json(here("data", "count-mutations.json"))
voc_muts <- unlist(sapply(VoC, function(x) 
    names(mutations)[which(grepl(x, names(mutations), fixed = TRUE))[1]]))
voc_muts <- unique(voc_muts[!is.na(voc_muts)])
mutations <- mutations[voc_muts]
mutations <- lapply(mutations, function(x) {
    x$mutations <- data.frame(mutation = names(x$mutations), number = as.numeric(x$mutation))
    row.names(x$mutations) <- NULL
    if (max(x$mutations$number) > 20) {
        x$mutations <- x$mutations[x$mutations$number > quantile(x$mutations$number, 0.9),]
    }
    x
})

print("Sanity Checks")
print("Number of columns")
length(unique(unlist(sapply(mutations, function(x) x$mutations$mutation))))
print("Variants with no mutations after removing bottom percentils")
which(sapply(mutations, function(x) nrow(x$mutations) == 0))

medvariantmatrix <- lapply(mutations, function(x) {
    muts <- x$mutations$mutation
    muts <- gsub("~", "m", gsub("-", "d", gsub("\\+", "a", gsub("\\|", "", x$mutations$mutation))))
    res <- as.data.frame(matrix(1, nrow = 1, ncol = length(muts)))
    names(res) <- muts
    res
}) %>% 
    bind_rows()

medvariantmatrix[is.na(medvariantmatrix)] <- 0
row.names(medvariantmatrix) <- names(mutations)

write.csv(medvariantmatrix, here("working", "variantmat_med.csv"))
