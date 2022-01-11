# Get the unique(5,95) mutations based on GISAID
library(stringr)
args <- commandArgs(TRUE)

# Load and fix headers (~+- all get replaced by X)
mutations_file <- args[0]
mutations <- read.csv(mutations_file, stringsAsFactors = FALSE, header = FALSE)
names(mutations)[1:3] <- mutations[1, 1:3]

new_names <- str_replace(str_replace(str_replace(
    mutations[1,-c(1:3)], "\\+", "I"), "-", "D"), "~", "M")
names(mutations)[-c(1:3)] <- new_names

mutations <- mutations[-1, ]
for(i in 4:ncol(mutations)) 
    mutations[, i] <- as.numeric(mutations[, i])

# Get high frequency voi mutations
voi_muts <- list()
for(voi in mutations$lineage[mutations$variant_of_interest == 1]){
    mut_freq <- mutations[which(mutations$lineage == voi), -c(1:3)]
    voi_muts[[voi]] <- names(mut_freq)[mut_freq > 0.7]
}

# Get high frequency of other lineages
max_bg <- apply(mutations[mutations$variant_of_interest == 0, -c(1:3)], 2, max)
remove_from_voi <- names(max_bg)[max_bg > 0.5]

for(i in seq_along(voi_muts)){
    voi_muts[[i]] <- voi_muts[[i]][!voi_muts[[i]] %in% remove_from_voi]
}

voi_muts
