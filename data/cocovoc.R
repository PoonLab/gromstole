# COCOVOC = COunts and COverage for VOC
suppressPackageStartupMessages(library(dplyr))
library(tidyr)
library(here)

# Gather file names
mfiles <- list.files(here("results"), 
    full.names = TRUE, pattern="*.mapped.csv$", recursive=TRUE) %>%
    sort()
cfiles <- list.files(here("results"), 
    full.names = TRUE, pattern="*.coverage.csv$", recursive=TRUE) %>%
    sort()

# Gather mutations list (BIG)
bigvariantmatrix <- read.csv(here("working", "bigvariantmatrix.csv"), 
    row.names = 1)
bigmuts <- names(bigvariantmatrix)
bigpos <- sapply(strsplit(bigmuts, "[A-Za-z]"), paste0, collapse = "")
rm(bigvariantmatrix) # Just need the names


cat("\n")

# Subset mapped by mutations list, get coverage
for (i in 1:length(mfiles)) {
    # Re-initialize dataframe
    mut_df <- data.frame(
        mutation2 = bigmuts, 
        position2 = as.numeric(bigpos))

    # Deal with mappings
    mapped <- read.csv(mfiles[i], stringsAsFactors = FALSE)
    mapped$type <- substr(mapped$label, 1, 1)
    mapped$pos <- mapped$position + 1
    mapped$alt <- sapply(
        strsplit(mapped$label, as.character(mapped$pos)), 
            function(x) {
                gsub("^\\.", "", x[2])
    })
    mapped$mutation2 <- paste0(case_when(
            mapped$type == "~" ~ "m",
            mapped$type == "+" ~ "a",
            mapped$type == "-" ~ "d"
        ),
        mapped$pos, mapped$alt)
    mapped <- dplyr::select(mapped, mutation2, frequency, coverage)
    mut_df <- left_join(mut_df, mapped, by = c("mutation2"))

    # Deal with Coverage
    cover <- read.csv(cfiles[i], stringsAsFactors = FALSE)
    names(cover) <- c("position2", "coverage2")
    mut_df <- left_join(mut_df, cover, by = "position2")

    # Sanity checking
    if(length(bigpos) != nrow(mut_df)) print("Something went wrong")

    mfile <- paste0(strsplit(mfiles[i], "/")[[1]][2:3], collapse = "/")
    cfile <- paste0(strsplit(cfiles[i], "/")[[1]][2:3], collapse = "/")
    if (mfile != cfile) print(paste0("Directory mismatch: ", mfile, ", ", cfile))

    # Record lab info
    mut_df$lab <- ifelse(grepl("guelph", mfiles[i]), "guelph", 
        ifelse(grepl("waterloo", mfiles[i]), "waterloo", "western"))
    mut_df$sample <- strsplit(mfiles[i], "/") %>%
        sapply(function(x) x[length(x)]) %>%
        strsplit("\\.") %>%
        sapply(function(x) x[1])
    cat(mut_df$lab[1], "/", mut_df$sample[1], "\n")

    # Record results
    if(i == 1) {
        cocovoc <- mut_df
    } else {
        cocovoc <- rbind(cocovoc, mut_df)
    }
}

cocovoc$frequency[is.na(cocovoc$frequency) & (!is.na(cocovoc$coverage2)) & cocovoc$coverage2 > 0] <- 0
cocovoc$count <- cocovoc$frequency * cocovoc$coverage2


if(FALSE) { # Desperately need to parse the metadata
# load metadata
metas <- list.files("uploads", 
    full.names = TRUE, pattern="*meta*", recursive=TRUE)
m <- lapply(metas, function(x) {
  # Reading everything in as a character to avoid type issues (logical v. character)
  meta <- read.csv(x, colClasses = "character")
  meta
})
m2 <- dplyr::bind_rows(m)
m3 <- m2 %>%
    select(sample = Specimen.collector.sample.ID, 
        coldate = sample.collection.date, 
        location = geolocation.name..region.,
        longitude = geolocation.latitude,
        latitude = geolocation.longitude)
counts <- left_join(counts, m3, by = "sample")
}

write.csv(cocovoc, file = here("data", "cocovoc.csv"), row.names = FALSE)


