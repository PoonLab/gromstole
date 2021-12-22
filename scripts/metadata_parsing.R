# Parsing Metadata
library(here)
library(lubridate)
library(dplyr)

metas <- list.files(here("uploads"), pattern = "metadata", 
    full.names = TRUE, recursive = TRUE)
metas <- metas[!grepl("template", metas)]
metas <- metas[!grepl("fixed", metas)]
metas

# Used google maps API to get this list; ask for the csv
locs <- read.csv(here("data", "maps", "sampling-locations.csv"))

# Put it all together
for(meta in metas) {
    metai <- read.csv(meta, stringsAsFactors = FALSE, colClasses = "character")
    names(metai) <- tolower(names(metai))
    metai$sample.collection.date <- parse_date_time(metai$sample.collection.date, c("mdy", "dmy", "ymd"))
    metai <- left_join(metai, locs, by = c("geolocation.name..region." = "loc"))
    metai$lab <- ifelse(grepl("guelph", meta), "guelph", ifelse(grepl("waterloo", meta), "waterloo", "western"))
    if(any(grepl("^x", names(metai)))) {
        metai <- metai[, -which(grepl("^x", names(metai)))]
    }
    metai$sample <- ifelse(length(metai$specimen.collector.sample.id) < 1, 
        metai$library.id, metai$specimen.collector.sample.id)

    if(meta == metas[1]) {
        all_meta <- metai
    } else {
        all_meta <- bind_rows(all_meta, metai)
    }
}

write.csv(all_meta, file = here("uploads", "metadata_fixed.csv"))

