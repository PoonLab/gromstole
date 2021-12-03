library(dplyr) # only for joining small dfs later - not for large dfs!
library(lubridate)

# Requires the unzipped metadata because I'm not good at zipped files.
metadata <- read.csv("metadata_2021-12-02.tsv", sep = "\t", stringsAsFactors=FALSE)
metadate <- metadata[, c("date", "pango_lineage", "country_exposure")]
rm(metadata) # File is big
metadate$ymd <- ymd(metadate$date)

for (month_lag in 1:12) {
    # Find all lineages more recent than month_lag
    index <- metadate$ymd > Sys.time() - dmonths(month_lag)
    # Subset and make a table
    last_few_months <- table(metadate$pango_lineage[index]) %>%
        as.data.frame()
    names(last_few_months) <- c("lineage", 
        paste0("last", month_lag, collapse = ""))
    if (!exists("recency")) {
        recency <- last_few_months
    } else {
        recency <- full_join(recency, last_few_months, by = "lineage")
    }
}

recency[is.na(recency)] <- 0
recency <- recency[order(-recency$last12), ]

write.csv(recency, file = "sequence_recency.csv")
