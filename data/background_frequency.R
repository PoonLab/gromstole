library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(here)

# Frequencies of each mutation in GISAID
freqs <- read.csv(here("data", "compare-out.csv"), 
    stringsAsFactors = FALSE, header = FALSE)
# Header row has special characters, substituting ~-+ with MDI
names(freqs) <- str_replace_all(freqs[1, ], fixed("~"), "M") %>%
    str_replace_all(fixed("-"), "D") %>%
    str_replace_all(fixed("+"), "I")
# remove existing row (headers)
freqs <- freqs[-1, ]
# Make it all numeric again (except lineage column)
freaks <- mutate_at(freqs, -1, as.numeric)
dim(freaks)

# Frequency of sequences within last n months
recency <- read.csv(here("data", "sequence_recency.csv"), 
    stringsAsFactors = FALSE)
last6 <- filter(recency, last6 > 10) %>% 
    select(lineage, last6) %>%
    mutate(pop_freq = last6/sum(last6)) %>%
    arrange(lineage)
head(last6)

# Frequency of sequenced omicron mutations
# Changing label to match column names of freaks
omicron <- read.csv(here("data", "omicron-mutations.csv")) %>%
    mutate(type2 = case_when(
        type == "~" ~ "M", 
        type == "-" ~ "D", 
        type == "+" ~ "I")) %>%
    unite(labels, type2, pos, alt, sep = "|") %>% 
    select(-type)
head(omicron)

# Frequency of Omicron mutations in the last 6 months
freak2 <- freaks[, c("lineage", "count", omicron$label)] %>%
    left_join(last6, by = "lineage") %>%
    select(lineage, pop_freq, everything()) %>%
    select(-last6) %>% 
    filter(!is.na(pop_freq))
freak2[1:10, 1:10]

freak3 <- as.matrix(select(freak2, -lineage, -pop_freq, -count))
rownames(freak3) <- freak2$lineage

dim(freak3)
length(freak2$pop_freq)
background_frequency <- data.frame(
        freq = t(freak2$pop_freq %*% freak3)) %>%
    rownames_to_column("label")


write.csv(background_frequency,
    file = here("data", "background_frequency.csv"),
    row.names = FALSE)

