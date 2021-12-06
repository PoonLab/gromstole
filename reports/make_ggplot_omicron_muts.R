library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

cvr <- read.csv(here("results/cvr.csv"), stringsAsFactors = FALSE)

cvr$X <- sapply(strsplit(cvr$X, split = "/"), function(x) {
  strsplit(x[length(x)], split = "\\.")[[1]][1]
})
cvr <- rename(cvr, sample = X)

cvr_long <- pivot_longer(cvr, -sample, names_to = "mutation", values_to = "coverage")

ggplot(cvr_long, aes(x = mutation, y = sample, fill = log(coverage))) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


counts <- read.csv(here("results", "counts.csv"), stringsAsFactors = FALSE)

names(counts)
counts$sample <- sapply(strsplit(counts$sample, split = "/"), function(x) {
  strsplit(x[length(x)], split = "\\.")[[1]][1]
})

str(counts)

counts_long <- counts %>%
  select(-X) %>% 
  pivot_longer(-c("sample", "lab", "coldate", "site"), 
    names_to = "mutation", values_to = "count") 

ggplot(counts_long, aes(x = mutation, y = sample, fill = log(count))) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


frequency <- left_join(counts_long, cvr_long, by = c("sample", "mutation")) %>%
  mutate(frequency = ifelse((!is.na(count)) & coverage > 0, count / coverage, 
    ifelse(coverage > 0, 0, NA)))
dim(frequency)

ggplot(frequency, aes(x = mutation, y = sample, fill = log(frequency))) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
