library(here)
library(tidyr)
library(dplyr)
library(ggplot2)

cvr <- read.csv(here("results/cvr.csv"), stringsAsFactors = FALSE)

cvr$X <- sapply(strsplit(cvr$X, split = "/"), function(x) {
  strsplit(x[length(x)], split = "\\.")[[1]][1]
})
cvr <- rename(cvr, sample = X)

cvr_long <- pivot_longer(cvr, -sample, names_to = "mutation", values_to = "coverage")
head(cvr_long)

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
  pivot_longer(-c("sample", "lab", "coldate", "site"), names_to = "mutation", values_to = "count")

ggplot(counts_long, aes(x = mutation, y = sample, fill = log(count))) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  