library(here)
library(ggplot2)
suppressPackageStartupMessages(library(dplyr))

# Open coco (first co: count, second co: coverage)
# Processes all of the mapped and coverage files on langley,
    # only includes omicron-specific mutations (currently)
# ssh dbecker7@langley "cd gromstole-aux; Rscript coco.R"
# scp dbecker7@langley:gromstole-aux/coco.csv ./data
coco <- read.csv(here("data", "coco.csv"))

head(coco)
table(coco$coldate)
sum(is.na(coco$coldate)); nrow(coco)  # Well, damn.
sum(is.na(coco$location)); nrow(coco) # Well, damn.
sum(is.na(coco$latitude)); nrow(coco) # Well, damn.
# So no metadata, then.

# Prepare for the results
res_df <- coco %>% 
    select(sample, lab, coldate, location, longitude, latitude) %>%
    distinct() %>%
    mutate(prop = rep(NA, n()), 
        lo = rep(NA, n()), hi = rep(NA, n()),
        status = rep(NA, n()))

for (i in seq_len(nrow(res_df))) {
    thiscoco <- coco[coco$sample == res_df$sample[i], ]
    success <- thiscoco$count
    failure <- thiscoco$coverage - success

    if (sum(success, na.rm = TRUE) < 3 | 
            sum(failure, na.rm = TRUE) < 10) {
        res_df$status[i] <- "Insufficient data"
        next
    }

    cocomod <- glm(cbind(success, failure) ~ 1, family = binomial)
    res_df$prop[i] <- exp(cocomod$coef) / (1+exp(cocomod$coef))
    ci <- confint(cocomod)
    res_df$lo[i] <- exp(ci[1]) / (1+exp(ci[1]))
    res_df$hi[i] <- exp(ci[2]) / (1+exp(ci[2]))
}

res_df

pdf(file = here("results", "all-Langley.pdf"), width = 8, height = 10)
ggplot(res_df) + theme_bw() +
    aes(x = sample, y = prop, ymin = lo, ymax = hi) +
    geom_errorbar() +
    geom_point() +
    geom_hline(yintercept = 0.01, linetype = "dashed",
        colour = "grey") +
    facet_wrap(~lab, scales = "free_y") +
    coord_flip()
dev.off()
