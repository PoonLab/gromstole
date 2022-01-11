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
sum(is.na(coco$coldate)); nrow(coco) 
sum(is.na(coco$location)); nrow(coco) 
sum(is.na(coco$latitude)); nrow(coco) 
# So no metadata, then.

# Prepare for the results
res_df <- coco %>% 
    select(sample, folder, lab, location) %>%
    distinct() %>%
    mutate(prop = rep(NA, n()), 
        lo = rep(NA, n()), hi = rep(NA, n()),
        status = rep(NA, n()),
        totol_succ = rep(NA, n()),
        totol_fail = rep(NA, n())
        )

for (i in seq_len(nrow(res_df))) {
    thiscoco <- coco[coco$sample == res_df$sample[i] & coco$lab == res_df$lab[i] & coco$folder == res_df$folder[i], ]
    success <- round(as.numeric(thiscoco$count))
    failure <- round(as.numeric(thiscoco$coverage - success))
    res_df$totol_fail[i] <- sum(failure, na.rm = TRUE)
    res_df$totol_succ[i] <- sum(success, na.rm = TRUE)

    errors <- c(sum(success, na.rm = TRUE) < 3, 
            sum(round(as.numeric(thiscoco$coverage)), na.rm = TRUE) < 10,
            any(success > round(as.numeric(thiscoco$count)), na.rm = TRUE))
    if (any(errors)) {
        print(c("<3 successes", "<10 coverage", "success > count")[which(errors)])
        print(res_df[i, c("sample", "lab", "folder")])
        res_df$status[i] <- "Insufficient data"
        next
    }

    cocomod <- tryCatch(glm(cbind(success, failure) ~ 1, family = binomial), error = function(e) e)
    if("error" %in% class(cocomod)) next
    res_df$prop[i] <- exp(cocomod$coef) / (1+exp(cocomod$coef))
    ci <- confint(cocomod)
    res_df$lo[i] <- exp(ci[1]) / (1+exp(ci[1]))
    res_df$hi[i] <- exp(ci[2]) / (1+exp(ci[2]))
}

res_df

pdf(file = here("results", "all-Langley.pdf"), width = 10, height = 10)
ggplot(res_df) + theme_bw() +
    aes(x = sample, y = prop, ymin = lo, ymax = hi) +
    geom_errorbar(aes(colour = lo > 0.01)) +
    scale_colour_manual(values = 1:2) +
    geom_point() +
    geom_hline(yintercept = 0.01, linetype = "dashed",
        colour = "grey") +
    facet_wrap(~lab, scales = "free_y") +
    coord_flip() +
    labs(colour = "CI above 0.01?")
dev.off()
