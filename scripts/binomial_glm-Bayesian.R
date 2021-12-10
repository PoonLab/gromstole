# Re-creating Art's binomial GLMs in a Bayesian framework.
# Why? Because more complicated models will come later.
    # Spoliers: Split p into p_BA1 and p_BA2 and p_BA1 + p_BA2

library(here)
library(rjags)
library(dplyr)
library(ggplot2)

p_model <- "
model{
    for(i in 1:N) {
        count[i] ~ dbinom(p, coverage[i])
    }

    p ~ dbeta(1,1)
    # spike and slab prior for p?
    #t0 ~ dbinom(0.5, 1)
    #t01 <- t0 + (1 - t0)*0.0001
    #p <- p0*t0
}
"

# run coco.R on langley to get coco.csv
coco <- read.csv(here("data", "coco.csv"), stringsAsFactors = FALSE)

head(coco)

coco1 <- dplyr::filter(coco, sample == "GB2point2")
head(coco1)

cococoda <- jags.model(file = textConnection(p_model), 
    data = list(count = coco1$count, coverage = coco1$coverage, N = nrow(coco1)),
    n.chains = 3, n.adapt = 500) %>%
    coda.samples(variable.names = c("p"), n.iter = 1000)

codadf <- bind_rows(lapply(1:length(cococoda), function(i) {
    x <- as.data.frame(cococoda[[i]])
    x$chain <- i
    x$sample <- 1:nrow(x)
    x
}))

ggplot(codadf) +
    aes(x = sample, y = p, colour = factor(chain)) + 
    geom_line()


resdf <- data.frame(sample = coco$sample, lab = coco$lab, 
        phat = NA, plo = NA, phi = NA) %>%
    distinct()

for(i in seq_len(nrow(resdf))){
    cocoi <- coco[coco$sample == resdf$sample[i], ]
    if(sum(is.na(cocoi$count)) > 6) next

    cococoda <- jags.model(
        file = textConnection(p_model), 
        data = list(
            count = cocoi$count, 
            coverage = cocoi$coverage, 
            N = nrow(cocoi)),
        n.chains = 3, 
        n.adapt = 1000,
        quiet = TRUE) %>%
    coda.samples(
        variable.names = c("p"), 
        n.iter = 2000,
        progress.bar = "none")

    cocodf <- bind_rows(lapply(1:length(cococoda), 
        function(i) {
            x <- as.data.frame(cococoda[[i]])
            x$chain <- i
            x$sample <- 1:nrow(x)
            x
        }
    ))

    resdf$phat[i] <- mean(cocodf$p)
    resdf$plo[i] <- quantile(cocodf$p, 0.05)
    resdf$phi[i] <- quantile(cocodf$p, 0.95)
}

head(resdf)

ggplot(data = filter(resdf, !is.na(phat))) +
    aes(x = sample, y = phat, ymin = plo, ymax = phi) + 
    geom_point(size = 2) + 
    geom_errorbar() +
    facet_wrap(~ lab, scales = "free_y") +
    coord_flip()
