library(rjags)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(here)
library(ggridges)

variantmat <- read.csv(here("data", "variantmatrix_med.csv"), row.names = 1)
variantmat <- as.matrix(variantmat)
dim(variantmat)
rankMatrix(variantmat) # Not full rank - some variants will be unidentifiable

cocovoc <- read.csv(here("data", "cocovoc_med.csv"), stringsAsFactors = FALSE)

cocovoc1 <- filter(cocovoc, sample == sample[1], lab == lab[1])
variantmat1 <- variantmat[, which(!is.na(cocovoc1$coverage2))]
cocovoc1 <- cocovoc1[!is.na(cocovoc1$coverage2), ]
cocovoc1$count <- round(cocovoc1$count)

p_model <- "
model{
    for(i in 1:N) {
        pvec_normalized[i] <- min(0.9999, pvec[i])
        count[i] ~ dbinom(pvec_normalized[i], coverage[i])
    }

    for(j in 1:P) {
        p1[j] ~ dbeta(1,1)
        p[j] <- 0.99*(p1[j] + 0.001)
    }
    pvec <- p %*% variantmat
}
"

t0 <- Sys.time()
coda1 <- jags.model(file = textConnection(p_model),
    data = list(count = cocovoc1$count, N = nrow(cocovoc1),
        coverage = cocovoc1$coverage2 + 1, # DANGER: Currently a hack to make it run
        P = nrow(variantmat), variantmat = variantmat),
    n.adapt = 500, n.chains = 3) %>%
    coda.samples(n.iter = 2000, variable.names = c("p"))
difftime(Sys.time(), t0, units = "mins")

summary(coda1)

codf <- bind_rows(lapply(1:length(coda1), 
    function(i) {
        x <- as.data.frame(coda1[[i]])
        x$chain <- i
        x$sample <- 1:nrow(x)
        x
    }
))

write.csv(codf, file = here("results", "cocovoccoda.csv"), row.names = FALSE)

voc_match <- data.frame(name = paste0("p.", 1:nrow(variantmat), "."), voc = rownames(variantmat))

codf_long <- tidyr::pivot_longer(codf, -c("chain", "sample")) %>%
    left_join(voc_match, by = "name")

ggplot(codf_long, aes(x = value, colour = voc)) + 
    geom_density() +
    xlim(c(0.005, 0.4))

codf_long %>%
    select(-chain, -sample, -name) %>%
    group_by(voc) %>%
    summarise(median = median(value)) %>%
    ggplot() +
        aes(y = median, x = voc, fill = voc) +
        theme_bw() +
        geom_col() + 
        theme(legend.position = "none") +
        coord_flip()

ggplot(codf_long) + 
    theme_bw() +
    aes(y = voc, x = value, colour = chain) +
    geom_density_ridges(scale = 2) +
    theme(legend.position = "none")
