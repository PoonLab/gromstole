library(rjags)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(here)
library(ggridges)

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

variantmat <- read.csv(here("data", "variantmat_med.csv"), row.names = 1)
variantmat <- as.matrix(variantmat)
dim(variantmat)
rankMatrix(variantmat) # Not full rank - some variants will be unidentifiable

cocovoc <- read.csv(here("data", "cocovoc_med.csv"), stringsAsFactors = FALSE)

labs <- unique(cocovoc$lab)
for(lab in labs) {
    samples <- unique(cocovoc$sample[cocovoc$lab == lab])

    for(sample in samples){
        print(paste0(lab, ", ", sample))
        cocovoc1 <- cocovoc[cocovoc$sample == sample & cocovoc$lab == lab,]
        variantmat1 <- variantmat[, which(!is.na(cocovoc1$coverage2))]
        cocovoc1 <- cocovoc1[!is.na(cocovoc1$coverage2), ]
        cocovoc1$count <- round(cocovoc1$count)

        

        t0 <- Sys.time()
        coda1 <- jags.model(file = textConnection(p_model),
            data = list(count = cocovoc1$count, N = nrow(cocovoc1),
                coverage = cocovoc1$coverage2 + 1, # DANGER: Currently a hack to make it run
                P = nrow(variantmat), variantmat = variantmat),
            n.adapt = 500, n.chains = 2) %>%
            coda.samples(n.iter = 3000, variable.names = c("p"))
        print(difftime(Sys.time(), t0, units = "mins"))

        codf <- bind_rows(lapply(1:length(coda1), 
            function(i) {
                x <- as.data.frame(coda1[[i]])
                x$chain <- i
                x$sample <- 1:nrow(x)
                x
            }
        ))
        voc_match <- data.frame(name = paste0("p.", 1:nrow(variantmat), "."), voc = rownames(variantmat))
        codf_long <- tidyr::pivot_longer(codf, -c("chain", "sample")) %>%
            left_join(voc_match, by = "name")

        write.csv(
            x = codf_long, 
            file = here("results", paste0(lab, "-",sample, "-", "cocovoccoda.csv")), 
            row.names = FALSE
        )

    }
}


