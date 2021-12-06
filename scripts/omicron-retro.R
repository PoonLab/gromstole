library(here)

# load mutations specific to Omicron variant
omicron <- read.csv(here('working/omicron/omicron-specific.csv'), 
  row.names = 1, stringsAsFactors = FALSE)
omicron$X <- NULL  # remove duplicate row names
omicron$label <- paste(omicron$type, omicron$pos, omicron$alt, sep="|")


mfiles <- list.files(here("results"), full.names = TRUE,
                     pattern="*.mapped.csv$", recursive=TRUE)
cfiles <- list.files(here("results"), full.names = TRUE,
                     pattern="*.coverage.csv$", recursive=TRUE)

## generate summary stats
#cover <- sapply(cfiles, function(f) {
#  coverage <- read.csv(f)
#  #coverage$coverage[coverage$coverage==0] <- 1
#  sum(coverage$coverage > 0)
#})

# determine coverage at omicron sites
cover <- sapply(cfiles, function(f) {
  coverage <- read.csv(f, stringsAsFactors = FALSE)
  # NOTE: coverage$position is 0-index, omicron$pos is 1-index
  coverage$coverage[which(is.element(coverage$position+1, omicron$pos))]
})
row.names(cover) <- omicron$label



maps <- sapply(mfiles, function(f) {
  mapped <- read.csv(f, stringsAsFactors = FALSE)
  mapped$type <- substr(mapped$label, 1, 1)
  mapped$pos <- mapped$position + 1
  mapped$alt <- sapply(strsplit(mapped$label, as.character(mapped$pos)), 
                       function(x) {
                         gsub("^\\.", "", x[2])
                       })
  # convert to omicron label type
  mapped$label <- paste(mapped$type, mapped$pos, mapped$alt, sep='|')
  
  index <- match(omicron$label, mapped$label)
  mapped$frequency[index]
})
row.names(maps) <- omicron$label

maps <- as.data.frame(maps)

# set column names to sample identifiers
names(maps) <- sapply(strsplit(colnames(maps), split = "/"), function(x) {
  strsplit(x[length(x)], split = "\\.")[[1]][1]
})

# set zeroes for sites with non-zero coverage
cover <- as.data.frame(cover)
mask <- (cover > 0)
for (i in 1:nrow(maps)) {
  for (j in 1:ncol(maps)) {
    if (is.na(maps[i,j])) {
      if (mask[i,j]) {
        maps[i,j] <- 0
      }
    }
  }
}

counts <- as.data.frame(t(maps*cover))
names(counts) <- ifelse(omicron$mutation=='None', 
                        gsub("~\\|", "", omicron$label), 
                        omicron$mutation)
counts$sample <- names(maps)
counts$lab <- sapply(cfiles, function(x) strsplit(x, "/")[[1]][1])
counts$coldate <- NA
counts$site <- NA

# load metadata
metas <- list.files(here("uploads"), full.names = TRUE,
                    pattern="*meta*", recursive=TRUE)
m <- lapply(metas, function(x) {
  # Reading everything in as a character to avoid type issues (logical v. character)
  meta <- read.csv(x, colClasses = "character")
  meta
})
#m <- do.call(rbind, m)
m <- dplyr::bind_rows(m)
idx <- match(counts$sample, m$Specimen.collector.sample.ID)
counts$coldate[!is.na(idx)] <- as.Date(
  m$sample.collection.date[idx[!is.na(idx)]], 
  format="%m/%d/%Y"
  )
counts$coldate <- as.Date(as.integer(counts$coldate), origin='1970-01-01')

write.csv(counts, here("results/counts.csv"))

# make a nicer coverage file too

cvr <- as.data.frame(t(cover))
names(cvr) <- ifelse(omicron$mutation=='None', 
                     gsub("~\\|", "", omicron$label), 
                     omicron$mutation)
row.names(cvr) <- counts$sample

write.csv(cvr, here("results/cvr.csv"))

# ===============================
# try out binomial regression test

y <- as.integer(counts[40, 1:16])
n <- as.integer(cvr[40, 1:16])
fit <- glm(cbind(y, n) ~ 1, family='binomial')
exp(fit$coef) / (1+exp(fit$coef))  # probability

probs <- c()
lo <- c()
hi <- c()
for (i in 1:nrow(counts)) {
  y <- as.integer(counts[i, 1:16])
  n <- as.integer(cvr[i, 1:16])
  p <- tryCatch({
    fit <- glm(cbind(y, n) ~ 1, family='binomial')
    exp(fit$coef) / (1+exp(fit$coef))  # probability
    }, 
    error = function(cond) { 
      return (NA) 
    })
  
  probs <- c(probs, p)
  if (is.na(p)) {
    lo <- c(lo, NA)
    hi <- c(hi, NA)
  } else {
    ci <- confint(fit)
    lo <- c(lo, exp(ci[1]) / (1+exp(ci[1])))
    hi <- c(hi, exp(ci[2]) / (1+exp(ci[2])))
  }
}
par(mar = c(3,5,1,1))
barplot(as.numeric(probs)[1:80], horiz=T, 
        names.arg=gsub("\\.[a-z0-9]+$", "", counts$sample[1:80]), 
        las=1, cex.names=0.5)
segments(x0=lo[1:80], x1=hi[1:80], y0=(1:80-0.5)*1.2)

# ============================

# make a nice plot


# ============================

pdf(here("results/detection.pdf"), width=6, height=6)

par(mar=c(5,7,3,1))
plot(NA, xlim=c(1e-5, 0.1), ylim=c(0.5, nrow(omicron)+0.5), 
     log='x', xaxt='n', yaxt='n', xlab='Frequency', ylab='')
abline(v=c(1e-4, 1e-3, 0.01, 0.1), col='grey80')
title(ylab="Mutation", line=5)
#axis(side=1, at=c(1e-5, 1e-4, 1e-3, 0.01, 0.1), 
#     label=c(0, expression(10^-4), expression(10^-3), "0.01", '0.1'))
axis(side=1, at=c(1e-5, 1e-4, 1e-3, 0.01, 0.1),
     label=c("Zero", "0.01%", "0.1%", "1%", "10%"))
axis(side=2, at=1:nrow(omicron), 
     label=ifelse(omicron$mutation=="None", omicron$label, omicron$mutation), 
     las=1, cex.axis=0.6)
abline(v=2.5e-5, lwd=1.5)

col <- as.integer(as.factor(
  sapply(mfiles, function(x) strsplit(x, "/")[[1]][1])
  ))
pal <- c(
  "#C20430",  # guelph
  "#E4B429",  # waterloo
  "#4F2683"  # western
  )

for (i in 1:nrow(omicron)) {
  x <- as.numeric(maps[i, ])
  z <- as.numeric(cover[i,])
  size <- x*z
  # only claim zero if we have coverage
  x[z > 10 & is.na(x)] <- 1e-5
  #print(x)
  #size <- sapply(log10(z), function(xx) max(0, xx))
  
  points(x=x, y=rep(i, length(x)), cex=ifelse(is.na(size), log10(z), sqrt(size)), 
         col=pal[col], pch=ifelse(x==1e-5, 0, 1), lwd=ifelse(x==1e-5, 1, 2))
  text(x=1.6e-5, y=i, label=sum(x==1e-5, na.rm=T), adj=0, cex=0.6)
}

par(xpd=NA)
points(x=c(5e-5, 1.5e-4, 5e-4), y=rep(18, 3), 
       cex=sqrt(c(1, 5, 20)))
text(x=c(6e-5, 2e-4, 8e-4), y=rep(18, 3), 
     label=c(1, 5, "20 reads"), cex=0.7, adj=0)
par(xpd=FALSE)

dev.off()



guelph <- as.data.frame(maps[,grepl("guelph", colnames(maps))])
names(guelph) <- gsub("guelph/.+/(.+)\\.mapped\\.csv", "\\1", names(guelph))

# modify this data frame so that NA mutations with decent coverage are zero
mask <- as.data.frame(cover[,grepl("guelph", colnames(cover))])
mask <- (mask > 0)
for (i in 1:nrow(guelph)) {
  for (j in 1:ncol(guelph)) {
    if (is.na(guelph[i,j])) {
      if (mask[i,j]) {
        guelph[i,j] <- 0
      }
    }
  }
}

gcounts <- guelph * cover[,grepl("guelph", colnames(cover))]

write.csv(gcounts, here("results/guelph-counts.csv"))


meta <- read.csv(here("uploads/guelph/20211105_151243/metadata_guelph_20211105_151243.csv"))

