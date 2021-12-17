library(here)

# load mutations specific to Omicron variant
omicron <- read.csv(here('data/omicron-specific.csv'), 
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
row.names(counts) <- NULL
names(counts) <- ifelse(omicron$mutation=='None', 
                        gsub("~\\|", "", omicron$label), 
                        omicron$mutation)
counts$sample <- names(maps)
# switching to `here` library broke this line
counts$lab <- sapply(cfiles, function(x) {
  tokens <- strsplit(x, "/")[[1]]
  tokens[length(tokens)-2]
  })
counts$coldate <- NA
counts$site <- NA

# load metadata
#metas <- list.files(here("uploads"), full.names = TRUE,
                    pattern="*meta*", recursive=TRUE)
#m <- lapply(metas, function(x) {
  # Reading everything in as a character to avoid type issues (logical v. character)
#  meta <- read.csv(x, colClasses = "character")
#  meta
#})
m <- read.csv("uploads/waterloo/run6/metadata.csv", header=T)

#m <- do.call(rbind, m)
m <- dplyr::bind_rows(m)
idx <- match(counts$sample, m$specimen.collector.sample.ID)  # they changed case!
counts$coldate[!is.na(idx)] <- as.Date(
  m$sample.collection.date[idx[!is.na(idx)]], 
  format="%b %d %Y"  #"%m/%d/%Y"
  )
counts$coldate <- as.Date(as.integer(counts$coldate), origin='1970-01-01')
counts$site <- m$geolocation.name..region.[idx]
#write.csv(counts, here("results/counts.csv"))

# make a nicer coverage file too

cvr <- as.data.frame(t(cover))
names(cvr) <- names(counts)[1:ncol(cvr)]
#row.names(cvr) <- counts$sample

write.csv(cvr, here("results/cvr.csv"))

# ===============================
# try out binomial regression test
b529 <- which(fin$lineage == 'B.1.1.529')
ba1 <- which(fin$lineage=='BA.1' | fin$lineage=='B.1.1.529')
ba2 <- which(fin$lineage=='BA.2' | fin$lineage=='B.1.1.529')

# sort by site and collection date
idx <- order(counts$site, counts$coldate)
counts <- counts[idx, ]
cvr <- cvr[idx, ]

est.freq <- function(midx) {
  probs <- c()
  lo <- c()
  hi <- c()
  for (i in 1:nrow(counts)) {
    y <- as.integer(counts[i, midx])  # number of "successes"
    n <- as.integer(cvr[i, midx])  # number of trials
    p <- tryCatch({
      fit <- glm(cbind(y, n-y) ~ 1, family='binomial')
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
      suppressMessages(ci <- confint(fit))
      lo <- c(lo, exp(ci[1]) / (1+exp(ci[1])))
      hi <- c(hi, exp(ci[2]) / (1+exp(ci[2])))
    }
  }
  return (list(probs=probs, lo=lo, hi=hi))
}



pdf(file="waterloo.pdf", width=15, height=5)
par(mar=c(5,8,1,1), mfrow=c(1,3), cex=1)

res <- est.freq(b529)
barplot(as.numeric(res$probs), horiz=T, main="B.1.1.529", adj=0, cex.main=1.5, 
        names.arg=paste(counts$site, format(counts$coldate, "%b %d"), counts$sample), 
        las=1, cex.names=0.6, xlim=c(0, 0.05),
        xlab="Estimated frequency", cex.lab=1.2,
        col=ifelse(res$lo > 0.01, 'salmon', 'grey'))
segments(x0=res$lo, x1=res$hi, y0=(1:length(res$probs)-0.45)*1.2, lwd=2)
abline(v=0.01, lty=2)

res <- est.freq(ba1)
barplot(as.numeric(res$probs), horiz=T, main="BA.1", adj=0, cex.main=1.5, 
        names.arg=paste(counts$site, format(counts$coldate, "%b %d"), counts$sample), 
        las=1, cex.names=0.6, xlim=c(0, 0.05),
        xlab="Estimated frequency", cex.lab=1.2,
        col=ifelse(res$lo > 0.01, 'salmon', 'grey'))
segments(x0=res$lo, x1=res$hi, y0=(1:length(res$probs)-0.45)*1.2, lwd=2)
abline(v=0.01, lty=2)

res <- est.freq(ba2)
barplot(as.numeric(res$probs), horiz=T, main="BA.2", adj=0, cex.main=1.5, 
        names.arg=paste(counts$site, format(counts$coldate, "%b %d"), counts$sample), 
        las=1, cex.names=0.6, xlim=c(0, 0.05),
        xlab="Estimated frequency", cex.lab=1.2,
        col=ifelse(res$lo > 0.01, 'salmon', 'grey'))
segments(x0=res$lo, x1=res$hi, y0=(1:length(res$probs)-0.45)*1.2, lwd=2)
abline(v=0.01, lty=2)

dev.off()




