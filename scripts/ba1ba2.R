## based on omicron-retro.R, priority analysis of BA1/BA2/B.1.1.529

#setwd('~/git/gromstole')

guelph <- "--guelph" %in% commandArgs(TRUE)

# load B.1.1.529/BA.1/BA.2 mutations list
omi <- read.csv('data/omicron-BA-fixed.csv')

# load background mutation frequencies (get-BA1BA2-uniques.py, from Nov 26)
bkgd <- read.csv('data/get-BA1BA2-uniques.csv')

# select mutations that are present in at most 5% of any other lineage
x <- apply(bkgd[,3:ncol(bkgd)], 2, max)
# mismatch: row 36 of omi is a duplicate; not present in columns of bkgd
# 
#x <- c(x[1:35], 0, x[36:length(x)])
omi <- omi[-36, ]
fin <- omi[which(x<0.05),]
#write.csv(fin, file="persei8.csv")
fin$label <- paste(fin$type, fin$pos, fin$alt, sep='|')


# locate data files
require(here)
if (guelph) {
  mfiles <- list.files(here("results/guelph/20211211_100707/"), full.names = TRUE,
                       pattern="*.mapped.csv$", recursive=TRUE)
  cfiles <- list.files(here("results/guelph/20211211_100707/"), full.names = TRUE,
                       pattern="*.coverage.csv$", recursive=TRUE)  
} else {
  mfiles <- list.files(here("results/waterloo/run6"), full.names = TRUE,
                       pattern="*.mapped.csv$", recursive=TRUE)
  cfiles <- list.files(here("results/waterloo/run6"), full.names = TRUE,
                       pattern="*.coverage.csv$", recursive=TRUE)  
}


# import coverage data
cover <- sapply(cfiles, function(f) {
  coverage <- read.csv(f, stringsAsFactors = FALSE)
  # NOTE: coverage$position is 0-index, omicron$pos is 1-index
  #coverage$coverage[which(is.element(coverage$position+1, fin$pos))]
  idx <- match(fin$pos, coverage$position+1)
  coverage$coverage[idx[!is.na(idx)]]
  #coverage$coverage[fin$pos-1]
})
row.names(cover) <- fin$label
cover <- as.data.frame(cover)
names(cover) <- sapply(strsplit(colnames(cover), split = "/"), function(x) {
  strsplit(x[length(x)], split = "\\.")[[1]][1]
})


# import mutation frequency data
maps <- sapply(mfiles, function(f) {
  mapped <- read.csv(f, stringsAsFactors = FALSE)
  mapped$type <- substr(mapped$label, 1, 1)
  mapped$pos <- mapped$position+1
  mapped$alt <- sapply(strsplit(mapped$label, as.character(mapped$pos)), 
                       function(x) {
                         gsub("^\\.", "", x[2])
                       })
  # convert to omicron label type
  mapped$label <- paste(mapped$type, mapped$pos, mapped$alt, sep='|')
  
  index <- match(fin$label, mapped$label)
  mapped$frequency[index]
})
row.names(maps) <- fin$label

maps <- as.data.frame(maps)

# set column names to sample identifiers
names(maps) <- sapply(strsplit(colnames(maps), split = "/"), function(x) {
  strsplit(x[length(x)], split = "\\.")[[1]][1]
})


# set zeroes for sites with non-zero coverage
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

# convert to integer counts
counts <- as.data.frame(t(maps*cover))


row.names(counts) <- NULL
names(counts) <- fin$mut_aa
counts$sample <- names(maps)
# switching to `here` library broke this line
counts$lab <- sapply(cfiles, function(x) {
  tokens <- strsplit(x, "/")[[1]]
  tokens[length(tokens)-2]
})


# load metadata
if (guelph) {
  m <- read.csv('uploads/guelph/20211211_100707/metadata_guelph_12Dec2021.csv')
  m <- m[1:30]
  m$Specimen.collector.sample.ID <- sapply(as.character(m$r1.fastq.filename), function(x) {
    strsplit(x, "_")[[1]][1]
  })  
} else {
  m <- read.csv("uploads/waterloo/run6/metadata.csv")  
}


if (guelph) {
  idx <- match(counts$sample, m$Specimen.collector.sample.ID) 
} else {
  idx <- match(counts$sample, m$specimen.collector.sample.ID)  
}


if (guelph) {
  counts$coldate <- as.Date(m$sample.collection.date[idx], 
                            format="%d-%b-%y"  #"%m/%d/%Y"
  )  
} else {
  counts$coldate[!is.na(idx)] <- as.Date(
    m$sample.collection.date[idx[!is.na(idx)]], 
    format="%b %d %Y")  #"%m/%d/%Y"
  counts$coldate <- as.Date(counts$coldate, origin="1970-01-01")  
}


counts$site <- m$geolocation.name..region.[idx]

cvr <- as.data.frame(t(cover))
names(cvr) <- names(counts)[1:ncol(cvr)]

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


b529 <- which(fin$lineage=='B.1.1.529')
ba1 <- which(fin$lineage=='BA.1' | fin$lineage=='B.1.1.529')
ba2 <- which(fin$lineage=='BA.2' | fin$lineage=='B.1.1.529')

# set output for PDF
if (guelph) {
  pdf(file="guelph2.pdf", width=15, height=5)  
} else {
  # waterloo
  pdf(file="waterloo.pdf", width=15, height=5)
}

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
