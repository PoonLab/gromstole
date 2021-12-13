## based on omicron-retro.R, priority analysis of BA1/BA2/B.1.1.529

setwd('~/git/gromstole')

# load B.1.1.529/BA.1/BA.2 mutations list
omi <- read.csv('data/omicron-BA-fixed.csv')

# load background mutation frequencies (get-BA1BA2-uniques.py)
bkgd <- read.csv('data/get-BA1BA2-uniques.csv')

# since we don't have frequency info for BA.*, just check background freqs
x <- apply(bkgd[3:ncol(bkgd)], 2, max)
par(mar=c(5,5,1,1))
barplot(x)


fin <- omi[which(x<0.05),]
#write.csv(fin, file="persei8.csv")

fin$label <- paste(fin$type, fin$pos, fin$alt, sep='|')

require(here)
mfiles <- list.files(here("results/waterloo/run6"), full.names = TRUE,
                     pattern="*.mapped.csv$", recursive=TRUE)
cfiles <- list.files(here("results/waterloo/run6"), full.names = TRUE,
                     pattern="*.coverage.csv$", recursive=TRUE)


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

counts <- as.data.frame(t(maps*cover))


row.names(counts) <- NULL
names(counts) <- fin$mut_aa
counts$sample <- names(maps)
# switching to `here` library broke this line
counts$lab <- sapply(cfiles, function(x) {
  tokens <- strsplit(x, "/")[[1]]
  tokens[length(tokens)-2]
})
counts$coldate <- NA
counts$site <- NA


