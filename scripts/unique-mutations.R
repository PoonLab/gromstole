require(jsonlite)

setwd("~/git/Untitled/data")
js <- read_json("retrieve.2022-11-27.json")

lineage.counts <- sapply(js$lineages, function(x) x$total)
sum(lineage.counts)  # total number of genomes in this database


hist(log10(lineage.counts), breaks=30)
rug(log10(lineage.counts))
table(lineage.counts<100)

# ignore lineages with fewer than 10 samples
js$lineages <- js$lineages[which(lineage.counts >= 10)]


# first, get mutation frequencies of focus lineage
focus <- "BA.5.2.1"

sublist <- js$lineages[[focus]]$mutations
total <- js$lineages[[focus]]$total

all.pos <- c()
all.alt <- c()
all.freq <- c()
for (pos in names(sublist)) {
  for (i in 1:length(sublist[[pos]])) {
    alt <- names(sublist[[pos]])[i]
    count <- sublist[[pos]][[i]]
    freq <- count / total
    if (freq > 0.5) {
      all.pos <- c(all.pos, pos)
      all.alt <- c(all.alt, alt)
      all.freq <- c(all.freq, freq)
    }
  }
}
foc.df <- data.frame(pos=all.pos, alt=all.alt, freq=all.freq)


# make a data frame for these mutations in all other lineages
others <- matrix(NA, nrow=length(js$lineages), ncol=nrow(foc.df))
row.names(others) <- names(js$lineages)
for (j in 1:length(js$lineages)) {
  if (names(js$lineages)[j] == focus) {
    next  # leave as NAs
  }
  ldata <- js$lineages[[j]]
  tot <- ldata$total
  for (i in 1:nrow(foc.df)) {
    pos <- foc.df$pos[i]
    alt <- foc.df$alt[i]
    if (is.null(ldata$mutations[[pos]]) || is.null(ldata$mutations[[pos]][[alt]])) {
      others[j, i] <- 0
    }
    else {
      others[j, i] <- ldata$mutations[[pos]][[alt]] / tot
    }
  }  
}

mean.freq <- apply(others, 2, function(x) sum(x, na.rm=T) / nrow(others))

require(ggfree)

pdf(file="BA-5-2-1.other-freq.pdf", height=10, width=5)
par(mar=c(5,5,1,1))
bpdata <- barplot(mean.freq, plot=F)
plot(NA, xlim=c(0, 1), ylim=range(bpdata), xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
add.grid(mode='y')
barplot(mean.freq, horiz=T, names.arg=paste0(foc.df$pos, foc.df$alt), 
        las=2, cex.names=0.5, add=T, border=NA, col='black')
dev.off()
