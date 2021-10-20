setwd('~/git/wastewater')

## load VCF from Virontus on Landgraff Nanopore data
#longshot <- read.table("data/results2/vcf/SRR13912812-MN908947.3.longshot.vcf")
#medaka <- read.table("data/results2/vcf/SRR13912812-MN908947.3.medaka.vcf")
#names(longshot) <- c('chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 
#                     'info', 'format', 'sample')
#names(medaka) <- c('chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 
#                     'info', 'format', 'sample')
#longshot$depth <- sapply(longshot$info, function(x) {
#  tokens <- strsplit(x, ";")[[1]]
#  dp <- as.integer(strsplit(tokens[1], "=")[[1]][2])
#})
#longshot$count <- sapply(longshot$info, function(x) {
#  tokens <- strsplit(x, ";")[[1]]
#  sum(as.integer(strsplit(strsplit(tokens[2], "=")[[1]][2], ",")[[1]]))
#})

## note these are substitutions only
#longshot$rel.freq <- longshot$count / longshot$depth
#longshot$label <- paste("~", longshot$pos, substr(longshot$alt, 1, 1), sep="")
#index <- match(longshot$pos, medaka$pos)
#longshot$medaka.q <- medaka$qual[index]


ivar <- read.csv('data/SRR13912812.ivar.tsv', sep='\t', header=T)

# load our results for Illumina data from same study
ont <- read.csv('data/ontario.csv')

ivar$label <- paste('~', ivar$POS, ivar$ALT, sep='')

index <- match(ont$label, ivar$label)
ont$ivar <- ivar$ALT_FREQ[index]

#longshot$illumina <- ont$frequency[index]


# Our try at Nanopore
nano <- read.csv('~/git/wastewater/data/SRR13912812.mapped.csv')

# coverage plot
temp <- cbind(nano$position, nano$coverage)
temp <- temp[!duplicated(temp), ]
temp <- temp[order(temp[,1]), ]
png("coverage2.png", width=6*150, height=5*150, res=150)
plot(temp[,1], temp[,2], type='s', log='y', xlab='Reference position', ylab='Read depth')
dev.off()

index <- match(ont$label, nano$label)
ont$nanopore <- nano$frequency[index]

res <- 72
png(filename="~/Desktop/nanopore.png", width=5*res, height=5*res, res=res)
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(ont$nanopore, ont$ivar, xlim=c(0,1), ylim=c(0,1),
     xlab='UWO Nanopore pipeline', ylab='iVar (Landgraff)')
abline(a=0, b=1, col=rgb(0,0,0,0.3))
dev.off()


png(filename="~/Desktop/platforms.png", width=5*res, height=5*res, res=res)
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(foo$frequency, foo$nanopore, xlab='UWO Illumina', ylab='UWO Nanopore', col=)
abline(a=0, b=1, col=rgb(0,0,0,0.3))
dev.off()

index <- match(longshot$label, nano$label)
longshot$our.nano <- nano$frequency[index]

plot(longshot$rel.freq, longshot$our.nano, type='n',
     xlab='Virontus frequency', ylab='Our estimates on same ONT data')
text(longshot$rel.freq, longshot$our.nano, label=longshot$label, cex=0.6)
