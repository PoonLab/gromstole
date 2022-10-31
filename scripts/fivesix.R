setwd("~/git/gromstole")
#ag <- read.csv(gzfile("data/aggregated-min5.csv.gz"), row.names=1)
ag <- read.csv("data/aggregate-mapped.csv")
ag <- ag[ag$count >= 5, ]

# use regular expressions to filter for mutations of interest
pat <- "aa:S:(R346[KT]|K444[MRT]|L452[MR]|N460K|F486[SV]|Q493R)"
df <- ag[grepl(pat, ag$amino), ]

# count is number of reads carrying mutation in any sample
df$freq <- df$count / df$coverage

# merge two northern regions, low number of samples
df$region[df$region=="North West"] <- "North"
df$region[df$region=="North East"] <- "North"
regions <- unique(df$region)

# restrict to current year, for convenience
df2 <- df[df$year==2022, ]
df2 <- df2[order(df2$epiweek), ]  # sort by week for clean line plots

# invert Q493R reversion (Q is reference amino acid)
# associated with compensatory effect on receptor affinity (BA.4/5)
df2$freq[df2$amino=="aa:S:Q493R"] <- 1 - df2$freq[df2$amino=="aa:S:Q493R"]
df2$amino[df2$amino=="aa:S:Q493R"] <- "aa:S:R493Q"

# get AA position (shared by multiple substitutions)
df2$pos <- gsub("aa:S:[A-Z]+([0-9]+)[A-Z]+", "\\1", df2$amino)
all.pos <- sort(unique(df2$pos))

# plotting code
require(ggfree)
pal <- gg.rainbow(3, l=60)

for (reg in regions) {
  df3 <- df2[df2$region==reg, ]
  
  # prepare plotting device
  pdf(paste(reg, "pdf", sep='.'), width=12, height=8)
  par(mfrow=c(2,3), mar=c(5,5,1,1), cex=0.9)
  
  for (i in 1:length(all.pos)) {
    pos <- all.pos[i]
    temp <- df3[df3$pos == pos, ]
    
    # skip empty data frame
    if (nrow(temp) == 0) {
      plot(NA, xlim=c(0, 1), ylim=c(0, 1), bty='n', xaxt='n', 
           yaxt='n', type='n', xlab='', ylab='')
      next
    }
    
    plot(df2$epiweek, df2$freq, type='n', ylim=c(0, 1.1), 
         xlab="Epiweek (2022)", cex.main=1.1, ylab="Mutation frequency")
    
    if (pos == all.pos[1]) title(main=reg, adj=0)
    
    mutations <- unique(df2$amino[df2$pos==pos])
    legend(x=26, y=1.15, xjust=0.5, xpd=NA, horiz=TRUE, legend=mutations, 
           cex=0.7, bty='n', col=pal[1:length(mutations)],
           pch=19, lty=1)
    
    for (j in 1:length(mutations)) {
      mut <- mutations[j]
      idx <- temp$amino == mut
      lines(temp$epiweek[idx], temp$freq[idx], col=pal[j])
      points(temp$epiweek[idx], temp$freq[idx], 
             col='white', bg=pal[j], pch=21, #cex=sqrt(temp$nsamples[idx])/2)
             cex=sqrt(temp$coverage[idx]/temp$nsamples[idx])/20)
    }
  }  
  dev.off()
}

