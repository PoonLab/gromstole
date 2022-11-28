require(lubridate)
require(ggfree)

setwd("~/git/gromstole/results/")

# gather CSV files
gather.csv <- function(variant) {
  files <- Sys.glob(paste("guelph/*/guelph*", variant, ".*.csv", sep=""))
  guelph <- do.call("rbind", lapply(files, read.csv))
  guelph$lab <- 'guelph'
  
  files <- Sys.glob(paste("waterloo/*/waterloo*", variant, ".*.csv", sep=""))
  waterloo <- do.call("rbind", lapply(files, read.csv))
  waterloo$lab <- 'waterloo'
  
  files <- Sys.glob(paste("western/*/western*", variant, ".*.csv", sep=""))
  western <- do.call("rbind", lapply(files, read.csv))
  western$lab <- 'western'
  western <- western[!grepl("BARR", western$X), ]
  
  df <- rbind(guelph, waterloo, western)
  #df <- df[!is.na(df$Estimate), ]
  df$coldate <- as.Date(df$Date.collected)
  # ignore samples with dates earlier than 2021-09-01, probably mis-parsed
  df[!is.na(df$coldate) & df$coldate > as.Date('2021-09-01') & df$coldate < today(),]
}

delta <- gather.csv("B16172")
ba1 <- gather.csv("BA1")
ba2 <- gather.csv("BA2")
ba4 <- gather.csv("BA4")
ba5 <- gather.csv("BA5")
ba275 <- gather.csv("BA275")
be1 <- gather.csv("BE1")

pal <- rep(add.alpha(c("#C20430", "#ffd54f", "#4F2683"), 0.5), 
           table(delta$lab))

months <- seq(as.Date('2021-09-01'), today(), by='months')

# ===============
png(paste("longitudinal", today(), "png", sep='.'), width=9*600, height=9*600, res=600)
par(mar=c(5,5,1,1), mfrow=c(3,3))

plot(delta$coldate, delta$Estimate, bg=pal, pch=21, main='Delta (B.1.617.2)',
     xlab='Sampling date', ylab='Frequency (%)', xaxt='n', lwd=0.75,
     cex=ifelse(delta$Lower.95..CI....<1 | is.na(delta$Lower.95..CI....), 0.5, 1))
axis.srt(side=1, at=months, label=format(months, "%b '%y"), srt=30, cex=0.8)
legend(x=as.Date('2022-04-01'), y=90, legend=c('Guelph', 'Waterloo', 'Western'), 
       fill=unique(pal), bty='n', cex=0.8)

plot(ba1$coldate, ba1$Estimate, bg=pal, pch=21, main='Omicron (BA.1)',
     xlab='Sampling date', ylab='Frequency (%)', xaxt='n', lwd=0.75,
     cex=ifelse(ba1$Lower.95..CI....<1 | is.na(ba1$Lower.95..CI....), 0.5, 1))
axis.srt(side=1, at=months, label=format(months, "%b '%y"), srt=30, cex=0.8)
legend(x=as.Date('2021-09-10'), y=90, legend=c('Guelph', 'Waterloo', 'Western'), 
       fill=unique(pal), bty='n', cex=0.8)

plot(ba2$coldate, ba2$Estimate, bg=pal, pch=21, main='Omicron (BA.2)',
     xlab='Sampling date', ylab='Frequency (%)', xaxt='n', lwd=0.75,
     cex=ifelse(ba2$Lower.95..CI....<1 | is.na(ba2$Lower.95..CI....), 0.5, 1))
axis.srt(side=1, at=months, label=format(months, "%b '%y"), srt=30, cex=0.8)
legend(x=as.Date('2021-09-10'), y=90, legend=c('Guelph', 'Waterloo', 'Western'), 
       fill=unique(pal), bty='n', cex=0.8)

plot(ba4$coldate, ba4$Estimate, bg=pal, pch=21, main='Omicron (BA.4)',
     xlab='Sampling date', ylab='Frequency (%)', xaxt='n', lwd=0.75, 
     cex=ifelse(ba4$Lower.95..CI....<1 | is.na(ba4$Lower.95..CI....), 0.5, 1),
     ylim=c(0, 100))
axis.srt(side=1, at=months, label=format(months, "%b '%y"), srt=30, cex=0.8)
legend(x=as.Date('2021-09-10'), y=90, legend=c('Guelph', 'Waterloo', 'Western'), 
       fill=unique(pal), bty='n', cex=0.8)

plot(ba5$coldate, ba5$Estimate, bg=pal, pch=21, main='Omicron (BA.5)',
     xlab='Sampling date', ylab='Frequency (%)', xaxt='n', lwd=0.75, 
     cex=ifelse(ba5$Lower.95..CI....<1 | is.na(ba5$Lower.95..CI....), 0.5, 1),
     ylim=c(0, 100))
axis.srt(side=1, at=months, label=format(months, "%b '%y"), srt=30, cex=0.8)
legend(x=as.Date('2021-09-10'), y=90, legend=c('Guelph', 'Waterloo', 'Western'), 
       fill=unique(pal), bty='n', cex=0.8)

okay <- ifelse(ba275$Lower.95..CI....<1 | is.na(ba275$Lower.95..CI....) | 
                 ba275$Positions.covered < 5, 0.5, 1)

plot(ba275$coldate, ba275$Estimate, bg=pal, pch=21, main='Omicron (BA.2.75)',
     xlab='Sampling date', ylab='Frequency (%)', xaxt='n', lwd=0.75, 
     cex=okay, ylim=c(0, 100))
axis.srt(side=1, at=months, label=format(months, "%b '%y"), srt=30, cex=0.8)
legend(x=as.Date('2021-09-10'), y=90, legend=c('Guelph', 'Waterloo', 'Western'), 
       fill=unique(pal), bty='n', cex=0.8)


plot(be1$coldate, be1$Estimate, bg=pal, pch=21, main='Omicron (BE.1)',
     xlab='Sampling date', ylab='Frequency (%)', xaxt='n', lwd=0.75, 
     cex=ifelse(be1$Lower.95..CI....<1 | is.na(be1$Lower.95..CI....), 0.5, 1),
     ylim=c(0, 100))
axis.srt(side=1, at=months, label=format(months, "%b '%y"), srt=30, cex=0.8)
legend(x=as.Date('2021-09-10'), y=90, legend=c('Guelph', 'Waterloo', 'Western'), 
       fill=unique(pal), bty='n', cex=0.8)

dev.off()

