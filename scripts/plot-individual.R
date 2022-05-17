library("stringr")

args <- commandArgs(trailingOnly=TRUE)
numArgs <- length(args)

mapfile <- file.path(args[1], str_c(args[2], '.mapped.csv'))
covfile <- file.path(args[1], str_c(args[2], '.coverage.csv'))

coverage <- read.csv(covfile)

# draw SARS-CoV-2 genome
draw.sc2 <- function() {
  rect(xleft=1, xright=265, ybottom=0.5e6, ytop=1.2e6)  # 5'UTR
  
  rect(xleft=266, xright=13468, ybottom=0.5e6, ytop=1.2e6)  # ORF1a
  text(x=mean(c(266, 13468)), y=0.9e6, cex=0.6, label='ORF1a')
  
  rect(xleft=13468, xright=21555, ybottom=0.5e6, ytop=1.2e6)  # ORF1b
  text(x=mean(c(13468, 21555)), y=0.9e6, cex=0.6, label='ORF1b')
  
  rect(xleft=21563, xright=25384, ybottom=0.5e6, ytop=1.2e6)  # S
  text(x=mean(c(21563, 25384)), y=0.9e6, cex=0.6, label='S')
  
  rect(xleft=25393, xright=26220, ybottom=0.5e6, ytop=1.2e6)  # ORF3a
  rect(xleft=26245, xright=26472, ybottom=0.5e6, ytop=1.2e6)  # E
  rect(xleft=26523, xright=27191, ybottom=0.5e6, ytop=1.2e6) # M
  rect(xleft=27202, xright=27387, ybottom=0.5e6, ytop=1.2e6)  # orf6
  rect(xleft=27394, xright=27759, ybottom=0.5e6, ytop=1.2e6)  # orf7a
  rect(xleft=27756, xright=27887, ybottom=0.5e6, ytop=1.2e6)  # orf7b
  rect(xleft=27894, xright=28259, ybottom=0.5e6, ytop=1.2e6)  # orf8
  rect(xleft=28274, xright=29533, ybottom=0.5e6, ytop=1.2e6)  # N
  rect(xleft=29558, xright=29674, ybottom=0.5e6, ytop=1.2e6)  # orf10
}

covfilename <- file.path(args[1], str_c(args[2], '.coverage.png'))

# coverage
res <- 300
png(file=covfilename, width=11*res, height=8*res, res=res)
par(mar=c(7,7,7,7))
df <- coverage
df[,1] <- df[,1] + 1  # shift from zero-index
df[,2][df[,2]==0] <- 0.1

plot(df[,1], df[,2], type='s', log='y', las=1, ylim=c(1, 1e6),
     xlab='Reference position', ylab='Read depth', yaxt='n')
axis(side=2, at=c(1, 10, 100, 1000, 1e4, 1e5, 1e6), las=1,
     labels=parse(text=paste("10^", 0:6)))
polygon(c(0.1, df[,1], 29903, 0), c(0, df[,2], 0.1, 0.1), col='grey')
title(main=args[2], adj=0)
draw.sc2()

invisible(dev.off())

# map <- read.csv(mapfile)

# # B.1.617.2 mutations
# delta <- c('aa:orf1b:P314L', 'aa:S:T19R', 'del:22029:6', 'aa:S:L452R', 
#            'aa:S:T478K', 'aa:S:D614G', 'aa:S:P681R', 'aa:S:D950N', 
#            'aa:orf3a:S26L', 'aa:M:I82T', 'aa:orf7a:V82A', 'aa:orf7a:T120I', 
#            'del:28271:1', 'aa:N:D63G', 'aa:N:R203M', 'aa:N:D377Y', 
#            'aa:orf1a:A1306S', 'aa:orf1a:P2046L', 'aa:orf1a:P2287S', 
#            'aa:orf1a:V2930L', 'aa:orf1a:T3255I', 'aa:orf1a:T3646A', 
#            'aa:orf1b:G662S', 'aa:orf1b:P1000L', 'aa:orf1b:A1918V', 
#            'aa:orf7b:T40I', 'del:28248:6', 'aa:N:G215C')

# alpha <- c('aa:orf1a:T1001I', 'aa:orf1a:A1708D', 'aa:orf1a:I2230T', 
#            'del:11288:9', 'aa:orf1b:P314L', 'del:21765:6', 'del:21991:3', 
#            'aa:S:N501Y', 'aa:S:A570D', 'aa:S:D614G', 'aa:S:P681H', 
#            'aa:S:T716I', 'aa:S:S982A', 'aa:S:D1118H', 'aa:orf8:Q27*', 
#            'aa:orf8:R52I', 'aa:orf8:Y73C', 'aa:N:D3H', 'aa:N:D3V', 'aa:N:D3E', 
#            'aa:N:R203K', 'aa:N:G204R', 'aa:N:S235F', 'del:28271:1')

# pal <- rev(hcl.colors(n=8))
# # assume a maximum of 500,000
# pal.idx <- exp(seq(0, log(5e5), length.out=8))
# pal.idx[1] <- 0

# mapfilename <- file.path(args[1], str_c(args[2], '.delta.png'))

# res <- 300
# png(filename=mapfilename, width=12*res, height=8*res, res=res)
# par(mar=c(7,7,7,7))

# df <- map
# temp <- df[df$coverage > 100, ]
# temp <- temp[temp$label != '~9053C', ]  # hack
# idx <- unlist(sapply(delta, function(m) {
#   if(is.element(m, temp$mutation)) {
#     which(temp$mutation == m)
#   } else {
#     NA
#   }
# }))
# freq <- temp$frequency[idx]
# names(freq) <- names(idx)
# count <- temp$coverage[idx]
# count[is.na(count)] <- 0

# barplot(freq, horiz=T, las=1, xlim=c(0,1), cex.names=0.75,
#         col=pal[as.integer(cut(count, breaks=pal.idx))],
#         xlab='Estimated frequency')
# text(x=ifelse(is.na(freq), 0, freq)+0.01, 
#      y=(1:length(freq)-0.5)*1.2, 
#      label=count, xpd=NA, adj=0, cex=0.8)
# title(main=args[2], adj=0)

# invisible(dev.off())

