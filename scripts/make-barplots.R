args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("\n\nUsage: Rscript make-barplots.R [input JSON] [output PDF]\n",
       "  input JSON:  JSON produced by estimate-freqs.R\n",
       "  output PDF:  path to write PDF\n\n")
}
require(jsonlite)
input <- jsonlite::read_json(args[1], simplifyVector = TRUE)
outfile <- args[2]

# parse JSON
run.dir <- input$run.dir
lineage <- input$lineage
cvr <- input$coverage[names(input$counts)]
estimate <- input$estimate
metadata <- input$metadata

if (is.null(metadata$coldate)) {
  metadata$coldate <- NA 
} else {
  metadata$coldate <- as.Date(metadata$coldate)
}
if (is.null(metadata$site)) {
  metadata$site <- NA
}

if (is.null(estimate$lower.95)) {
  estimate$lower.95 <- NA
} 

if (is.null(input$estimate$upper.95)) {
  estimate$upper.95 <- NA
} 

if (all(is.na(metadata$site))) {
  names.arg <- metadata$sample
} else {
  names.arg <- paste(substr(metadata$site, 1, 10), format(metadata$coldate, "%b %d"), 
                     metadata$sample)  
}

pdf(file=outfile, width=7, height=max(5, nrow(input$counts)/5))

n.muts <- apply(cvr, 1, function(x) sum(x>0, na.rm=T))
pal <- colorRampPalette(c("white", "steelblue"))(ncol(cvr))

par(mar=c(5,8,2,8), mfrow=c(1,1), cex=1)
mp <- barplot(as.numeric(estimate$est), horiz=T, main=NA, cex.main=1.3, 
        names.arg=names.arg,
        las=1, cex.names=0.6, xlim=c(0, 1),  #max(res$hi, na.rm=T)),
        xlab="Estimated frequency", cex.lab=1.2,
        col=pal[n.muts], 
        border=ifelse(estimate$lower.95 > 0.01, "black", "grey60")
        )
        #col=ifelse(res$lo > 0.01, 'salmon', 'grey'))
title(main=paste(" ", lineage), font.main=1, cex.main=1.6, 
      adj=0, outer=TRUE, line=-1.5)
segments(x0=estimate$lower.95, x1=estimate$upper.95, y0=mp, 
         lwd=2, col=ifelse(estimate$upper.95 - estimate$lower.95 >= 0.9, "red", 
                           ifelse(estimate$lower.95 > 0.01, "black", "grey60")))
abline(v=0.01, lty=2)

u <- par('usr')  # user coordinates of plot region (x1, x2, y1, y2)

# add heatmap
pal2 <- colorRampPalette(c("white", "steelblue4"))(4)
dy <- diff(mp)[1]
par(xpd=NA)
for (i in 1:nrow(cvr)) {
  for (j in 1:ncol(cvr)) {
    rect(xleft=1.05+0.3*(j-1)/ncol(cvr), xright=1.05+0.3*j/ncol(cvr),
         ybottom=dy*(i-1), ytop=dy*(i), 
         col=pal2[cut(cvr[i,j], breaks = c(-Inf, 0, 10, 100, Inf))],
         border='white', lwd=0.5)
  }
}

text(x=dy, y=dy*(nrow(cvr)+2), adj=0.5, label="Coverage", cex=0.8)
legend(x=dy, y=dy*(nrow(cvr)+2), legend=c("1-9", "10-99", "100+"), 
       fill=pal2[2:4], bty='n', horiz=T, adj=-0.1, xjust=0.5, cex=0.7, x.intersp=0)
title(xlab=paste(" ", run.dir), outer=TRUE, adj=0, cex.lab=0.5, 
      font=1, line=-1)
par(xpd=FALSE)

abline(h = 0:ncol(cvr)*dy[1]+0.1, col = "gray88", lty = 2)

dev.off()

