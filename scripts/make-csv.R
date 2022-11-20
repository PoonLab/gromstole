args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("\n\nUsage: Rscript make-csv.R [input JSON] [output CSV]\n",
       "  input JSON:  JSON produced by estimate-freqs.R\n",
       "  output CSV:  path to write CSV\n\n")
}
require(jsonlite)
input <- jsonlite::read_json(args[1], simplifyVector = TRUE)
outfile <- args[2]

# deal with columns of all NAs in counts
coverage <- as.data.frame(sapply(unique(input$results$nucleotide), function(d) {
  input$results$coverage[which(input$results$nucleotide == d)]
}))
row.names(coverage) <- unique(input$results$sample)

counts <- as.data.frame(sapply(unique(input$results$nucleotide), function(d) {
  input$results$count[which(input$results$nucleotide == d)]
}))
row.names(counts) <- unique(input$results$sample)

# deal with missing metadata
coldate <- input$metadata$coldate
if (is.null(coldate)) {
  coldate <- NA
}
site <- input$metadata$site
if (is.null(site)) {
  site <- NA
}

fq <- counts / coverage
fq[coverage<10] <- NA

if (is.null(input$estimate$lower.95)) {
  lower.95 <- NA
} else {
  lower.95 <- 100*input$estimate$lower.95
}
if (is.null(input$estimate$upper.95)) {
  upper.95 <- NA
} else {
  upper.95 <- 100*input$estimate$upper.95
}
if (is.null(input$estimate$est)) {
  input$estimate$est <- NA
}

df <- data.frame(
  coldate=coldate,
  site=site,
  n.muts=apply(counts, 1, function(x) sum(x>0, na.rm=T)),
  cover=apply(coverage, 1, function(x) sum(x>0)),
  #min.freq=apply(fq, 1, function(x) 100*min(x, na.rm=T)),
  #med.freq=apply(fq, 1, function(x) 100*median(x, na.rm=T)),
  #max.freq=apply(fq, 1, function(x) 100*max(x, na.rm=T))
  est=100*input$estimate$est,
  lower.95=lower.95,
  upper.95=upper.95
)

write.table(t(c("", "Date collected", "Site", 
            "Positions with mutations detected", "Positions covered",
            "Estimated frequency (%)", "Lower 95% CI (%)", 
            "Upper 95% CI (%)")), outfile, sep=',', 
            quote=TRUE, col.names=FALSE, row.names=F)
write.table(df, outfile, sep=',', quote=T, col.names=FALSE, append=TRUE)
