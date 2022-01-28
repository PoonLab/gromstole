args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("\n\nUsage: Rscript make-csv.R [input JSON] [output CSV]\n",
       "  input JSON:  JSON produced by estimate-freqs.R\n",
       "  output PDF:  path to write CSV\n\n")
}
require(jsonlite)
input <- jsonlite::read_json(args[1], simplifyVector = TRUE)
outfile <- args[2]

# date collected, site, sequencing platform, number of mutations, frequency range, interpretation
counts <- input$counts[names(input$coverage)]
fq <- counts/input$coverage
fq[input$coverage<10] <- NA

df <- data.frame(
  coldate=input$metadata$coldate,
  site=input$metadata$site,
  platform="Illumina",
  n.muts=apply(counts, 1, function(x) sum(x>0, na.rm=T)),
  cover=apply(input$coverage, 1, function(x) sum(x>0)),
  min.freq=apply(fq, 1, function(x) 100*min(x, na.rm=T)),
  med.freq=apply(fq, 1, function(x) 100*median(x, na.rm=T)),
  max.freq=apply(fq, 1, function(x) 100*max(x, na.rm=T))
)

write.table(t(c("", "Date collected", "Site", "Sequencing platform",
            "Positions with mutations detected", "Positions covered",
            "Minimum frequency (%)", "Median frequency (%)", 
            "Maximum frequency (%)")), outfile, sep=',', 
            quote=TRUE, col.names=FALSE, row.names=F)
write.table(df, outfile, sep=',', quote=T, col.names=FALSE, append=TRUE)
