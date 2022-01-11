suppressPackageStartupMessages(library(dplyr))
library(tidyr)
library(here)

# load mutations specific to Omicron variant
omicron <- read.csv(here("data", 'omicron-specific.csv'), 
  row.names = 1, stringsAsFactors = FALSE)
omicron <- omicron[!is.na(omicron$pos),]

omicron$X <- NULL  # remove duplicate row names
omicron$label <- paste(omicron$type, omicron$pos, omicron$alt, sep="|")


mfiles <- list.files("results", 
    full.names = TRUE, pattern="*.mapped.csv$", recursive=TRUE)
cfiles <- list.files("results", 
    full.names = TRUE, pattern="*.coverage.csv$", recursive=TRUE)

# determine coverage at omicron sites
cover <- sapply(cfiles, function(f) {
  coverage <- read.csv(f, stringsAsFactors = FALSE)
  # NOTE: coverage$position is 0-index, omicron$pos is 1-index
  coverage$coverage[which(is.element(coverage$position+1, omicron$pos))]
})
row.names(cover) <- omicron$label



maps <- sapply(mfiles, function(f) {
  mapped <- read.csv(f, stringsAsFactors = FALSE)
  mapped$type <- substr(mapped$label, 1, 1)
  mapped$pos <- mapped$position + 1
  mapped$alt <- sapply(
    strsplit(mapped$label, as.character(mapped$pos)), 
    function(x) {
        gsub("^\\.", "", x[2])
    })
  # convert to omicron label type
  mapped$label <- paste(mapped$type, mapped$pos, mapped$alt, 
    sep='|')
  
  index <- match(omicron$label, mapped$label)
  mapped$frequency[index]
})
row.names(maps) <- omicron$label

maps <- as.data.frame(maps)

# set column names to sample identifiers
names(maps) <- sapply(strsplit(colnames(maps), split = "/"), 
    function(x) {
        strsplit(x[length(x)], split = "\\.")[[1]][1]
})

# set zeroes for sites with non-zero coverage
cover <- as.data.frame(cover)
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

names(counts) <- ifelse(omicron$mut_aa=='None', 
    omicron$label, 
    omicron$mut_aa)

counts$sample <- names(maps)
# switching to `here` library broke this line
counts$lab <- sapply(cfiles, function(x) {
  tokens <- strsplit(x, "/")[[1]]
  tokens[length(tokens)-2]
})
counts$folder<- sapply(strsplit(cfiles, "/"), function(x) x[length(x) - 1])

# load metadata
metas <- list.files("uploads", 
    full.names = TRUE, pattern="*meta*", recursive=TRUE)
metas <- metas[!grepl("template", metas)]
metas <- metas[!grepl("fixed", metas)]
m <- lapply(metas, function(x) {
  # Reading everything in as a character to avoid type issues (logical v. character)
  meta <<- read.csv(x, colClasses = "character")
  names(meta) <- tolower(names(meta))

  fldrtmp <- strsplit(x, "/")[[1]]
  meta$folder <- fldrtmp[length(fldrtmp) - 1]
  meta$lab <- fldrtmp[length(fldrtmp) - 2]
  meta$sample <- sapply(strsplit(meta$r1.fastq.filename, "\\_"), function(x) x[1])
  if(any(meta$sample == "NE1-1001", na.rm = TRUE)) {
    meta$sample[meta$sample == "NE1-1001"] <- meta$specimen.collector.sample.id[meta$sample == "NE1-1001"]
  }

  meta
})
#m <- do.call(rbind, m)
m2 <- dplyr::bind_rows(m)
m3 <- m2 %>%
    select(sample, folder,
        coldate = sample.collection.date, 
        location = geolocation.name..region.,
        longitude = geolocation.latitude,
        latitude = geolocation.longitude) %>%
    filter(!is.na(sample)) %>%
    distinct()
counts <- left_join(counts, m3, by = c("sample", "folder"))


# make a nicer coverage file too
cvr <- as.data.frame(t(cover))

names(cvr) <- ifelse(omicron$mut_aa=='None', 
    omicron$label, 
    omicron$mut_aa)

cvr$sample <- counts$sample
cvr$lab <- counts$lab

# Long (tidy) data format
counts_long <- pivot_longer(counts, 
    -c(sample, lab, folder, coldate, location, latitude, longitude),
    names_to = "mutation", 
    values_to = "count")
cvr_long <- pivot_longer(cvr, -c(sample, lab),
    names_to = "mutation", 
    values_to = "coverage")

coco <- full_join(counts_long, cvr_long, by = c("sample", "lab", "mutation"))

write.csv(coco, file = here("data", "coco.csv"), row.names = FALSE)


