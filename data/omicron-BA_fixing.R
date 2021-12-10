# Fixing Omicron data

omicron <- read.csv("https://github.com/cov-lineages/pango-designation/files/7668225/Omicron_BA.1_BA.2_mutations.csv",
    stringsAsFactors = FALSE)
names(omicron)[c(1,4)] <- c("mut_nuc", "mut_aa")
omicron <- omicron[, -3] # column name is "X", all entries are NA
write.csv(omicron, file = "data/omicron-BA-raw.csv", row.names = FALSE)
head(omicron)

# Finding the Type
omicron$type <- NA
# Short names means mutation
omicron$type[nchar(omicron$mut_nuc) <= 7] <- "~"
omicron$type[grepl(", ", omicron$mut_nuc)] <- "~"
# Deletions and insertions are (mostly) marked
omicron$type[grepl("del", omicron$mut_nuc)] <- "-"
omicron$type[grepl("_", omicron$mut_nuc) & !grepl("(main)", omicron$mut_nuc)] <- "-"
omicron$type[grepl("ins", omicron$mut_nuc)] <- "+"

omicron[, c("type", "mut_nuc", "mut_aa")]

# Finding the position
omicron$pos <- NA
easy_muts <- which(nchar(omicron$mut_nuc) <= 7)
omicron$pos[easy_muts] <- unlist(substr(omicron$mut_nuc[easy_muts], 2, 
    nchar(omicron$mut_nuc[easy_muts]) - 1))
dels <- which(omicron$type == "-")
omicron$pos[dels] <- sapply(strsplit(omicron$mut_nuc[dels], split = "_"), `[`, 1)
inss <- which(omicron$type == "+")
omicron$pos[inss] <- sapply(strsplit(omicron$mut_nuc[inss], "[GATCins]"), `[`, 1)

omicron[, c("pos", "mut_nuc", "type")]



# Finding Alt (the hard one)
omicron$alt <- NA
omicron$alt[easy_muts] <- substr(omicron$mut_nuc[easy_muts], 
    nchar(omicron$mut_nuc[easy_muts]), 
    nchar(omicron$mut_nuc[easy_muts]))
multi_muts <- which(grepl(", ", omicron$mut_nuc))
omicron$alt[multi_muts] <- sapply(strsplit(omicron$mut_nuc[multi_muts], split = ", "),
    function(x) paste0(sapply(x,
        function(y) {
            if(nchar(y) == 0) {
                NULL
            } else {
                substr(y, nchar(y), nchar(y))
            }
        }), collapse = ""))
omicron$alt[dels] <- sapply(strsplit(omicron$mut_nuc[dels], split = "[_del]"), 
    function(x) {
        as.numeric(x[2]) - as.numeric(x[1])
    })
omicron$alt[inss] <- sapply(strsplit(omicron$mut_nuc[inss], split = "[(0-9)ins]"), 
    function(x) {
        x[nchar(x) > 0]
    })

omicron[, c("alt", "mut_nuc")]


# Dealing with (main)
omicron[grepl("(main)", omicron$mut_nuc), c("type", "pos", "alt")] <- c("-", 21987, 8)

# One column for lineage
omicron$lineage <- apply(omicron[, c("BA.1", "BA.2", "B.1.1.529")], 1, 
    function(x) {
        c("BA.1", "BA.2", "B.1.1.529")[which(x == "Y")]
    }
)

write.csv(omicron, file = "data/omicron-BA.csv", row.names = FALSE)

