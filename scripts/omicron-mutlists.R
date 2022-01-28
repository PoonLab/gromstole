
# load B.1.1.529/BA.1/BA.2 mutations list
omi <- read.csv('data/omicron-BA-final.csv')

# load background mutation frequencies (get-BA1BA2-uniques.py, from Nov 26)
bkgd <- read.csv('data/get-BA1BA2-uniques.csv')

# select mutations that are present in at most 5% of any other lineage
x <- apply(bkgd[,3:ncol(bkgd)], 2, max)
fin <- omi[which(x<0.05),]
fin$label <- paste(fin$type, fin$pos, fin$alt, sep='|')
#write.csv(fin, "~/Desktop/fin.csv")


# TODO: generalize to mutation list
b529 <- fin[which(fin$lineage=='B.1.1.529'), c("type", "pos", "alt", "mut_aa")]
names(b529) <- c("type", "pos", "alt", "label")
write.csv(b529, "lineages/B.1.1.529.csv", row.names=F)

ba1 <- fin[which(fin$lineage=='BA.1'),  #| fin$lineage=='B.1.1.529'), 
            c("type", "pos", "alt", "mut_aa")]
names(ba1) <- c("type", "pos", "alt", "label")
write.csv(ba1, "lineages/BA.1.csv", row.names=F)

ba2 <- fin[which(fin$lineage=='BA.2'),  #| fin$lineage=='B.1.1.529'), 
            c("type", "pos", "alt", "mut_aa")]
names(ba2) <- c("type", "pos", "alt", "label")
write.csv(ba2, "lineages/BA.2.csv", row.names=F)

