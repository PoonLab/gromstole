#awk -F "," '{print $1}' lineages.csv | awk 'NR!=1' | less -S > lineage_names.txt
#xzgrep -f lineage_names.txt metadata.tsv > test.txt # Killed (maxed out the memory)

import gzip
import csv

max_iter = 31000

listnames = open("data/lineages.csv")

lins = csv.reader(listnames)
_ = next(lins)

# Populate with all lineages
pangomiss = []
counter = 0
for taxon, lineage in lins:
    pangomiss.append(taxon)
    counter += 1
    if counter > max_iter:
        break
print(len(pangomiss))

# If the lineage exists in metadata, remove it from pangomiss
counter = 0
with gzip.open("data/metadata.tsv.gz", "rt") as metadata:
    for line in metadata:
        taxon = line.split("\t")[0]
        if taxon in pangomiss:
            pangomiss.remove(taxon)
        counter += 1
        if counter > max_iter:
            break

print(len(pangomiss))