import lzma
import json

handle = lzma.open("../data/wastewater.json.xz")
outfile = open("../data/wastewater_names.txt", "w")

for row in handle:
    json_dict = json.loads(row)
    seq_name = json_dict["covv_virus_name"].split("/")
    seq_name.remove("hCoV-19")

    seq_lineage = json_dict["covv_lineage"]

    if len(seq_lineage) == 0:
        print("/".join(seq_name))

    if len(seq_lineage) > 0:
        outfile.writelines("/".join(seq_name)+"\n")
    
"""
#!/bin/bash

cd data
# Delete and re-create so that it's empty
touch sequences_pangolin2.fasta && rm sequences_pangolin2.fasta && touch sequences_pangolin2.fasta

while read line; do
    echo "$line"
    xzgrep -A 1 -m 1 "$line" sequences.fasta.xz >> sequences_pangolin2.fasta 
done < wastewater_names.txt

# How big is the difference?
grep -c ">" sequences_pangolin2.fasta
wc -l wastewater_names.txt

xz sequences_pangolin2.fasta

"""