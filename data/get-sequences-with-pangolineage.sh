# Get data:
    # cd data
    # ./download.sh

# Unzip 
# (-k means keep original)
# Unfortunately, awk and seqtk subseq both need unzipped files
#unxz -kv sequences.fasta.xz 
gunzip -kv metadata.tsv.gz

# Get correct rows 
# (column 20 is pangolineage, 1 is fasta description, row 1 is column headers)
echo "Dumping fasta descriptors with known pangolineages to metadata_pango-names.txt"
awk -F '\t' '{ if ($20 != "?") && ($16 == "Homo sapiens") { print $1} }' metadata.tsv | awk 'NR!=1' > metadata_pango-names.txt
echo "taxon,lineage" > lineages2.csv
awk -F '\t' '{ if (($20 != "?") && ($16 == "Homo sapiens"))  { print $1","$20} }' metadata.tsv | awk 'NR!=1' >> lineages2.csv

# Find seqs
echo "Extracting sequences with labelled pangolineages."
#seqtk subseq sequences.fasta metadata_pango-names.txt > sequences_pangolin.fasta


# Clean up
echo "Zipping pangolineages and removing unzipped files"
#xz -v sequences_pangolin.fasta # don't keep original
#rm sequences.fasta
rm metadata.tsv
