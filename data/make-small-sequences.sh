# unit test-able data
# Goal: Make files small enough that I can test commands without applying to the whole dang file

#cd data

# Gather metadata
head -200 metadata_pango-names.txt > metadata_small.txt
zcat metadata.tsv.gz | head -1 > metadata_pango-names_small.tsv
zgrep -f metadata_small.txt metadata.tsv.gz >> metadata_pango-names_small.tsv

# Extract relevant sequences
xzgrep -A1 -f metadata_small.txt sequences.fasta.xz > sequences_small_tmp.fasta
grep -v -e "^-" sequences_small_tmp.fasta > sequences_small.fasta
#rm sequences_small_tmp.fasta
xz -k sequences_small.fasta


