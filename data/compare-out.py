import csv
import json
import gzip

# read in mutations observed in Omicron genomes
omicron = {}
for row in csv.DictReader(open("omicron-mutations.csv")):
    key = "{type}|{pos}|{alt}".format(**row)
    omicron.update({key: float(row['freq'])})

# read in JSON of mutations in all lineages
with gzip.open("count-mutations.json.gz", 'rt') as handle:
    mutations = json.load(handle) 

outfile = open("compare-out.csv", 'w')
outfile.write("lineage,count," + ",".join(omicron.keys()) + '\n')

for lineage, ldata in mutations.items():
    denom = ldata['count']
    outfile.write("{},{}".format(lineage, denom))
    
    for key, freq in omicron.items():
        count = ldata['mutations'].get(key, 0)
        outfile.write(",{}".format(count/denom)) 
 
    outfile.write('\n')

outfile.close()
