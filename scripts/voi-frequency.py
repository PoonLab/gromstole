import csv
import json
import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description='Get mutation frequencies from a list of variants.')
parser.add_argument('variant_list', metavar='variant_list', type=str, nargs='+',
                    help='list of variants of interest')
args = parser.parse_args()


# read in JSON of mutations in all lineages
with gzip.open("data/count-mutations_nsgb.json.gz", 'rt') as handle:
    mutations = json.load(handle)

#variant_list = ["B.1.1.529", "BA.1", "BA.2"]
variant_list = args.variant_list
for key in variant_list:
    if key not in mutations.keys():
        sys.exit("Error: " + key + " is not in nsgb mutations.")

variant = []
for key in variant_list:
    variant.extend(mutations[key]['mutations'])

outfile = open("data/voi-frequency_" + "_".join(variant_list).replace(".", "-") + ".csv", 'w')
outfile.write("lineage,variant_of_interest,count," + ",".join(variant) + '\n')
for lineage, ldata in mutations.items():
    denom = ldata['count']
    voi = 1*(lineage in variant_list)
    outfile.write("{},{},{}".format(lineage, voi, denom))

    for key in variant:
        count = ldata['mutations'].get(key, 0)
        outfile.write(",{}".format(count/denom))

    outfile.write('\n')

outfile.close()
