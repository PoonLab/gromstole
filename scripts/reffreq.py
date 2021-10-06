import gzip
import json

with gzip.open('../data/recoded.json.gz') as handle:
    recoded = json.load(handle)

for lineage, ldata in recoded.items():
    # convert feature dict to list
    features = [k for k, v in ldata['union'].items()]

    # first pass, generate individual mutation frequencies
    singles = {}
    for vidx, variant in enumerate(ldata['indexed']):
        labels = ldata['labels'][vidx]  # feature vector index to list of sequence names
        size = len(labels)

