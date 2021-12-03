import json

print('loading JSON')
recoded = json.load(open("recoded.json"))
"""
{'union': {'~|503|A': 0, '~|2527|A': 1,
'labels': {'0': ['hCoV-19/Philippines/PH-PGC-39888/2021|EPI_ISL_4715831|2021-04-16'], 
'indexed': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], 
"""

mutations = {}
for lineage, ldata in recoded.items():
    print(lineage)
    mutations.update({lineage: {'mutations': {}, 'count': 0}})

    # invert union dict
    features = {}
    for f, fidx in ldata['union'].items():
        typ, pos, alt = f.split('|')
        key = '|'.join(map(str, [typ, int(pos)+1, alt]))
        features.update({fidx: key})
    
    for vidx, variant in enumerate(ldata['indexed']):
        # number of genomes in variant
        count = len(ldata['labels'][str(vidx)])
        mutations[lineage]['count'] += count  # denominator
        
        for fidx in variant:
            feature = features[fidx]  # map index to feature
            if feature not in mutations[lineage]['mutations']:
                mutations[lineage]['mutations'].update({feature: 0})
            mutations[lineage]['mutations'][feature] += 1
        
with open('count-mutations.json', 'w') as outfile:
    json.dump(mutations, outfile, indent=2)
