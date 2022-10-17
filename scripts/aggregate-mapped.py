import gzip
import csv
from datetime import datetime
from epiweeks import Week
import json


def parse_date(dt, formats=('%Y-%m-%d', '%d-%b-%y', '%m/%d/%Y')):
    """ Try multiple date formats """
    if dt in ['UNK']:
        return None
    for fmt in formats:
        try:
            return datetime.strptime(dt, fmt)
        except ValueError:
            pass
    raise ValueError(f"No supported format detected for date string {dt}")


# import location to region map as a dict
regions = json.load(open("data/regions.json"))

print("Aggregating data from input CSV...")
data = {}  # region / week / mutation
handle = gzip.open("data/collate-mapped.csv.gz", 'rt')
reader = csv.DictReader(handle)
for i, row in enumerate(reader):
    if i > 1000:
        break
    location = row['region']
    region = regions.get(location, None)
    if region is None:
        continue
    if region not in data:
        data.update({region: {}})

    # convert collection date to epiweek
    coldate = parse_date(row['coldate'])
    epiweek = Week.fromdate(coldate.date())
    if epiweek not in data[region]:
        data[region].update({epiweek: {}})

    mutation = (row['label'], row['mutation'])
    if mutation not in data[region][epiweek]:
        data[region][epiweek].update({
            mutation: {'count': 0, 'coverage': 0, 'samples': 0}
        })
    freq = float(row['freq'])
    coverage = int(row['coverage'])
    data[region][epiweek][mutation]['count'] += freq*coverage
    data[region][epiweek][mutation]['coverage'] += coverage
    data[region][epiweek][mutation]['samples'] += 1


# output aggregated data
print("Writing outputs...")
outfile = open("aggregate-mapped.csv", 'w')
outfile.write("region,year,epiweek,nuc,amino,nsamples,count,coverage\n")
for region, rdata in data.items():
    for epiweek, edata in rdata.items():
        for mutation, mdata in edata.items():
            outfile.write(
                f"{region},{epiweek.year},{epiweek.week},{mutation[0]},"
                f"{mutation[1]},{mdata['samples']},{mdata['count']},"
                f"{mdata['coverage']}\n"
            )
outfile.close()
