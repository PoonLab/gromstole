from glob import glob
import subprocess
import os
import sys
import csv
import gzip

fieldnames = [
    'specimen collector sample ID', 'host specimen voucher', 'sample collected by', 
    'sample collection date', 'sample collection time', 'geolocation name (region)', 
    'CT value', 'geolocation latitude', 'geolocation longitude', 
    'environmental material', 'environmental site', 'collection device', 
    'collection protocol', 'specimen processing method', 'specimen processed by', 
    'library preparation kit', 'library preparation by', 'library ID', 
    'sequencing date', 'sequencing by', 'sequencing instrument', 
    'flowcell barcode', 'sequencing protocol name', 'sequencing kit number', 
    'run folder name', 'raw sequence data processing method', 'r1 fastq filename', 
    'r1 SHA-1 checksum', 'r2 fastq filename', 'r2 SHA-1 checksum', 
    'InterOp file checksum', 'Note'
]

labs = ['western', 'waterloo', 'guelph']

# gather all metadata CSV files and extract sample collection dates
stdout = subprocess.check_output(['find', '/home/wastewater/uploads/', '-name', 'metadata.csv'])
paths = stdout.split()
metadata = {}
for path in paths:
    dirname = os.path.dirname(path).decode()
    
    toks = dirname.split('/')
    lab, runname = toks[-2], toks[-1] 

    if lab not in labs:
        continue  # mapped CSV in wrong location

    if lab not in metadata:
        metadata.update({lab: {}})
    if runname not in metadata[lab]:
        metadata[lab].update({runname: {}})

    handle = open(path, encoding='latin-1')
    for row in csv.DictReader(handle):
        sample_key = 'r1 fastq filename'  #list(row.keys())[0]
        try:
            if len(row) == 1 or row[sample_key] == '':
                # skip blank row, i.e., ",,,,,,,"
                continue
        except:
            print(row)
            raise

        sample = row[sample_key].split('_')[0]

        if lab == 'western':
            sample = sample.replace('_', '-')
        try:
            metadata[lab][runname].update({sample: {
                'coldate': row['sample collection date'],
                'region': row['geolocation name (region)'],
                'latitude': row['geolocation latitude'],
                'longitude': row['geolocation longitude']
            }})
        except:
            print(path)
            continue


# gather all mapped.csv files and import contents, storing lab and run name
stdout = subprocess.check_output(['find', '/home/wastewater/results/', '-name', '*.mapped.csv'])
paths = stdout.split()

outfile = gzip.open("collate-mapped.csv.gz", mode='wt')
writer = csv.writer(outfile)
writer.writerow(["lab", "runname", "sample", "coldate", "region", "latitude",
                 "longitude", "label", "mutation", "freq", "coverage"])

for path in paths:
    filename = os.path.basename(path).decode()
    sample = filename.split('.')[0]
    print(sample)

    # retrieve lab and run name from path
    dirname = os.path.dirname(path).decode()
    
    toks = dirname.split('/')
    lab, runname = toks[-2], toks[-1]
    if lab not in labs:
        continue  # mapped CSV in wrong location

    # troubleshoot failed matches
    if lab not in metadata:
        print(f"lab {lab} not in metadata: {metadata.keys()}")
        continue
    if runname not in metadata[lab]:
        print(f"run {runname} not in {lab} metadata: {metadata[lab].keys()}")
        continue
    if sample not in metadata[lab][runname]:
        print(f"sample {sample} not in {lab}/{runname} metadata: {metadata[lab][runname].keys()}")
        continue

    md = metadata[lab][runname][sample]

    handle = open(path, 'rb')
    for line in handle:
        pos, label, mut, freq, cover = line.decode('latin-1').split(',')
        if freq == 'frequency':
            continue
        freq = float(freq)
        cover = int(cover)
        if cover < 10 or freq < 1e-3:
            continue

        writer.writerow([lab, runname, sample, md['coldate'], md['region'],
                         md['latitude'], md['longitude'], label, mut, freq, cover])

outfile.close()

