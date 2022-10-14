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

# gather all metadata CSV files and extract sample collection dates
stdout = subprocess.check_output(['find', '/data/wastewater/uploads/', '-name', 'metadata.csv'])
paths = stdout.split()
metadata = {}
for path in paths:
    dirname = os.path.dirname(path).decode()
    try:
        _, _, _, _, lab, runname = dirname.split('/')
    except:
        continue  # mapped CSV in wrong location

    if lab not in metadata:
        metadata.update({lab: {}})
    if runname not in metadata[lab]:
        metadata[lab].update({runname: {}})

    handle = open(path, 'rb')
    for line in handle:
        line2 = line.decode('latin-1')  #.encode("utf-8")
        tokens = line2.split(',')
        if tokens[0] == '' or tokens[0] == 'specimen collector sample ID':
            # skip header row or blank row, i.e., ",,,,,,,"
            continue
        sample = tokens[0]
        metadata[lab][runname].update({sample: {
            'coldate': tokens[3],
            'region': tokens[5],
            'latitude': tokens[7],
            'longitude': tokens[8]
        }})


# gather all mapped.csv files and import contents, storing lab and run name
stdout = subprocess.check_output(['find', '/data/wastewater/results/', '-name', '*.mapped.csv'])
paths = stdout.split()

outfile = gzip.open("collate-mapped.csv.gz", mode='wt')
outfile.write("lab,runname,sample,coldate,region,latitude,longitude,label,"\
              "mutation,freq,coverage\n")

for path in paths:
    filename = os.path.basename(path).decode()
    sample = filename.split('.')[0]
    print(sample)

    # retrieve lab and run name from path
    dirname = os.path.dirname(path).decode()
    try:
        _, _, _, _, lab, runname = dirname.split('/')
    except:
        continue  # mapped CSV in wrong location

    if lab not in metadata or runname not in metadata[lab] or \
        sample not in metadata[lab][runname]:

        print(f"failed to map {lab}/{runname}/{sample} to coldate dict")
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

        outfile.write(
            f"{lab},{runname},{sample},{md['coldate']},{md['region']},"\
            f"{md['latitude']},{md['longitude']},{label},{mut},{freq},{cover}\n"
        )

outfile.close()

