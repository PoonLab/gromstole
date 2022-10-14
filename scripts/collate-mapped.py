from glob import glob
import subprocess
import os
import sys
import csv
import gzip

# gather all metadata CSV files and extract sample collection dates
stdout = subprocess.check_output(['find', '/data/wastewater/uploads/', '-name', 'metadata.csv'])
paths = stdout.split()
coldates = {}
for path in paths:
    dirname = os.path.dirname(path).decode()
    try:
        _, _, _, _, lab, runname = dirname.split('/')
    except:
        continue  # mapped CSV in wrong location

    if lab not in coldates:
        coldates.update({lab: {}})
    if runname not in coldates[lab]:
        coldates[lab].update({runname: {}})

    handle = open(path, 'rb')
    for line in handle:
        line2 = line.decode('latin-1')  #.encode("utf-8")
        tokens = line2.split(',')
        if tokens[0] == '' or tokens[0] == 'specimen collector sample ID':
            # skip header row or blank row, i.e., ",,,,,,,"
            continue
        sample = tokens[0]
        coldate = tokens[3]
        coldates[lab][runname].update({sample: coldate})


# gather all mapped.csv files and import contents, storing lab and run name
stdout = subprocess.check_output(['find', '/data/wastewater/results/', '-name', '*.mapped.csv'])
paths = stdout.split()

outfile = gzip.open("collate-mapped.csv", mode='wt')

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

    if lab not in coldates or runname not in coldates[lab] or sample not in coldates[lab][runname]:
        print(f"failed to map {lab}/{runname}/{sample} to coldate dict")
        continue
    coldate = coldates[lab][runname][sample]

    handle = open(path, 'rb')
    for line in handle:
        pos, label, mut, freq, cover = line.decode('latin-1').split(',')
        if freq == 'frequency':
            continue
        freq = float(freq)
        cover = int(cover)
        if cover < 10 or freq < 1e-3:
            continue

        outfile.write(f"{lab},{runname},{sample},{coldate},{label},{mut},{freq},{cover}\n")

outfile.close()
