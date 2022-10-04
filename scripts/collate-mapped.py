from glob import glob
import subprocess
import os
import csv

# gather all mapped.csv files and import contents, storing lab and run name
stdout = subprocess.check_output(['find', '/data/wastewater/results/', '-name', '*.mapped.csv'])
paths = stdout.split()

mapped = {}
for path in paths:
    filename = os.path.basename(path).decode()
    sample = filename.split('.')[0]
    # retrieve lab and run name from path
    dirname = os.path.dirname(path).decode()
    try:
        _, _, _, _, lab, runname = dirname.split('/')
    except:
        continue  # mapped CSV in wrong location

    if lab not in mapped:
        mapped.update({lab: {}})
    if runname not in mapped[lab]:
        mapped[lab].update({runname: {}})
    mapped[lab][runname].update({sample: {'coldate': None, 'mutations': {}}})

# gather all metadata CSV files and extract sample collection dates
stdout = subprocess.check_output(['find', '/data/wastewater/uploads/', '-name', 'metadata.csv'])
paths = stdout.split()

for path in paths:
    dirname = os.path.dirname(path).decode()
    try:
        _, _, _, _, lab, runname = dirname.split('/')
    except:
        continue  # mapped CSV in wrong location

    if runname not in mapped[lab]:
        print(f"failed to locate {runname}")
        continue

    handle = open(path, 'rb')
    for line in handle:
        line2 = line.decode('latin-1')  #.encode("utf-8")
        tokens = line2.split(',')
        if tokens[0] == '' or tokens[0] == 'specimen collector sample ID':
            continue
        sample = tokens[0]
        coldate = tokens[3]
