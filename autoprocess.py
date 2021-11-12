from datetime import datetime
from gromstole.scripts.progress_utils import Callback
import sqlite3
import argparse
import sys
import os
import glob
import subprocess


def open_connection(database, callback=None):
    if not os.path.exists(database):
        if callback:
            callback("Failed to open sqlite3 connection, path {} does not exist".format(database), level="ERROR")
        sys.exit()

    conn = sqlite3.connect(database, check_same_thread=False)
    cur = conn.cursor()

    # create table if it don't exist
    records_table = 'CREATE TABLE IF NOT EXISTS RECORDS (file VARCHAR(255) ' \
                    'PRIMARY KEY, creation_date DATE, location VARCHAR(255), checksum VARCHAR(255));'
    cur.execute(records_table)

    conn.commit()
    return cur, conn


def ignore_file(file, ignore):
    for dir in ignore:
        if dir in file:
            return True
    return False


def get_files(curr, files, ignore_list, callback=None):
    curr.execute("SELECT file FROM RECORDS;")
    results = [x[0] for x in curr.fetchall()]
    new_files = []
    ignored = []
    for file in files:
        if ignore_file(file, ignore_list):
            ignored.append(file)
            continue
        path, filename = os.path.split(file)
        if filename not in results:
            new_files.append(file)

    if len(ignored) > 0 and callback:
        callback("Ignoring {} files".format(len(ignored)))
        for f in ignored:
            callback("\t {}".format(f), level="DEBUG")

    return new_files


def insert_record(curr, filepath):
    path, file = os.path.split(filepath)

    # Get creation date
    file_timestamp = os.path.getmtime(filepath)
    date_obj_r1 = datetime.fromtimestamp(file_timestamp)
    formatted_date = date_obj_r1.strftime('%Y-%m-%d %H:%M:%S')

    # Calculate SHA-1 Checksum
    stdout = subprocess.getoutput("cd {}; sha1sum {}".format(path, file))
    checksum = stdout.split(' ')[0]

    # Insert into Database
    curr.execute("INSERT INTO RECORDS (file, creation_date, location, checksum)"
                 " VALUES(?, ?, ?, ?)",
                 [file, formatted_date, path, checksum])


def process_files(curr, new_files, binpath="minimap2", callback=None):
    ignored = []
    for r1 in new_files:
        r2 = r1.replace('_R1_', '_R2_')
        if not os.path.exists(r2):
            ignored.append(r1)
            continue

        path, filename = os.path.split(r1)
        file_prefix = filename.split('_')[0]

        try:
            subprocess.check_call(['python3', 'scripts/minimap2.py', r1, r2,
                                   '-o', 'results',
                                   '-p', file_prefix,
                                   '-t', '8',
                                   '--ref', 'data/NC_045512.fa',
                                   '-x', binpath])
        except subprocess.CalledProcessError:
            if callback:
                callback("Error running minimap2.py for {} and {}".format(r1, r2), level="ERROR")
            continue

        # Insert file data to the database
        insert_record(curr, r1)
        insert_record(curr, r2)

        # Copy files to the results directory

        # Generate Plot


    if len(ignored) > 0 and callback:
        callback("Ignoring {} files without an R2 pair".format(len(ignored)))
        for f in ignored:
            callback("\t {}".format(f), level="DEBUG")


def parse_args():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
      description="Runs minimap2.py on new data"
    )
    parser.add_argument('--db', type=str, default="data/gromstole.db",
                        help="Path to the database")
    parser.add_argument('--files', type=str, default="/data/wastewater/uploads",
                      help="Path to the uploads directory")
    parser.add_argument('-i', '--ignore-list', nargs="*", default=[],
                        help="Directories to ignore or keywords to ignore in the filename")
    parser.add_argument('-x', '--binpath', type=str, default='minimap2',
                        help="<option> path to minimap2 executable")
    parser.add_argument('-c', '--cbinpath', type=str, default='cutadapt',
                        help="<option> path to cutadapt executable")

    return parser.parse_args()


if __name__ == '__main__':
    cb = Callback()
    args = parse_args()

    cursor, conn = open_connection(args.db, callback=cb.callback)

    files = glob.glob("{}/**/*_R1_*.fastq.gz".format(args.files), recursive=True)
    new_files = get_files(cursor, files, args.ignore_list, callback=cb.callback)
    process_files(cursor, new_files, binpath=args.binpath, callback=cb.callback)

    conn.commit()
    conn.close()