from datetime import datetime
from gromstole.scripts.progress_utils import Callback
import sqlite3
import argparse
import shutil
import os
import glob
import subprocess


def open_connection(database, callback=None):
    """
    Open connection to the database. Creates db file if it doesn't exist

    :param database: str, name of the database
    :param callback: function, option to print messages to the console
    :return: objects, connection object and cursor object
    """
    if not os.path.exists(database):
        if callback:
            callback("Database doesn't exist. Creating {}".format(database))

    conn = sqlite3.connect(database, check_same_thread=False)
    cur = conn.cursor()

    # create table if it don't exist
    records_table = 'CREATE TABLE IF NOT EXISTS RECORDS (file VARCHAR(255) ' \
                    'PRIMARY KEY, creation_date DATE, location VARCHAR(255), checksum VARCHAR(255));'
    cur.execute(records_table)

    conn.commit()
    return cur, conn


def insert_record(curr, filepath):
    """
    Inserts the file name, file creation date, path and checksum of the file into the database

    :param curr: object, cursor object
    :param filepath: str, path to the file
    :return: None
    """
    path, file = os.path.split(filepath)

    # Get creation date
    file_timestamp = os.path.getmtime(filepath)
    date_obj_r1 = datetime.fromtimestamp(file_timestamp)
    formatted_date = date_obj_r1.strftime('%Y-%m-%d %H:%M:%S')

    # Calculate SHA-1 Checksum
    stdout = subprocess.getoutput("sha1sum {}".format(filepath))
    checksum = stdout.split(' ')[0]

    # Insert into Database
    curr.execute("INSERT INTO RECORDS (file, creation_date, location, checksum)"
                 " VALUES(?, ?, ?, ?)",
                 [file, formatted_date, path, checksum])


def get_files(curr, paths, ignore_list, callback=None):
    """
    Detects if there are any new data files that have been uploaded by comparing the list of files to those that
    are already in the database

    :param curr: object, cursor object
    :param paths: list, paths to all R1 files in the uploads directory
    :param ignore_list: list, directories that the user would like to avoid processing
    :param callback: function, option to print messages to the console
    :return: list, paths to all files that have not been inserted into the database
    """
    curr.execute("SELECT file FROM RECORDS;")
    results = [x[0] for x in curr.fetchall()]
    unentered = []
    ignore = []
    for file in paths:
        if ignore_file(file, ignore_list):
            ignore.append(file)
            continue
        path, filename = os.path.split(file)
        if filename not in results:
            unentered.append(file)

    if len(ignore) > 0 and callback:
        callback("Ignoring {} files".format(len(ignore)))
        for f in ignore:
            callback("\t {}".format(f), level="DEBUG")

    return unentered


def ignore_file(file, ignore):
    """
    Ignores file if the file is within any of the user specified directories

    :param file: str, path to a file
    :param ignore: list, directories to ignore
    :return: bool, True if the directory is in the file path, False otherwise
    """
    for directory in ignore:
        if directory in file:
            return True
    return False


def missing_file(prefix, suffixes):
    """
    Checks if any of the output files are missing

    :param prefix: str, file prefix
    :param suffixes: list, file extensions
    :return: bool, True if the file is missing, False otherwise
    """
    for suffix in suffixes:
        if not os.path.exists("results/{}.{}".format(prefix, suffix)):
            return True
    return False


def process_files(curr, indir, outdir, paths, binpath="minimap2", cutabin="cutadapt", callback=None):
    """
    Process new R1 and R2 data files and copy output files to the results directory

    :param curr: object, cursor object
    :param indir: str, path to the uploads directory where all the FASTAQ data files are location
    :param outdir: str, path to the results directory
    :param paths: list, paths to all the new files that need to be processed
    :param binpath: str, path to the minimap2 binary executable
    :param cutabin: str, path to the cutadapt binary executable
    :param callback: function, option to print messages to the console
    :return: None
    """
    suffixes = ["coverage.csv", "mapped.csv", "coverage.png", "delta.png"]
    ignore = []
    for r1 in paths:
        r2 = r1.replace('_R1_', '_R2_')

        # ignore files that don't have an R2 file
        if not os.path.exists(r2):
            ignore.append(r1)
            continue

        path, filename = os.path.split(r1)
        prefix = filename.split('_')[0]

        try:
            subprocess.check_call(['python3', 'scripts/minimap2.py', r1, r2,
                                   '-o', 'results',
                                   '-p', prefix,
                                   '-t', '8',
                                   '--ref', 'data/NC_045512.fa',
                                   '-x', binpath])
        except subprocess.CalledProcessError:
            if callback:
                callback("Error running minimap2.py for {} and {}".format(r1, r2), level="ERROR")
            continue

        # Generate plots
        try:
            subprocess.check_call(['Rscript', 'scripts/plot-individual.R',
                                   os.path.join(os.getcwd(), "results"),
                                   prefix])
        except subprocess.CalledProcessError:
            if callback:
                callback("Error generating plot for {}.mapped.csv and {}.coverage.csv".format(prefix, prefix),
                         level="ERROR")
            continue

        result_dir = outdir + path.split(indir)[1]
        os.makedirs("results/{}".format(prefix), exist_ok=True)
        os.makedirs(result_dir, exist_ok=True)

        # Verifies that all the output files exist before moving them to the results directory
        if missing_file(prefix, suffixes):
            if callback:
                callback("Output file missing for {}".format(prefix), level="ERROR")
            continue

        # Rename files with checksum and copy files to the results directory
        for suffix in suffixes:
            filepath = "results/{}.{}".format(prefix, suffix)
            stdout = subprocess.getoutput("sha1sum {}".format(filepath))
            checksum = stdout.split(' ')[0]
            new_filepath = "results/{}.{}.{}".format(prefix, checksum, suffix)
            os.rename(filepath, new_filepath)
            shutil.copy(new_filepath, result_dir)
            shutil.move(new_filepath, "results/{}".format(prefix))

        # Insert file data to the database
        insert_record(curr, r1)
        insert_record(curr, r2)

    if len(ignore) > 0 and callback:
        callback("Ignoring {} files without an R2 pair".format(len(ignore)))
        for f in ignore:
            callback("\t {}".format(f), level="DEBUG")


def parse_args():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description="Runs minimap2.py on new data"
    )
    parser.add_argument('--db', type=str, default="data/gromstole.db",
                        help="Path to the database")
    parser.add_argument('--indir', type=str, default="/data/wastewater/uploads",
                        help="Path to the uploads directory")
    parser.add_argument('--outdir', type=str, default="/data/wastewater/results",
                        help="Path to the results directory")
    parser.add_argument('-i', '--ignore-list', nargs="*", default=[],
                        help="Directories to ignore or keywords to ignore in the filename")
    parser.add_argument('-x', '--binpath', type=str, default='minimap2',
                        help="<option> path to minimap2 executable")
    parser.add_argument('-c', '--cutabin', type=str, default='cutadapt',
                        help="<option> path to cutadapt executable")

    return parser.parse_args()


if __name__ == '__main__':
    cb = Callback()
    args = parse_args()

    cursor, connection = open_connection(args.db, callback=cb.callback)

    files = glob.glob("{}/**/*_R1_*.fastq.gz".format(args.indir), recursive=True)
    new_files = get_files(cursor, files, args.ignore_list, callback=cb.callback)
    process_files(cursor, args.indir, args.outdir, new_files, binpath=args.binpath, cutabin=args.cutabin,
                  callback=cb.callback)

    connection.commit()
    connection.close()
