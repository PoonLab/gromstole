from datetime import datetime
from progress_utils import Callback
from trim import cutadapt, minimap2, freyja
import sqlite3
import argparse
import shutil
import re
import os
import glob
import subprocess
import sys

import smtplib
from dotenv import dotenv_values
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

error_msgs = []

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
    curr.execute("INSERT OR REPLACE INTO RECORDS (file, creation_date, location, checksum)"
                 " VALUES(?, ?, ?, ?)",
                 [filepath, formatted_date, path, checksum])


def get_runs(paths, ignore_list, runs_list, callback=None):
    """
    Gets the list of runs to process

    :param paths: list, paths to all R1 files in the uploads directory
    :param ignore_list: list, directories that the user would like to avoid processing
    :param runs_list: list, directories that the user would like to processing, processes everything if none specified
    :param callback: function, option to print messages to the console
    :return: list, paths to all files that have not been inserted into the database and the filepaths to all runs
    """
    runs = set()
    ignored = set()
    for file in paths:
        path, filename = os.path.split(file)
        if contains_run(path, ignore_list):
            ignored.add(path)
            continue
        if len(runs_list) == 0 or contains_run(path, runs_list):
            runs.add(path)

    for path in ignored:
        callback("Ignoring '{}'".format(path))

    return runs


def get_files(curr, paths, ignore_list, check_processed, callback=None):
    """
    Detects if there are any new data files that have been uploaded by comparing the list of files to those that
    are already in the database

    :param curr: object, cursor object
    :param paths: list, paths to all R1 files in the uploads directory
    :param ignore_list: list, directories that the user would like to avoid processing
    :param callback: function, option to print messages to the console
    :return: list, paths to all files that have not been inserted into the database and the filepaths to all runs
    """
    runs = set()
    curr.execute("SELECT file, checksum FROM RECORDS;")
    results = {x[0]: x[1] for x in curr.fetchall()}

    unentered = []
    entered = []
    ignore = []
    for file in paths:
        if contains_run(file, ignore_list):
            ignore.append(file)
            continue
        path, filename = os.path.split(file)
        if file not in results:
            unentered.append(file)
            runs.add(path)
        else:
            entered.append(file)

    if check_processed:
        # Check if R1 or R2 files have changed since the last run
        if len(entered) > 0:
            r2_files = [file.replace('_R1_', '_R2_') for file in entered]
            stdout = (subprocess.getoutput("sha1sum {} {}".format(' '.join(entered), ' '.join(r2_files)))).split()
            checksums = {stdout[i]: stdout[i - 1] for i in range(1, len(stdout), 2)}

        for file in entered:
            path, filename = os.path.split(file)
            r2 = file.replace('_R1_', '_R2_')
            _, r2_filename = os.path.split(r2)

            if results[filename] != checksums[file] or results[r2_filename] != checksums[r2]:
                unentered.append(file)
                runs.add(path)
                cb.callback("Re-processing {} and its R2 pair from the database".format(file))

    if len(ignore) > 0 and callback:
        callback("Ignoring {} files".format(len(ignore)))
        for f in ignore:
            callback("\t {}".format(f), level="DEBUG")

    return unentered, runs


def contains_run(run, usr_list):
    """
    Ignores/Includes run if its in the user specified directories

    :param run: str, path to a run
    :param usr_list: list, directories to ignore
    :return: bool, True if the directory is in the file path, False otherwise
    """
    for directory in usr_list:
        if directory in run:
            return True
    return False


def process_files(curr, indir, outdir, paths, binpath="minimap2", cutabin="cutadapt", fbin="freyja", threads=4, ref="data/NC_045512.fa", callback=None):
    """
    Process new R1 and R2 data files and copy output files to the results directory

    :param curr: object, cursor object
    :param indir: str, path to the uploads directory where all the FASTAQ data files are location
    :param outdir: str, path to the results directory
    :param paths: list, paths to all the new files that need to be processed
    :param binpath: str, path to the minimap2 binary executable
    :param cutabin: str, path to the cutadapt binary executable
    :param fbin: str, path to the freyja binary executable
    :param threads: int, number of threads to run minimap2
    :param ref: str, path to the reference sequence
    :param callback: function, option to print messages to the console
    :return: None
    """

    ignore = []
    for r1 in paths:
        r2 = r1.replace('_R1_', '_R2_')

        # ignore files that don't have an R2 file
        if not os.path.exists(r2):
            ignore.append(r1)
            continue

        path, filename = os.path.split(r1)
        prefix = filename.split('_')[0]

        os.makedirs(os.path.join(os.getcwd(), "results"), exist_ok=True)

        if callback:
            callback("starting {} from {}".format(prefix, path))

        result_dir = os.path.join(outdir + path.split(indir)[1], prefix)
        local_results_dir = os.path.join("results" + path.split(indir)[1], prefix)
        os.makedirs(local_results_dir, exist_ok=True)
        os.makedirs(result_dir, exist_ok=True)

        tf1, tf2 = cutadapt(fq1=r1, fq2=r2, ncores=2, path=cutabin)
        bamsort = minimap2(fq1=tf1, fq2=tf2, ref=ref, nthread=threads,
                            path=binpath)
        frey = freyja(bamsort=bamsort, ref=ref, sample=prefix, outpath=local_results_dir,
                        path=fbin)

        # Remove cutadapt output temporary files
        for tmpfile in [tf1, tf2, bamsort]:
            cb.callback("Removing {}".format(tmpfile))
            try:
                os.remove(tmpfile)
            except FileNotFoundError:
                pass
        
        # Rename files with checksum and move to the output directory
        for f in os.listdir(local_results_dir):
            current_filepath = os.path.join(local_results_dir, f)
            stdout = subprocess.getoutput("sha1sum {}".format(current_filepath))
            checksum = stdout.split(' ')[0][:10]
            filename_toks = f.split('.')
            new_filename = '{}.{}.{}'.format('.'.join(filename_toks[:-1]), checksum,filename_toks[-1])
            new_filepath = os.path.join(local_results_dir, new_filename)
            os.rename(current_filepath, new_filepath)
            shutil.copy(new_filepath, result_dir)

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
    parser.add_argument('--indir', type=str, default="/home/wastewater/uploads",
                        help="Path to the uploads directory")
    parser.add_argument('--outdir', type=str, default="/home/wastewater/results",
                        help="Path to the results directory")
    parser.add_argument('-i', '--ignore-list', nargs="*", default=[],
                        help="Directories to ignore or keywords to ignore in the filename")
    parser.add_argument('--minimap2', type=str, default='minimap2',
                        help="<option> path to minimap2 executable")
    parser.add_argument('--cutadapt', type=str, default='cutadapt',
                        help="<option> path to cutadapt executable")
    parser.add_argument('--freyja', type=str, default="freyja",
                        help="<option> path to freyja executable")
    parser.add_argument('--threads', type=int, default=4,
                        help="Number of threads (default 4) for minimap2")
    parser.add_argument('-e', '--email', type=str, default=None,
                        help="<option> recipient email address to send error messages")                        
    parser.add_argument('--check', dest='check', action='store_true',
                        help="Check processed files to see if the files have changed since last processing them")
    parser.add_argument('--replace', dest='replace', action='store_false',
                        help="Remove previous result files")
    parser.add_argument('--rerun', dest='rerun', action='store_true',
                        help="Process result files again (generate JSON, csv and barplots again)")
    parser.add_argument('--runs', nargs="*", default=[],
                        help="Runs to process again")
    parser.set_defaults(check=False, replace=True, rerun=False)

    return parser.parse_args()


if __name__ == '__main__':
    cb = Callback()
    args = parse_args()
    cb.callback("Starting script")

    files = glob.glob("{}/**/*_R1_*.fastq.gz".format(args.indir), recursive=True)

    if args.rerun:
        runs = get_runs(files, args.ignore_list, args.runs, callback=cb.callback)
    else:
        cursor, connection = open_connection(args.db, callback=cb.callback)
        new_files, runs = get_files(cursor, files, args.ignore_list, args.check, callback=cb.callback)
        process_files(cursor, args.indir, args.outdir, new_files, binpath=args.minimap2, cutabin=args.cutadapt,
                      callback=cb.callback)
        connection.commit()
        connection.close()  

    # Generate result files
    # run_scripts(runs, args.indir, args.outdir, args.replace, callback=cb.callback)

    if len(error_msgs) > 0 and args.email is not None:
        cb.callback("Sending error message")

        config = dotenv_values(".env")
        try:
            server = smtplib.SMTP_SSL(config["HOST"], int(config["PORT"]))
            server.ehlo()
            server.login(config["EMAIL_ADDRESS"], config["EMAIL_PASSWORD"])
        except:
            cb.callback("There was a problem initializing a connection with the server")
            exit(-1)

        msg = MIMEMultipart("related")
        msg['Subject'] = "ATTENTION: Error running gromstole pipeline"
        msg['From'] = "Gromstole Notification <{}>".format(config["EMAIL_ADDRESS"])
        msg['To'] = args.email

        body = '\r\n'.join(error_msgs)
        msg.attach(MIMEText(body, 'plain'))
        server.sendmail(config["EMAIL_ADDRESS"], args.email, msg.as_string())
        server.quit()
    
    if not args.rerun and len(new_files) == 0:
        cb.callback("No new data files")
    else:
        cb.callback("All Done!")     
