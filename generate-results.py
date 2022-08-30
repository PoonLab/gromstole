from datetime import datetime
from scripts.progress_utils import Callback
import argparse
import shutil
import re
import os
import glob
import subprocess

import smtplib
from dotenv import dotenv_values
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

error_msgs = []

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


def run_scripts(runs, indir, outdir, replace, callback=None):
    """
    Run estimate-freqs.R, make-barplots.R and make-csv.R for each run

    :param runs: list, location of the run files for the estimate-freqs.R script
    :param indir: str, location of the uploads folder
    :param outdir: str, location of the results folder
    :param replace: bool, True if old files need to be replaced, False otherwise
    :param callback: function, option to print messages to the console
    :return: None
    """

    lineages = ['BA.1', 'BA.2', 'B.1.617.2', 'BA.4', 'BA.5', 'BA.2.75', 'BE.1']
    suffixes = ['json', 'barplot.pdf', 'csv']
    for run in runs:
        result_dir = outdir + run.split(indir)[1]
        for lineage in lineages:
            constellation = 'constellations/constellations/definitions/c{}.json'.format(lineage)
            if not os.path.exists(constellation):
                # Initialize submodule if the constellation file cannot be located
                try:
                    subprocess.check_call("git submodule init; git submodule update", shell=True)
                except:
                    callback("Error adding the required submodules")
                    sys.exit()

            filename = "{}-{}-{}".format(os.path.basename(os.path.dirname(result_dir)),
                                         os.path.basename(result_dir), lineage.replace('.',''))
            cmd = ['Rscript', 'scripts/estimate-freqs.R', result_dir, constellation, '{}.json'.format(filename)]

            # Get metadata filename, not all files are named "metadata.csv"
            metadata = 'metadata*'
            for file in os.listdir(run):
                if re.search(metadata, file):
                    cmd.append(os.path.join(run, file))
                    break

            callback("Running '{}'".format(' '.join(cmd)))
            try:
                subprocess.check_call(cmd)
            except subprocess.CalledProcessError:
                if callback:
                    callback("Error running estimate-freqs.R: {}, {}".format(constellation, run),
                             level="ERROR")
                    error_msgs.append("Error running estimate-freqs.R: {}, {}".format(constellation, run))
                continue

            cmd = ['Rscript', 'scripts/make-barplots.R', '{}.json'.format(filename), '{}.barplot.pdf'.format(filename)]
            callback("Running '{}'".format(' '.join(cmd)))
            try:
                subprocess.check_call(cmd)
            except subprocess.CalledProcessError:
                if callback:
                    callback("Error running make-barplots.R: {}, {}".format(constellation, run),
                             level="ERROR")
                    error_msgs.append("Error running make-barplots.R: {}, {}".format(constellation, run))
                continue

            cmd = ['Rscript', 'scripts/make-csv.R', '{}.json'.format(filename), '{}.csv'.format(filename)]
            callback("Running '{}'".format(' '.join(cmd)))
            try:
                subprocess.check_call(cmd)
            except subprocess.CalledProcessError:
                if callback:
                    callback("Error running make-csv.R: {}, {}".format(constellation, run),
                             level="ERROR")
                    error_msgs.append("Error running make-csv.R: {}, {}".format(constellation, run))
                continue

            # Remove prior output files
            if replace:
                resfiles = glob.glob("{}/**/{}.*".format(result_dir, filename), recursive=True)
                for resfile in resfiles:
                    callback("Removing '{}'".format(resfile))
                    os.remove(resfile)

            for suffix in suffixes:
                stdout = subprocess.getoutput("sha1sum {}.{}".format(filename, suffix))
                checksum = stdout.split(' ')[0][:10]
                new_filename = "{}.{}.{}".format(filename, checksum, suffix)
                os.rename('{}.{}'.format(filename, suffix), new_filename)
                shutil.copy(new_filename, result_dir)
                shutil.move(new_filename, "results/{}/{}".format(os.path.basename(result_dir), new_filename))


def parse_args():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description="Runs R scripts to generate JSON, barplots and csv files"
    )
    parser.add_argument('--indir', type=str, default="/data/wastewater/uploads",
                        help="Path to the uploads directory")
    parser.add_argument('--outdir', type=str, default="/data/wastewater/results",
                        help="Path to the results directory")
    parser.add_argument('-i', '--ignore-list', nargs="*", default=[],
                        help="Directories to ignore or keywords to ignore in the filename")
    parser.add_argument('-e', '--email', type=str, default=None,
                        help="<option> recipient email address to send error messages")
    parser.add_argument('--runs', nargs="*", default=[],
                        help="Runs to process again")
    parser.add_argument('--dont-replace', dest='replace', action='store_false',
                        help="Do not remove previous result files")
    parser.add_argument('--replace', dest='replace', action='store_true',
                        help="Do not remove previous result files")
    parser.set_defaults(replace=True)

    return parser.parse_args()


if __name__ == '__main__':
    cb = Callback()
    args = parse_args()
    cb.callback("Starting script")

    # update submodules
    try:
        subprocess.check_call("git submodule foreach git pull origin main", shell=True)
    except:
        cb.callback("Could not update submodules", level='ERROR')

    files = glob.glob("{}/**/*_R1_*.fastq.gz".format(args.indir), recursive=True)
    runs = get_runs(files, args.ignore_list, args.runs, callback=cb.callback)
    run_scripts(runs, args.indir, args.outdir, args.replace, callback=cb.callback)

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

    cb.callback("All Done!")