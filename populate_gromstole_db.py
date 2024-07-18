import argparse
import os
import sys
import glob
import csv
from datetime import datetime
from epiweeks import Week
import json
from scripts.progress_utils import Callback
import psycopg2
import psycopg2.extras
from psycopg2 import sql
from psycopg2.errors import DuplicateDatabase


def parse_args():
    parser = argparse.ArgumentParser(description="Populate database with gromstole results.")

    parser.add_argument('--indir', type=str, default="/home/wastewater/results/gromstole",
                        help="Path to the gromstole results directory")
    parser.add_argument('--upload_dir', type=str, default="/home/wastewater/uploads",
                        help="Path to the uploads directory")
    parser.add_argument('--dbname', type=str, default=os.environ.get("POSTGRES_DB", "gromstole_db"),
                        help="Postgresql database name")
    parser.add_argument('--dbhost', type=str, default=os.environ.get("POSTGRES_HOST", "localhost"),
                        help="Postgresql database host address")
    parser.add_argument('--dbport', type=str, default=os.environ.get("POSTGRES_PORT", "5432"),
                        help="Connection to port number")
    parser.add_argument('--dbuser', type=str, default=os.environ.get("POSTGRES_USER", None),
                        help="Postgresl user")
    parser.add_argument('--dbpswd', type=str, default=os.environ.get("POSTGRES_PASSWORD", None),
                        help="Postgresl password")
    
    return parser.parse_args()


def open_connection(connection_parameters):
    """ open connection to database, initialize tables if they don't exist
        :out:
        :cursor: interactive sql object containing tables
    """
    conn = psycopg2.connect(**connection_parameters)
    cur = conn.cursor(cursor_factory = psycopg2.extras.RealDictCursor)

    # create tables if they don't exist
    results_table = '''CREATE TABLE IF NOT EXISTS RESULTS (
                    ID SERIAL PRIMARY KEY, 
                    LAB VARCHAR(255), 
                    RUN VARCHAR(255),
                    SAMPLE VARCHAR(255),
                    COLDATE DATE,
                    REGION VARCHAR(255),
                    LATITUDE FLOAT NULL,
                    LONGITUDE FLOAT NULL,
                    LABEL VARCHAR(255),
                    MUTATION VARCHAR(255),
                    FREQUENCY VARCHAR(255), 
                    COVERAGE VARCHAR(255),
                    PATH VARCHAR(255))'''
    cur.execute(results_table)

    cur.execute('''CREATE INDEX IF NOT EXISTS frequency_index ON RESULTS (frequency)''')

    aggregate_table = '''CREATE TABLE IF NOT EXISTS AGGREGATE_MAPPED (
                    ID SERIAL PRIMARY KEY,
                    REGION VARCHAR(255),
                    YEAR VARCHAR(255),
                    EPIWEEK VARCHAR(255),
                    NUC VARCHAR(255),
                    AMINO VARCHAR(255),
                    NSAMPLES INTEGER,
                    COUNT INTEGER,
                    COVERAGE INTEGER )'''
    cur.execute(aggregate_table)

    conn.commit()
    return cur, conn


def parse_date(dt, formats=('%y-%m-%d', '%d-%m-%y', '%d-%m-%Y', '%m-%d-%Y', '%Y-%m-%d', '%d-%b-%y', '%m/%d/%Y', '%d/%m/%y')):
    """ Try multiple date formats """
    if dt in ['UNK', '', ' ', 'null', None]:
        return None
    for fmt in formats:
        try:
            return datetime.strptime(dt, fmt)
        except ValueError:
            pass
    raise ValueError(f"No supported format detected for date string '{dt}'")


def retrieve_metadata(runs, upload_dir, callback=None):
    """
    Retrieve metadata from the metadata.csv file

    :param runs: set, runs to retrieve metadata for
    :param upload_dir: str, path to the upload directory
    :return: dict, metadata
    """
    metadata = {}
    for lab, run in runs:
        if lab not in metadata:
            metadata.update({lab: {}})
        if run not in metadata[lab]:
            metadata[lab].update({run: {}})
        
        upload_path = os.path.join(upload_dir, lab, run, 'metadata.csv')
        if not os.path.exists(upload_path):
            callback('Metadata file not found: {}'.format(upload_path))
            continue

        with open(upload_path, encoding='latin-1') as handle:
            for row in csv.DictReader(handle):
                sample_key = 'r1 fastq filename'
                try:
                    if len(row) == 1 or row[sample_key] == '':
                        continue
                except:
                    callback('Empty metadata file: {}'.format(row))
                    raise

                sample = row[sample_key].split('_')[0]

                if lab == 'western':
                    sample = sample.replace('_', '-')
                try:
                    metadata[lab][run].update({sample: {
                        'coldate': row['sample collection date'],
                        'region': row['geolocation name (region)'],
                        'latitude': row['geolocation latitude'],
                        'longitude': row['geolocation longitude']
                    }})
                except:
                    callback('Column names incorrect in {}'.format(os.path.join(upload_path)))
                    continue
    return metadata


def new_mapped_files(cur, files, callback=None):
    """
    Get the set of new mapped files to process

    :param curr: object, cursor object
    :param filepath: str, path to the file
    :return: set, new files, new runs
    """
    new_files = set()
    new_runs = set()
    for lab, run, sample, path in files:
        cur.execute('SELECT * FROM RESULTS WHERE LAB = %s AND RUN = %s AND SAMPLE = %s', (lab, run, sample))
        if cur.fetchone() is not None:
            continue
        new_files.add((lab, run, sample, path))
        new_runs.add((lab, run))
    return new_files, new_runs
    

def get_files(files):
    f = set()
    labs = ['western', 'waterloo', 'guelph']
    for file in files:
        normfile = os.path.normpath(file)
        found_keywords = [keyword for keyword in labs if keyword in normfile]
        if len(found_keywords) == 0:
            print("No lab found in file path: {}".format(normfile))
            continue
        runpath, sample = os.path.split(normfile)
        run = os.path.basename(runpath)
        lab = found_keywords[0]
        f.add((lab, run, sample.split('.')[0], normfile))
    return f


def insert_files(cur, files, metadata, callback=None):
    """ 
    Insert files into the database

    :param curr: object, cursor object
    :param files: set, files to insert
    :param metadata: dict, metadata
    :return: None
    """

    # import location to region map as a dict
    # TODO: Need a better way to map locations to regions
    regions = json.load(open("data/regions.json"))

    for lab, run, sample, path in files:
        if lab not in metadata:
            callback("Lab {} not in metadata".format(lab))
            continue
        if run not in metadata[lab]:
            callback("Run {} not in metadata".format(run))
            continue
        if sample not in metadata[lab][run]:
            callback("Sample {} from {}/{} not in metadata".format(sample, lab, run))
            continue

        md = metadata[lab][run][sample]
        callback("Inserting results from file: {}".format(path))
        with open(path, 'r') as f:
            mapped = csv.DictReader(f)
            last_fail = None
            failed_latlong = set()
            for line in mapped:
                coverage = float(line['coverage'])
                frequency = float(line['frequency'])
                if coverage < 10 or frequency < 1e-3:
                    continue
                
                coldate = parse_date(md['coldate'].strip())
                if coldate is None:
                    continue

                try:
                    latitude = float(md['latitude'])
                    longitude = float(md['longitude'])
                except ValueError:
                    if (lab, run, sample) not in failed_latlong:
                        failed_latlong.add((lab, run, sample))
                        callback(f"Setting latitude and longitude to None for {lab}/{run}/{sample}")
                    latitude = None
                    longitude = None

                cur.execute("INSERT INTO RESULTS (LAB, RUN, SAMPLE, COLDATE, REGION, LATITUDE, LONGITUDE, LABEL, MUTATION, FREQUENCY, COVERAGE, PATH)"
                            " VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
                            " ON CONFLICT DO NOTHING",
                            (lab, run, sample, coldate, md['region'], latitude, longitude,
                              line['position'], line['label'], line['frequency'], line['coverage'], path))

                
                # Aggregate Data
                location = md['region'].strip()
                region = regions.get(location, None)
                if region is None and location != last_fail:
                    callback(f"Failed to map location {location} to region")
                    last_fail = location
                    continue

                epiweek = Week.fromdate(coldate.date())
                mutation = (line['label'], line['mutation'])

                # Check to see if the value exists in database
                cur.execute("SELECT * FROM AGGREGATE_MAPPED WHERE REGION = %s AND YEAR = '%s' AND EPIWEEK = '%s' AND NUC = %s AND AMINO = %s",
                            (region, epiweek.year, epiweek.week, mutation[0], mutation[1]))
                if cur.fetchone() is None:
                    cur.execute("INSERT INTO AGGREGATE_MAPPED (REGION, YEAR, EPIWEEK, NUC, AMINO, NSAMPLES, COUNT, COVERAGE)"
                                " VALUES(%s, %s, %s, %s, %s, %s, %s, %s)",
                                (region, epiweek.year, epiweek.week, mutation[0], mutation[1], 1, frequency*coverage, coverage))
                else:
                    cur.execute("UPDATE AGGREGATE_MAPPED SET NSAMPLES = NSAMPLES + 1, COUNT = COUNT + %s, COVERAGE = COVERAGE + %s"
                                " WHERE REGION = %s AND YEAR = '%s' AND EPIWEEK = '%s' AND NUC = %s AND AMINO = %s",
                                (frequency*coverage, coverage, region, epiweek.year, epiweek.week, mutation[0], mutation[1]))


if __name__ == "__main__":
    args = parse_args()
    cb = Callback()

    # Check if database exists
    connection_parameters = {
        "host": args.dbhost,
        "port": args.dbport,
        "user": args.dbuser,
        "password": args.dbpswd,
    }

    connection = None
    try:
        connection = psycopg2.connect(**connection_parameters)
        connection.autocommit = True

        cursor = connection.cursor()
        cursor.execute(sql.SQL('CREATE DATABASE {}').format(sql.Identifier(args.dbname)))
        cb.callback("Database {} created successfully.".format(args.dbname))
    except DuplicateDatabase:
        cb.callback("Database {} already exists.".format(args.dbname))
    except psycopg2.Error as e:
        cb.callback("Error initiating connection to database: {}".format(e))
        sys.exit()
    finally:
        if connection is not None:
            cursor.close()
            connection.close()

    connection_parameters['dbname'] = args.dbname
    cur, conn = open_connection(connection_parameters)
    files = glob.glob("{}/**/*.mapped.csv".format(args.indir), recursive=True)
    mapped_files = get_files(files)
    unprocessed_mapped_files, unprocessed_runs = new_mapped_files(cur, mapped_files, callback=cb.callback)
    metadata = retrieve_metadata(unprocessed_runs, args.upload_dir, callback=cb.callback)
    insert_files(cur, unprocessed_mapped_files, metadata, callback=cb.callback)
    conn.commit()
    conn.close()