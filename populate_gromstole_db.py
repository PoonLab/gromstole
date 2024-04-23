import argparse
import os
import sys
import glob
import csv
from scripts.progress_utils import Callback
import psycopg2
import psycopg2.extras
from psycopg2 import sql
from psycopg2.errors import DuplicateDatabase


def parse_args():
    parser = argparse.ArgumentParser(description="Populate database with gromstole results.")

    parser.add_argument('--indir', type=str, default="/home/wastewater/results/gromstole",
                        help="Path to the gromstole results directory")
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
                    POSITION VARCHAR(255), 
                    LABEL VARCHAR(255),
                    FREQUENCY VARCHAR(255), 
                    COVERAGE VARCHAR(255),
                    PATH VARCHAR(255))'''
    cur.execute(results_table)

    cur.execute('''CREATE INDEX IF NOT EXISTS frequency_index ON RESULTS (frequency)''')

    conn.commit()
    return cur, conn


def insert_files(cur, files, callback=None):
    """
    Inserts the file name, file creation date, path and checksum of the file into the database

    :param curr: object, cursor object
    :param filepath: str, path to the file
    :return: None
    """
    for lab, run, sample, path in files:
        cur.execute('SELECT * FROM RESULTS WHERE LAB = %s AND RUN = %s AND SAMPLE = %s', (lab, run, sample))
        if cur.fetchone() is not None:
            continue
        callback("Inserting results from file: {}".format(path))
        with open(path, 'r') as f:
            mapped = csv.DictReader(f)
            for line in mapped:
                cur.execute("INSERT INTO RESULTS (LAB, RUN, SAMPLE, POSITION, LABEL, FREQUENCY, COVERAGE, PATH)"
                            " VALUES(%s, %s, %s, %s, %s, %s, %s, %s)"
                            " ON CONFLICT DO NOTHING",
                            (lab, run, sample, line['position'], line['label'], line['frequency'], line['coverage'], path))
    

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
        f.add((lab, run, sample, normfile))
    return f


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
    clean_files = get_files(files)
    insert_files(cur, clean_files, callback=cb.callback)
    conn.commit()
    conn.close()