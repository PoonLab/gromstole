import subprocess
import argparse
import sys
import os
import tempfile
import fnmatch
from progress_utils import Callback

import smtplib
import freyja as freyja_module
from dotenv import dotenv_values
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

def send_error_notification(message): 
    """
    Send an email notification to a user
    :param message:  str, message to send
    :return:
    
    """
    config = dotenv_values(".env")
    try:
        server = smtplib.SMTP_SSL(config["HOST"], int(config["PORT"]))
        server.ehlo()
        server.login(config["EMAIL_ADDRESS"], config["EMAIL_PASSWORD"])
    except:
        sys.stderr.write(f"There was a problem initializing a connection with the server")
        exit(-1)

    msg = MIMEMultipart("related")
    msg['Subject'] = "ATTENTION: Error running freyja pipeline"
    msg['From'] = "Gromstole Notification <{}>".format(config["EMAIL_ADDRESS"])
    msg['To'] = config["SEND_TO"]

    body = '{}\r\n'.format(message)
    msg.attach(MIMEText(body, 'plain'))
    server.sendmail(config["EMAIL_ADDRESS"], config["SEND_TO"], msg.as_string())
    server.quit()


def cutadapt(fq1, fq2, path="cutadapt", adapter="AGATCGGAAGAGC", ncores=1, sendemail=False, minlen=10):
    """
    Wrapper for cutadapt
    :param fq1:  str, path to FASTQ R1 file
    :param fq2:  str, path to FASTQ R2 file
    :param adapter:  adapter sequence, defaults to universal Illumina TruSeq
    :param ncores:  int, number of cores to run cutadapt
    :param minlen:  int, discard trimmed reads below this length
    :return:
    """
    of1 = tempfile.NamedTemporaryFile('w', delete=False)
    of2 = tempfile.NamedTemporaryFile('w', delete=False)
    cmd = [path, '-a', adapter, '-A', adapter, '-o', of1.name, '-p', of2.name,
           '-j', str(ncores), '-m', str(minlen), '--quiet', fq1, fq2]
    
    try:
        _ = subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        sys.stderr.write(f"Error running cutadapt command:\n")
        sys.stderr.write(f"{' '.join(cmd)}\n")
        if sendemail:
            send_error_notification(message='Error running {}'.format(' '.join(cmd)))
        raise

    of1.close()
    of2.close()
    return of1.name, of2.name


def minimap2(fq1, fq2, ref, path='minimap2', nthread=1, sendemail=False):
    """
    Wrapper for minimap2 for processing paired-end read data.

    :param fq1:  str, path to FASTQ R1 file (uncompressed or gzipped)
    :param fq2:  str, path to FASTQ R2 file (uncompressed or gzipped); if None,
                 treat fq1 as an unpaired (single) FASTQ
    :param ref:  str, path to FASTA file with reference sequence
    :param path:  str, path to minimap2 binary executable
    :param nthread:  int, number of threads to run

    :return:  str, path to sorted BAM file
    """
    tempsam = tempfile.NamedTemporaryFile('w', delete=False)
    tempbam = tempfile.NamedTemporaryFile('w', delete=False)
    tempbamsort = tempfile.NamedTemporaryFile('w', delete=False)

    # assume we are working with paired-end Illumina
    cmd = [path, '-t', str(nthread), '-ax', 'sr', '--eqx', '--secondary=no', ref, fq1, fq2, '>', tempsam.name]
    os.system(' '.join(cmd))
    tempsam.close()

    cmd = ['samtools', 'view', '-S', '-b', tempsam.name, '>', tempbam.name]
    os.system(' '.join(cmd))
    tempbam.close()

    cmd = ['samtools', 'sort', tempbam.name, '-o', tempbamsort.name]
    os.system(' '.join(cmd))
    tempbamsort.close()

    # remove temporary files
    for tmpfile in [tempsam.name, tempbam.name]:
        try:
            os.remove(tmpfile)
        except FileNotFoundError:
            sys.stderr.write(f"Failed to remove file {tmpfile}, file not found")
            if sendemail:
                send_error_notification(message='Check logs. Could not locate temporary file')
            raise

    return tempbamsort.name 


def freyja(bamsort, ref, sample, outpath, email=None, path='freyja', sendemail=False, callback=None):
    """
    Wrapper function for andersen-lab/Freyja

    :param bamsort:  str, path to sorted BAM file
    :param ref:  str, path to reference FASTA
    :param sample:  str, sample name
    :param outpath:  str, path to write outputs
    :param path:  str, path to Freyja executable
    """
    cmd = [path, 'variants', bamsort, '--variants', '{}/var.{}'.format(outpath, sample),
           '--depths', '{}/depth.{}.csv'.format(outpath, sample), '--ref', ref]
    try:
        _ = subprocess.check_call(cmd)
    except FileNotFoundError:
        sys.stderr.write(f"Failed to find executable file at {path}\n")
        if sendemail:
            send_error_notification(message='Check error logs. Could not find {}'.format(path))
        raise
    except subprocess.CalledProcessError:
        sys.stderr.write(f"Error running Freyja variants command:\n")
        sys.stderr.write(f"{' '.join(cmd)}\n")
        if sendemail:
            send_error_notification(message='Error running {}'.format(' '.join(cmd)))
        raise

    # NOTE: Temporarily disabiling Bootstrap step
    # bootstrap = [path, 'boot', '{}/var.{}.tsv'.format(outpath, sample),
    #              '{}/depth.{}.csv'.format(outpath, sample),
    #              '--nt', '4', '--nb', '10', '--output_base',
    #              '{}/boostrap.{}'.format(outpath, sample), '--boxplot', 'pdf']
    # try:
    #     _ = subprocess.check_call(bootstrap)
    # except subprocess.CalledProcessError:
    #     sys.stderr.write(f"Error running Freyja boot command:")
    #     sys.stderr.write(' '.join(bootstrap))
    #     raise

    demix = [path, 'demix', '{}/var.{}.tsv'.format(outpath, sample),
             '{}/depth.{}.csv'.format(outpath, sample),
             '--output', '{}/lin.{}.tsv'.format(outpath, sample)]

    p = subprocess.Popen(demix, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdoutput, error = p.communicate()

    if p.returncode != 0: 
        if b'cvxpy.error.SolverError' in error:
            sys.stderr.write(f"cvxpy.error.SolverError when running Freyja demix. Check coverage - ")
            sys.stderr.write(f"{' '.join(demix)}\n")
        elif b'demix: Solver error' in stdoutput:
            sys.stderr.write(f"cvxpy.error.SolverError when running Freyja demix. Check coverage - ")
            sys.stderr.write(f"{' '.join(demix)}\n")
        else:
            if sendemail:
                send_error_notification(message="Error running {}".format(' '.join(demix)))
            sys.stderr.write(f"Error running Freyja demix command - ")
            sys.stderr.write(f"{' '.join(demix)}\n")
            raise Exception(f"Subprocess failed with return code {p.returncode} when running {' '.join(demix)}")

    p.stdout.close()
    p.kill()

    return True if p.returncode != 0 else False


def rename_file(filepath):
    """
    Rename a file to include a checksum of the file contents

    :param filepath:  str, path to file to rename
    :return:  str, path to renamed file
    """

    filename = os.path.basename(filepath)
    stdout = subprocess.getoutput("sha1sum {}".format(filepath))
    checksum = stdout.split(' ')[0][:10]
    filename_toks = filename.split('.')
    return filepath.replace(filename, '{}.{}.{}'.format('.'.join(filename_toks[:-1]), checksum,filename_toks[-1]))


if __name__ == '__main__':
    # command line interface
    parser = argparse.ArgumentParser(
        description="Run Freyja Pipeline"
    )
    parser.add_argument('paths', type=argparse.FileType('r'),
                        help="Path to file with list of files to process")
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help="Path to file to write processed files")
    parser.add_argument('--outdir', type=str, default="/home/wastewater/results/freyja",
                        help="Path to the results directory")
    parser.add_argument('--indir', type=str, default="/home/wastewater/uploads",
                        help="Path to the uploads directory")
    parser.add_argument('--ref', type=str, default=os.path.join(freyja_module.__path__[0], "data", "NC_045512_Hu-1.fasta"),
                        help="<input> path to target FASTA (reference)")
    parser.add_argument('--sendemail', dest='sendemail', action="store_true",
                        help="<option> send email notification when there is an error") 
    parser.add_argument('--cutadapt', type=str, default="cutadapt",
                        help="Path to cutadapt")
    parser.add_argument('--minimap2', type=str, default="minimap2",
                        help="Path to minimap2")
    parser.add_argument('--freyja', type=str, default="freyja",
                        help="Path to Freyja")
    parser.add_argument('--threads', type=int, default=4,
                        help="Number of threads (default 4)")
    parser.set_defaults(sendemail=False)

    args = parser.parse_args()

    cb = Callback()
    cb.callback("Starting Script")

    try:
        from mpi4py import MPI
    except ModuleNotFoundError:
        print("Script requires mpi4py - https://pypi.org/project/mpi4py/")
        sys.exit()
    
    args.indir = os.path.normpath(args.indir)
    args.outdir = os.path.normpath(args.outdir)

    comm = MPI.COMM_WORLD
    my_rank = comm.Get_rank()
    nprocs = comm.Get_size()

    r1_paths = []
    # with open(args.paths, 'r') as handle:
    for line in args.paths:
        r1_paths.append(os.path.normpath(line.strip()))
    
    error_messages = []
    records = []
    for num, r1 in enumerate(r1_paths):
        if num % nprocs != my_rank:
            continue
        
        r2 = r1.replace('_R1_', '_R2_')

        # ignore files that don't have an R2 file
        if not os.path.exists(r2):
            sys.stderr.write(f"Missing R2 file for {r1}\n")
            send_error_notification("Missing R2 file for {}".format(r1))
            continue

        path, filename = os.path.split(r1)
        prefix = filename.split('_')[0]

        # os.makedirs(os.path.join(os.getcwd(), "results"), exist_ok=True)

        cb.callback("starting {} from {}".format(prefix, path))

        result_dir = os.path.join(args.outdir + path.split(args.indir)[1], prefix)
        # local_results_dir = os.path.join("results" + path.split(args.indir)[1], prefix)
        # os.makedirs(local_results_dir, exist_ok=True)
        os.makedirs(result_dir, exist_ok=True)

        # Check to see if the results directory has freyja run data
        if len(fnmatch.filter(os.listdir(result_dir), '*.*')) > 0:
            cb.callback("Files already exist in the folder")
            continue

        tf1, tf2 = cutadapt(fq1=r1, fq2=r2, ncores=2, path=args.cutadapt, sendemail=args.sendemail)
        bamsort = minimap2(fq1=tf1, fq2=tf2, ref=args.ref, nthread=args.threads,
                            path=args.minimap2, sendemail=args.sendemail)
        res = freyja(bamsort=bamsort, ref=args.ref, sample=prefix, outpath=result_dir,
                    path=args.freyja, sendemail=args.sendemail, callback=cb.callback)
        
        # cvxpy error
        if res:
            error_messages.append(result_dir)

        # Remove cutadapt output temporary files
        for tmpfile in [tf1, tf2, bamsort]:
            cb.callback("Removing {}".format(tmpfile))
            try:
                os.remove(tmpfile)
            except FileNotFoundError:
                pass

        # Rename files with checksum and move to the output directory
        for f in os.listdir(result_dir):
            current_filepath = os.path.join(result_dir, f)
            new_filepath = rename_file(current_filepath)
            os.rename(current_filepath, new_filepath)

        records.append(r1)
        records.append(r2)


    comm.Barrier()
    all_records = comm.gather(records, root=0)
    all_errors = comm.gather(error_messages, root=0)

    if my_rank == 0: 
        gather_records = [r for record in all_records for r in record]
        gather_errors = [e for errs in all_errors for e in errs]
        for err in gather_errors:
            sys.stderr.write(f"cvxpy error running Freyja demix command - {err} \n")

        for rec in gather_records:
            args.outfile.write(f"{rec}\n")
    
