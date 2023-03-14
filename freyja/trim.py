import csv
import subprocess
import argparse
import re
import sys
import os
import tempfile
import pandas as pd
import glob
from progress_utils import Callback


def generate_csv(outpath):
    """

    :param outpath:  str, path to write output TSV files
    :return:
    """
    tsvfiles = glob.glob("{}/**/lin.*.tsv".format(outpath), recursive=True)
    results = {}
    lineage_set = set()
    sample_set = set()
    for tsv in tsvfiles:
        _ , filename = os.path.split(tsv)
        sample = filename.split('.')[1]
        sample_set.add(sample)
        results.update({sample: {}})
        result_file = pd.read_csv(tsv, sep='\t', index_col=0).to_dict()
        values = result_file[list(result_file.keys())[0]]

        # Empty results file
        if values['summarized'] == '[]':
            continue

        for lineage, abundance in zip(values['lineages'].split(), values['abundances'].split()):
            lineage_set.add(lineage)
            results[sample].update({lineage: float(abundance)})

    lineages = sorted(list(lineage_set))
    samples = sorted(list(sample_set))

    df = pd.DataFrame(columns=lineages, index=samples)

    for sample, lineages in results.items():
        for lineage, abundance in lineages.items():
            df[lineage][sample] = abundance

    lab, run = outpath.split('/')
    outfile = '{}-{}'.format(lab, run)
    df.to_csv('{}/{}.freyja.csv'.format(outpath, outfile), encoding='utf-8')


def cutadapt(fq1, fq2, path="cutadapt", adapter="AGATCGGAAGAGC", ncores=1, minlen=10):
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
    _ = subprocess.check_call(cmd)
    of1.close()
    of2.close()
    return of1.name, of2.name


def minimap2(fq1, fq2, ref, path='minimap2', nthread=1):
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
            raise

    return tempbamsort.name 


def freyja(bamsort, ref, sample, outpath, path='freyja'):
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
        sys.stderr.write(f"Failed to find executable file at {path}")
        raise
    except subprocess.CalledProcessError:
        sys.stderr.write(f"Error running Freyja variants command:")
        sys.stderr.write(' '.join(cmd))
        raise

    bootstrap = [path, 'boot', '{}/var.{}.tsv'.format(outpath, sample),
                 '{}/depth.{}.csv'.format(outpath, sample),
                 '--nt', '4', '--nb', '10', '--output_base',
                 '{}/boostrap.{}'.format(outpath, sample), '--boxplot', 'pdf']
    try:
        _ = subprocess.check_call(bootstrap)
    except subprocess.CalledProcessError:
        sys.stderr.write(f"Error running Freyja boot command:")
        sys.stderr.write(' '.join(bootstrap))
        raise

    demix = [path, 'demix', '{}/var.{}.tsv'.format(outpath, sample),
             '{}/depth.{}.csv'.format(outpath, sample),
             '--output', '{}/lin.{}.tsv'.format(outpath, sample)]
    try:
        _ = subprocess.check_call(demix)
    except subprocess.CalledProcessError:
        sys.stderr.write(f"Error running Freyja demix command:")
        sys.stderr.write(' '.join(demix))
        raise


if __name__ == '__main__':
    # command line interface
    parser = argparse.ArgumentParser(
        description="Run Freyja Pipeline"
    )
    parser.add_argument('infile', type=str,
                        help="Path to the R1 FASTQ file of the sample")
    parser.add_argument('lab', type=str, default="freyja",
                        help="Name of the lab (western/guelph/waterloo)")
    parser.add_argument('run', type=str, default="freyja",
                        help="Name of the run that the sample belongs to")

    parser.add_argument('--ref', type=str, default="data/NC_045512.fa",
                        help="<input> path to target FASTA (reference)")
    parser.add_argument('--barcodes', type=str, default="data/usher_barcodes.csv",
                        help="<input> path to target FASTA (reference)")
    parser.add_argument('--cutadapt', type=str, default="cutadapt",
                        help="Path to cutadapt")
    parser.add_argument('--minimap2', type=str, default="minimap2",
                        help="Path to minimap2")
    parser.add_argument('--freyja', type=str, default="freyja",
                        help="Path to Freyja")
    parser.add_argument('--threads', type=int, default=4,
                        help="Number of threads (default 4)")

    args = parser.parse_args()

    cb = Callback()
    cb.callback("Starting Script")

    # prepare output directory
    outpath = '{}/{}'.format(args.lab, args.run)
    os.makedirs(outpath, exist_ok=True)

    fq1 = args.infile
    cb.callback("Running {}".format(fq1))
    fq2 = fq1.replace('_R1_', '_R2_')  # determine R2 filename
    _ , filename = os.path.split(fq1)
    prefix = filename.split('_')[0]
    
    tf1, tf2 = cutadapt(fq1=fq1, fq2=fq2, ncores=2, path=args.cutadapt)
    bamsort = minimap2(fq1=tf1, fq2=tf2, ref=args.ref, nthread=args.threads,
                       path=args.minimap2)
    frey = freyja(bamsort=bamsort, ref=args.ref, sample=prefix, outpath=outpath,
                  path=args.freyja)

    # Remove cutadapt output temporary files
    for tmpfile in [tf1, tf2, bamsort]:
        cb.callback("Removing {}".format(tmpfile))
        try:
            os.remove(tmpfile)
        except FileNotFoundError:
            pass
            
    # cb.callback("Generating CSV for {}".format('/'.join(outpath)))
    # generate_csv(outpath)
    
    cb.callback("All done")
