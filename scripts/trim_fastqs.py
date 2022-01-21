#!/usr/bin/env python3.6

"""
Censor reads based on phiX quality data, and also trim adapter sequences.
adapted from https://github.com/cfe-lab/MiCall/blob/master/micall/core/trim_fastqs.py
as distributed under the AGPLv3 license
"""

import argparse
import csv
import shutil
from gzip import GzipFile
import itertools
import math
import os
from tempfile import NamedTemporaryFile
import subprocess

from io import TextIOWrapper
from itertools import chain
from pathlib import Path

from Bio import SeqIO


def trim(fq1, fq2, out1, out2, badfile=None, use_gzip=True, datadir=os.getcwd()):
    """
    :param fq1:  str, path to first FASTQ file
    :param fq2:  str, path to second FASTQ file
    :param datadir:  str, path to directory containing adapter and primer FASTAs
    :param out1:  str, path to write first trimmed FASTQ file
    :param out2:  str, path to write second trimmed FASTQ file
    :param badfile:  str, path to CSV file with 'tile' and 'cycle' columns
    :param use_gzip: True if FASTQ files are gzip-compressed
    """
    if badfile:
        with open(badfile, 'r') as bad_cycles:
            bad_cycles = list(csv.DictReader(bad_cycles))

        src1 = open(fq1, 'rb')
        dest1 = NamedTemporaryFile(delete=False)
        censor(src1, bad_cycles, dest1.name, use_gzip, cycle_sign=1)
        src2 = open(fq2, 'rb')
        dest2 = NamedTemporaryFile(delete=False)

        censor(src2, bad_cycles, dest2.name, use_gzip, cycle_sign=-1)
        cut_all(dest1, dest2, out1, out2, datadir=datadir)
        purge_temp_files([dest1, dest2])
    else:
        cut_all(fq1, fq2, out1, out2, datadir=datadir)


def cut_all(fq1, fq2, out1, out2, datadir):
    dedapted = [NamedTemporaryFile(delete=False), NamedTemporaryFile(delete=False)]
    cut_adapters(fq1, fq2, dedapted[0].name, dedapted[1].name, datadir)

    ltrimmed = [NamedTemporaryFile(delete=False), NamedTemporaryFile(delete=False)]
    cut_lr_primers(dedapted[0].name, dedapted[1].name, ltrimmed[0].name, ltrimmed[1].name,
                   datadir=datadir, left=True)

    rtrimmed = [NamedTemporaryFile(delete=False), NamedTemporaryFile(delete=False)]
    cut_lr_primers(dedapted[0].name, dedapted[1].name, rtrimmed[0].name, rtrimmed[1].name,
                   datadir=datadir, left=False)

    combine_primer_trimming(
        dedapted[0].name, dedapted[1].name, ltrimmed[0].name, ltrimmed[1].name,
        rtrimmed[0].name, rtrimmed[1].name, out1, out2
    )
    purge_temp_files(dedapted)
    purge_temp_files(ltrimmed)
    purge_temp_files(rtrimmed)


def purge_temp_files(filenames):
    for filename in filenames:
        try:
            os.remove(filename)
        except OSError:
            # We tried to tidy up a temporary file, but it's not critical.
            pass


def cut_adapters(original_fastq1, original_fastq2, trimmed_fastq1, trimmed_fastq2, datadir):
    adapter_files = [os.path.join(datadir, f'adapters_read{i}.fasta') for i in (1, 2)]
    subprocess.check_call([
        'cutadapt', '--quiet',
        '-a', 'file:' + adapter_files[0],
        '-A', 'file:' + adapter_files[1],
        '-o', str(trimmed_fastq1),
        '-p', str(trimmed_fastq2),
        str(original_fastq1),
        str(original_fastq2)
    ])


def cut_lr_primers(original_fastq1, original_fastq2, trimmed_fastq1, trimmed_fastq2, datadir, left=True):
    """ Trim left/right primers off both sets of reads.

    The filtering options are a little tricky. --trimmed-only means that only
    the reads that got filtered are written into the ltrimmed.fastq file.
    --pair-filter=both means that both reads have to be untrimmed before they
    will get excluded from the ltrimmed.fastq file. In other words, if either
    read had a primer trimmed off, then both reads will get written to the
    ltrimmed.fastq file.
    """
    lr1 = 'left' if left else 'right_end'
    primer_file1 = os.path.join(datadir, f'file:primers_sarscov2_{lr1}.fasta')
    lr2 = 'left_end' if left else 'right'
    primer_file2 = os.path.join(datadir, f'file:primers_sarscov2_{lr2}.fasta')
    subprocess.check_call(['cutadapt', '--quiet',
                           '-g' if left else '-a', primer_file1,
                           '-A' if left else '-G', primer_file2,
                           '-o', str(trimmed_fastq1),
                           '-p', str(trimmed_fastq2),
                           '--overlap=8', '--trimmed-only', '--pair-filter=both',
                           str(original_fastq1),
                           str(original_fastq2)])


def combine_primer_trimming(original_fastq1, original_fastq2, ltrimmed_fastq1, ltrimmed_fastq2,
                            rtrimmed_fastq1, rtrimmed_fastq2, trimmed_fastq1, trimmed_fastq2):
    trimmed_sequences1 = cut_primer_dimer_sequences(
        original_fastq1,
        start_trimmed_fastq=ltrimmed_fastq1,
        end_trimmed_fastq=rtrimmed_fastq1)
    trimmed_sequences2 = cut_primer_dimer_sequences(
        original_fastq2,
        start_trimmed_fastq=rtrimmed_fastq2,
        end_trimmed_fastq=ltrimmed_fastq2)

    with trimmed_fastq1.open('w') as f1, trimmed_fastq2.open('w') as f2:
        for seq1, seq2 in zip(trimmed_sequences1, trimmed_sequences2):
            if len(seq1) > 0 and len(seq2) > 0:
                SeqIO.write([seq1], f1, 'fastq')
                SeqIO.write([seq2], f2, 'fastq')


def cut_primer_dimer_sequences(original_fastq,
                               start_trimmed_fastq,
                               end_trimmed_fastq):
    start_trimmed_seqs = chain(SeqIO.parse(start_trimmed_fastq, 'fastq'), [None])
    end_trimmed_seqs = chain(SeqIO.parse(end_trimmed_fastq, 'fastq'), [None])
    next_start_trimmed_seq = next(start_trimmed_seqs)
    next_end_trimmed_seq = next(end_trimmed_seqs)
    for original_seq in SeqIO.parse(original_fastq, 'fastq'):
        if (next_start_trimmed_seq is None or
                next_start_trimmed_seq.id != original_seq.id):
            start = 0
        else:
            start = len(original_seq) - len(next_start_trimmed_seq)
            next_start_trimmed_seq = next(start_trimmed_seqs)

        if (next_end_trimmed_seq is None or
                next_end_trimmed_seq.id != original_seq.id):
            end = None
        else:
            end = len(next_end_trimmed_seq)
            next_end_trimmed_seq = next(end_trimmed_seqs)
        yield original_seq[start:end]


def censor(original_file, bad_cycles_reader, censored_file, use_gzip=True,
           summary_writer=None, cycle_sign=1):
    """ Censor bases from a FASTQ file that were read in bad cycles.

    @param original_file: an open FASTQ file to read from
    @param bad_cycles_reader: an iterable collection of bad cycle entries:
        {'tile': tile, 'cycle': cycle}
    @param censored_file: an open FASTQ file to write to: censored bases will
        be written as 'N' with a quality '#'.
    @param use_gzip: True if the original file should be unzipped
    @param summary_writer: an open CSV DictWriter to write to: write a single row
        with the average read quality for the whole sample
    @param cycle_sign: +1 or -1, shows which direction to run cycles when
        checking bad_cycles.
    """
    bad_cycles = set()
    for cycle in bad_cycles_reader:
        bad_cycles.add((cycle['tile'], int(cycle['cycle'])))

    src = original_file
    dest = open(censored_file, 'w')
    base_count = 0
    score_sum = 0.0
    if use_gzip:
        src = GzipFile(fileobj=original_file)
    src = TextIOWrapper(src)

    for ident, seq, opt, qual in itertools.zip_longest(src, src, src, src):
        # returns an aggregate of 4 lines per call
        dest.write(ident)
        if not bad_cycles:
            dest.write(seq)
            tile = None
        else:
            ident_fields, read_fields = map(str.split, ident.split(' '), '::')
            tile = ident_fields[4]
            bad_count = 0
            for cycle, base in enumerate(seq.rstrip(), start=1):
                cycle = math.copysign(cycle, cycle_sign)
                if (tile, cycle) in bad_cycles:
                    bad_count += 1
                else:
                    if bad_count:
                        dest.write('N' * bad_count)
                        bad_count = 0
                    dest.write(base)
            dest.write('\n')
        dest.write(opt)

        bad_count = 0
        for cycle, score in enumerate(qual.rstrip(), start=1):
            float_score = ord(score) - 33
            score_sum += float_score
            base_count += 1
            cycle = math.copysign(cycle, cycle_sign)
            if (tile, cycle) in bad_cycles:
                bad_count += 1
            else:
                if bad_count:
                    dest.write('#' * bad_count)
                    bad_count = 0
                dest.write(score)
        dest.write('\n')

    if summary_writer is not None:
        avg_quality = score_sum/base_count if base_count > 0 else None
        summary = dict(base_count=base_count,
                       avg_quality=avg_quality)
        summary_writer.writerow(summary)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Trim adapter and primer sequences, and censor bad tiles '
                    'and cycles from paired FASTQ files.')

    parser.add_argument('fq1', help='<input> FASTQ.gz containing original reads (read 1)')
    parser.add_argument('fq2', help='<input> FASTQ.gz containing original reads (read 2)')
    parser.add_argument('out1', help='<output> uncompressed FASTQ containing trimmed reads (read 1)')
    parser.add_argument('out2', help='<output> uncompressed FASTQ containing trimmed reads (read 2)')
    parser.add_argument('--datadir', default=os.getcwd(),
                        help='<input> path to directory containing adapter and primer FASTA files')

    parser.add_argument('--bad_cycles_csv', default=None,
                        help='<input> List of tiles and cycles rejected for poor quality')
    parser.add_argument('--unzipped', '-u', action='store_true',
                        help='Set if the original FASTQ files are not compressed')

    args = parser.parse_args()

    trim(fq1=args.fq1, fq2=args.fq2, out1=args.out1, out2=args.out2,
         badfile=args.bad_cycles_csv, datadir=args.datadir,
         use_gzip=not args.unzipped)
