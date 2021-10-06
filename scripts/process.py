import argparse
import os
import sys
import json
from datetime import datetime
import lzma
import subprocess

import seq_utils
from progress_utils import Callback

from minimap2 import extract_features


def parse_args():
    parser = argparse.ArgumentParser(
        description="Pipeline automation for execution on local files"
    )

    parser.add_argument("--outdir", type=str, default='data/',
                        help="option, path to write output files")

    parser.add_argument('--poisson-cutoff', type=float, default=0.001,
                        help='filtering outlying genomes whose distance exceeds the upper '
                             'quantile of Poisson distribution (molecular clock).  Default 0.001 '
                             'corresponds to 99.9%% cutoff.')
    parser.add_argument('--minlen', type=int, default=29000, help='minimum genome length (nt)')

    parser.add_argument('--batchsize', type=int, default=2000,
                        help='number of records to batch process with minimap2')

    parser.add_argument("--ref", type=str,
                        default="data/NC_045512.fa",
                        help="path to FASTA file with reference genome")
    parser.add_argument('--mmbin', type=str, default='minimap2',
                        help="path to minimap2 binary executable")
    parser.add_argument('-mmt', "--mmthreads", type=int, default=16,
                        help="number of threads for minimap2.")

    parser.add_argument('--misstol', type=int, default=300,
                        help="maximum tolerated number of missing bases per "
                             "genome (default 300).")
    parser.add_argument("--vcf", type=str,
                        default=os.path.join(covizu.__path__[0], "data/problematic_sites_sarsCov2.vcf"),
                        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
                             "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")

    parser.add_argument('--clock', type=float, default=8e-4,
                        help='specify molecular clock rate for '
                             'constraining Treetime analysis (default 8e-4).')

    parser.add_argument('--datetol', type=float, default=0.1,
                        help='exclude tips from time-scaled tree '
                             'with high discordance between estimated and '
                             'known sample collection dates (year units,'
                             'default: 0.1)')

    return parser.parse_args()



def sort_by_lineage(records, callback=None):
    """
    Resolve stream into a dictionary keyed by Pangolin lineage
    :param records:  generator, return value of extract_features()
    :return:  dict, lists of records keyed by lineage
    """
    result = {}
    for i, record in enumerate(records):
        if callback and i % 1000 == 0:
            callback('aligned {} records'.format(i))
        lineage = record['lineage']
        if lineage is None or lineage == '':
            continue  # discard unclassified genomes

        if lineage not in result:
            result.update({lineage: []})
        result[lineage].append(record)

    return result


def analyze_feed(handle, args, callback=None):
    """
    :param handle:  file stream in read mode, from lzma.open()
    :param args:  Namespace
    :param callback:  optional progress monitoring, see progress_utils.py
    """
    # check that user has loaded openmpi module
    if args.machine_file or args.np:
        try:
            subprocess.check_call(['mpirun', '-np', '2', 'ls'], stdout=subprocess.DEVNULL)
        except FileNotFoundError:
            if callback:
                callback("mpirun not loaded - run `module load openmpi/gnu`", level='ERROR')
            sys.exit()

    # pre-processing feed
    feed = map(json.loads, handle)
    batcher = seq_utils.batch_fasta(feed, size=args.batchsize)
    aligned = extract_features(batcher, ref_file=args.ref, binpath=args.mmbin,
                               nthread=args.mmthreads, minlen=args.minlen)
    filtered = seq_utils.filter_problematic(aligned, vcf_file=args.vcf, cutoff=args.poisson_cutoff,
                                            callback=callback)
    by_lineage = sort_by_lineage(filtered, callback=callback)


def get_mutations(by_lineage):
    """
    Extract common mutations from feature vectors for each lineage
    :param by_lineage:  dict, return value from process_feed()
    :return:  dict, common mutations by lineage
    """
    result = {}
    for lineage, samples in by_lineage.items():
        # enumerate features
        counts = {}
        for sample in samples:
            for diff in sample['diffs']:
                feat = tuple(diff)
                if feat not in counts:
                    counts.update({feat: 0})
                counts[feat] += 1
        # filter for mutations that occur in at least half of samples
        common = [feat for feat, count in counts.items() if count/len(samples) >= 0.5]
        result.update({lineage: common})
    return result


if __name__ == '__main__':
    cb = Callback()
    args = parse_args()
    handle = lzma.open(args.infile, 'rb')
    analyze_feed(handle, args, callback=cb.callback)
