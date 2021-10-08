from datetime import date
import bisect
import pkg_resources


complement_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
                    'W':'S', 'R':'Y', 'K':'M', 'Y':'R', 'S':'W', 'M':'K',
                    'B':'V', 'D':'H', 'H':'D', 'V':'B',
                    '*':'*', 'N':'N', '-':'-'}


def revcomp(seq):
    rseq = seq[::-1]
    rcseq = ''
    for i in rseq:  # reverse order
        rcseq += complement_dict[i]
    return rcseq

complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'W': 'S', 'R': 'Y', 'K': 'M', 'Y': 'R', 'S': 'W', 'M': 'K',
                   'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
                   '*': '*', 'N': 'N', '-': '-'}


def revcomp(seq):
    """ Reverse-complement input nucleotide sequence """
    rseq = seq[::-1]
    rcseq = ''
    for i in rseq:  # reverse order
        rcseq += complement_dict[i]
    return rcseq


def read_seq(handle):
    """
    Read sequence from plain text file (no format).  Used for importing
    reference sequence.

    :param handle:
    :return:  str, sequence
    """
    seq = ''
    for line in handle:
        seq += line.strip()
    return seq


def parse_label(label):
    """
    Extract country and date of sample collection from GISAID label.
    Return date as None if collection date is ambiguous.

    :param label:  str, GISAID sequence label
    :return: (country, date)
    """
    info, epi_id, ymd = label.split('|')
    country = info.split('/')[1]
    try:
        year, month, day = list(map(int, ymd.split('-')))
        return country, date(year, month, day)
    except ValueError:
        return country, None
    except:
        raise


def iter_fasta(handle):
    """
    Parse open file as FASTA.  Returns a generator
    of handle, sequence tuples.

    :param handle:  open stream to FASTA file in read mode
    :yield tuples, (header, sequence)
    """
    h, sequence = None, ''
    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                yield h, sequence
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()
    yield h, sequence


def convert_fasta(handle):
    """
    Parse FASTA file as a list of header, sequence list objects
    :param handle:  open file stream
    :return:  List of [header, sequence] records
    """
    result = []
    h, sequence = None, ''
    for line in handle:
        if line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h, sequence])
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()
    result.append([h, sequence])  # handle last entry
    return result


def total_missing(row):
    """ Calculate the total number of missing sites from closed-open interval annotations """
    res = 0
    if type(row) is dict:
        missing = row['missing']
    else:
        _, _, missing = row

    for left, right in missing:
        res += right-left
    return res


def apply_features(diffs, missing, refseq):
    """
    Reconstitute genome sequence from feature vector (genetic differences) and
    missing data vector.

    :param row:  list, entry from features list returned by import_json()
    :param refseq:  str, reference genome
    :return:  str, aligned genome
    """
    result = list(refseq)  # strings are not mutable

    # apply missing intervals
    for left, right in missing:
        for i in range(left, right):
            result[i] = 'N'

    # apply substitutions and deletions (skip insertions)
    for dtype, pos, diff in diffs:
        if dtype == '~':
            result[pos] = diff
        elif dtype == '-':
            for i in range(pos, pos+diff):
                result[i] = '-'

    return ''.join(result)


def fromisoformat(dt):
    """ Convert ISO date to Python datetime.date object to support Python <3.7 """
    year, month, day = map(int, dt.split('-'))
    return date(year, month, day)


def load_vcf(vcf_file="data/problematic_sites_sarsCov2.vcf"):
    """
    Load VCF of problematic sites curated by Nick Goldman lab
    NOTE: The curators of this VCF used MN908947.3, which is identical to NC_045512.
    *** It is very important to check that your reference is compatible! ***
    TODO: align user's reference to NC_045512 to generate custom coordinate system

    :param vcf_file:  str, path to VCF file
    :return:  dict, tuples keyed by reference coordinate
    """
    vcf = open(vcf_file)
    mask = {}
    for line in vcf.readlines():
        if line.startswith('#'):
            continue
        try:
            _, pos, _, ref, alt, _, filt, info = line.strip().split('\t')
        except ValueError:
            raise
        if filt == 'mask':
            mask.update({int(pos)-1: {  # convert to 0-index
                'ref': ref, 'alt': alt, 'info': info}
            })
    return mask


def filter_problematic_sites(obj, mask, callback=None):
    """
    Apply problematic sites annotation from de Maio et al.,
    https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
    which are published and maintained as a VCF-formatted file.
    FIXME: this duplicates some functionality of filter_problematic(), #290

    :param obj:  list, entries are (1) dicts returned by import_json or (2) tuples
    :param mask:  dict, problematic site index from load_vcf()
    :param vcf_file:  str, path to VCF file
    :return:
    """
    # apply filters to feature vectors
    count = 0
    result = []
    for row in obj:
        if type(row) is dict:
            qname, diffs, missing = row['qname'], row['diffs'], row['missing']
        else:
            qname, diffs, missing = row  # unpack tuple

        filtered = []
        for typ, pos, alt in diffs:
            if typ == '~' and int(pos) in mask and alt in mask[pos]['alt']:
                continue
            if typ != '-' and 'N' in alt:
                # drop substitutions and insertions with uncalled bases
                continue
            filtered.append(tuple([typ, pos, alt]))

        count += len(diffs) - len(filtered)
        result.append([qname, filtered, missing])

    if callback:
        callback('filtered {} problematic features'.format(count))
    return result


def batch_fasta(gen, size=100):
    """
    Concatenate sequence records in stream into FASTA-formatted text in batches of
    <size> records.
    :param gen:  generator, return value of load_gisaid()
    :param size:  int, number of records per batch
    :yield:  str, list; FASTA-format string and list of records (dict) in batch
    """
    stdin = ''
    batch = []
    for i, record in enumerate(gen, 1):
        qname, sequence = record
        stdin += '>{}\n{}\n'.format(qname, sequence)
        if i > 0 and i % size == 0:
            yield stdin
            stdin = ''

    if stdin:
        yield stdin


def filter_problematic(records, origin='2019-12-01', rate=0.0655, cutoff=0.005,
                       maxtime=1e3, vcf_file='data/problematic_sites_sarsCov2.vcf',
                       misstol=300, callback=None):
    """
    Apply problematic sites annotation from de Maio et al.,
    https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
    which are published and maintained as a VCF-formatted file.

    :param records:  generator, records from extract_features()
    :param origin:  str, date of root sequence in ISO format (yyyy-mm-dd)
    :param rate:  float, molecular clock rate (subs/genome/day), defaults
                  to 8e-4 * 29900 / 365
    :param cutoff:  float, use 1-cutoff to compute quantile of Poisson
                    distribution, defaults to 0.005
    :param maxtime:  int, maximum number of days to cache Poisson quantiles
    :param vcf_file:  str, path to VCF file
    :param misstol:  int, maximum tolerated number of uncalled bases
    :param callback:  function, option to print messages to console
    :yield:  generator, revised records
    """
    # load resources
    mask = load_vcf(vcf_file)
    qp = QPois(quantile=1-cutoff, rate=rate, maxtime=maxtime, origin=origin)

    n_sites = 0
    n_outlier = 0
    n_ambig = 0
    for record in records:
        # exclude problematic sites
        filtered = []
        diffs = record['diffs']
        for typ, pos, alt in diffs:
            if typ == '~' and int(pos) in mask and alt in mask[pos]['alt']:
                continue
            if typ != '-' and 'N' in alt:
                # drop substitutions and insertions with uncalled bases
                continue
            filtered.append(tuple([typ, pos, alt]))

        ndiffs = len(filtered)
        n_sites += len(diffs) - ndiffs
        record['diffs'] = filtered

        # exclude genomes with excessive divergence from reference
        coldate = record['coldate']
        if qp.is_outlier(coldate, ndiffs):
            n_outlier += 1
            continue

        # exclude genomes with too much missing data
        if total_missing(record) > misstol:
            n_ambig += 1
            continue

        yield record

    if callback:
        callback("filtered {} problematic features".format(n_sites))
        callback("         {} genomes with excess missing sites".format(n_ambig))
        callback("         {} genomes with excess divergence".format(n_outlier))
