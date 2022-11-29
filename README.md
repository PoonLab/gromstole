Gromstole is collection of Python and R scripts to estimate the relative frequencies of different SARS-CoV-2 [lineages](https://cov-lineages.org/), including specific "variants of concern", from next-generation sequence data (FASTQ) derived from wastewater samples.  This includes:
* a wrapper script (`minimap2.py`) to rapidly stream output from the reference mapping program [minimap](https://github.com/lh3/minimap2) to extract coverage and mutation frequency statistics for each sample;
* an R script (`estimate-freqs.R`) that uses quasibinomial regression to estimate the frequency of a variant of concern from the frequencies of mutations based on the associated [constellation file](https://github.com/cov-lineages/constellations/)
* a second R script (`make-barplots.R`) for visualizing these variant frequency estimates across a set of samples as a barplot


## Dependencies
* [cutadapt](https://github.com/marcelm/cutadapt) 1.18+
* [minimap2](https://github.com/lh3/minimap2) version 2.17+
* [Python](https://www.python.org/) version 3.6+
* [R](https://cran.r-project.org/)

## Usage

The following summarizes a workflow based on a pair of FASTQ files associated with a [study of wastewater by the Nevada State Public Health Laboratory](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP354147).  These data can be obtained from the NCBI Sequence Read Archive using the command line tool `fasterq-dump` that is distributed with the [sra-tools](https://github.com/ncbi/sra-tools) package.

The `consellations` submodule needs to be initialized prior to running the scripts for the first time. 

```console
git submodule init; git submodule update
```

Prior to running the scripts, the submodule should be updated to use the latest constellation JSON files.

```console
git submodule foreach git pull origin main
```

```console
art@Langley:~/git/gromstole$ ~/src/sratoolkit.2.11.3-ubuntu64/bin/fasterq-dump -p SRR17724299
join   :|-------------------------------------------------- 100%   
concat :|-------------------------------------------------- 100%   
spots read      : 614,374
reads read      : 1,228,748
reads written   : 1,228,748
art@Langley:~/git/gromstole$ python3 scripts/minimap2.py SRR17724299_1.fastq SRR17724299_2.fastq -o testrun
Defaulting output prefix stem to SRR17724299_1.fastq
50000.0 reads, 46769.0 (94%) mapped
100000.0 reads, 93735.0 (94%) mapped
150000.0 reads, 140548.0 (94%) mapped
200000.0 reads, 187560.0 (94%) mapped
250000.0 reads, 234509.5 (94%) mapped
300000.0 reads, 281501.0 (94%) mapped
350000.0 reads, 328490.0 (94%) mapped
400000.0 reads, 375566.0 (94%) mapped
450000.0 reads, 422431.5 (94%) mapped
500000.0 reads, 469568.5 (94%) mapped
550000.0 reads, 516504.5 (94%) mapped
600000.0 reads, 563587.5 (94%) mapped
art@Langley:~/git/gromstole$ Rscript scripts/estimate-freqs.R testrun constellations/constellations/definitions/cB.1.617.2.json \
testrun/SRR17724299_1.delta.json 
Loading required package: lubridate

Attaching package: ‘lubridate’

The following objects are masked from ‘package:base’:

    date, intersect, setdiff, union

Loading required package: jsonlite
Loading required package: seqinr
Loading required package: here
here() starts at /home/art/git/gromstole
```

This yields the following JSON file:
```json
{
  "results": [
    {
      "sample": "SRR17724299_1",
      "mutation": "S:T19R",
      "nucleotide": "C21618G",
      "count": 643,
      "coverage": 653
    },
    {
      "sample": "SRR17724299_1",
      "mutation": "S:P681R",
      "nucleotide": "C23604G",
      "count": 2010,
      "coverage": 2025
    },
    {
      "sample": "SRR17724299_1",
      "mutation": "S:D950N",
      "nucleotide": "G24410A",
      "count": 1629,
      "coverage": 2078
    },
    {
      "sample": "SRR17724299_1",
      "mutation": "ORF3a:S26L",
      "nucleotide": "C25469T",
      "count": 1539,
      "coverage": 1556
    },
    {
      "sample": "SRR17724299_1",
      "mutation": "M:I82T",
      "nucleotide": "T26767C",
      "count": 1915,
      "coverage": 1922
    },
    {
      "sample": "SRR17724299_1",
      "mutation": "ORF7a:V82A",
      "nucleotide": "T27638C",
      "count": 501,
      "coverage": 502
    },
    {
      "sample": "SRR17724299_1",
      "mutation": "ORF7a:T120I",
      "nucleotide": "C27752T",
      "count": 563,
      "coverage": 567
    },
    {
      "sample": "SRR17724299_1",
      "mutation": "N:D63G",
      "nucleotide": "A28461G",
      "count": 2056,
      "coverage": 2078
    },
    {
      "sample": "SRR17724299_1",
      "mutation": "N:R203M",
      "nucleotide": "G28881T",
      "count": 845,
      "coverage": 852
    },
    {
      "sample": "SRR17724299_1",
      "mutation": "N:D377Y",
      "nucleotide": "G29402T",
      "count": 646,
      "coverage": 655
    }
  ],
  "metadata": [
    {
      "sample": "SRR17724299_1",
      "lab": "gromstole"
    }
  ],
  "estimate": [
    {
      "est": 0.9576,
      "lower.95": 0.8888,
      "upper.95": 0.9903,
      "_row": "SRR17724299_1"
    }
  ],
  "lineage": ["B.1.617.2"],
  "run.dir": ["testrun"]
}
```

## Description

Our current modelling strategy consists of a quasibinomial GLM for the count and coverage of the mutations that define a VOC, which are parsed from constellation JSON files curated by [cov-lineages](https://github.com/cov-lineages/constellations/).

The script `minimap2.py` currently accepts paired-end reads in separate FASTQ files and outputs `[prefix]-mapped.csv` and `[prefix]-coverage.csv` into the specified output directory. It uses `data/NC043312.fa` as a reference genome, but an alternative FASTA file can be specified by the user.  By default, it uses cutadapt to trim adapter sequences from the data.

Use `python3 scripts/minimap2.py -h` to see all of the options.

```sh
python3 scripts/minimap2.py r1.fastq r2.fastq --outdir results --prefix name-of-sample
```

Given a directory containing one or more subdirectories that each contain paired-end Illumina FASTQ fiiles, the following R script runs the binomial regression and outputs a json file that contains all relevant information:

- the counts of each mutation of the `lineages/` file
- the coverage at every position on the reference genome, 
- the metadata that was used as input (optional), 
- the estimate of the proportion (including 95% confidence interval), 
- the lineage name, and 
- the name of the input directory.

```sh
Rscript scripts/estimate-freqs.R results/name-of-sample constellations/constellations/definitions/cBA.1.csv results/outfile.json path/to/metadata.csv
```

