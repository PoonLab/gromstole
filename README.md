## Description

Gromstole is collection of Python and R scripts to estimate the relative frequencies of different SARS-CoV-2 [lineages](https://cov-lineages.org/), including specific "variants of concern", from next-generation sequence data (FASTQ) derived from wastewater samples.  This includes:
* a wrapper script (`minimap2.py`) to rapidly stream output from the reference mapping program [minimap](https://github.com/lh3/minimap2) to extract coverage and mutation frequency statistics for each sample;
* a script (`retrieve-nsgb.py`) for streaming an [open feed](https://nextstrain.org/blog/2021-07-08-ncov-open-announcement) of SARS-CoV-2 genome data and metadata, which is curated from the NCBI Genbank database by the [Nextstrain](https://nextstrain.org/) team and comprises over 3 million genomes (as of January 11, 2022);
* additional scripts for generating a set of lineage-specific mutations, by extracting the frequencies of mutations in association with specific lineages, and comparing those results to all other lineages as background;
* screening outputs from `minimap2.py` to estimate the frequencies of a specific lineage in a set of samples (`estimate-freqs.R`);
* generating visualizations, *e.g.*, `make-barplots.R`

## Dependencies
* [cutadapt](https://github.com/marcelm/cutadapt) 1.18+
* [minimap2](https://github.com/lh3/minimap2) version 2.17+
* [Python](https://www.python.org/) version 3.6+
* [R](https://cran.r-project.org/)

## Usage

Our current modelling strategy consists of a quasibinomial GLM for the count and coverage of the mutations that define a VOC, which are parsed from constellation JSON files curated by [cov-lineages](https://github.com/cov-lineages/constellations/).

The script `minimap2.py` is currently accepts paired-end reads in separate FASTQ files and outputs `[prefix]-mapped.csv` and `[prefix]-coverage.csv` into the specified output directory. It uses `data/NC043312.fa` as a reference genome, but an alternative FASTA file can be specified by the user.  By default, it uses cutadapt to trim adapter sequences from the data.

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

To get a nice summary of the results, the following scripts produces a stacked bar plot with the specified filename.
```sh
Rscript scripts/make-barplots.R results/outfile_B-1-1-529.json results/barplot_B-1-1-529.pdf
```
