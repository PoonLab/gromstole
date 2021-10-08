# Wastewater Surveillance of SARS-CoV-2

Please expand this list!

- Chrystal Landgraff's paper: 
    - Landgraff, C., Wang, L. Y. R., Buchanan, C., Wells, M., Schonfeld, J., Bessonov, K., Ali, J., Robert, E., & Nadon, C. (2021). Metagenomic sequencing of municipal wastewater provides a near-complete SARS-CoV-2 genome sequence identified as the B.1.1.7 variant of concern from a Canadian municipality concurrent with an outbreak [Preprint]. Public and Global Health. https://doi.org/10.1101/2021.03.11.21253409

- Emmanuel's journal club paper:
    - Huisman, J. S., Scire, J., Caduff, L., Fernandez, X., Ganesanandamoorthy, P., Kull, A., Scheidegger, A., Boehm, A. B., Hughes, B., Knudson, A., Topol, A., Wigginton, K. R., Wolfe, M. K., Kohn, T., Ort, C., & Julian, T. R. (n.d.). Wastewater-based estimation of the effective reproductive number of SARS-CoV-2. 39.

- Rios G, Lacoux C, Leclercq V, Diamant A, Lebrigand K, Lazuka A, Soyeux E, Lacroix S, Fassy J, Couesnon A, Thiery R. Monitoring SARS-CoV-2 variants alterations in Nice neighborhoods by wastewater nanopore sequencing. The Lancet Regional Health-Europe. 2021 Aug 17:100202.
  * "Available bioinformatics workflows for SARS-CoV-2 sequencing data were not directly suitable for wastewater data since they were designed for the analysis of individual patient data, where each sample is assigned to one single lineage."
  * adapted pipeline - iVar to generate mutation frequency data from BAM outputs
  * filtered for positions with >100 read depth
  * used R to fix "known issue" with iVar output (incorrect assignment of deletion frequencies)
  * compiled database of mutations associated with known lineages from https://covidcg.org
  * 

- Wurtz N, Revol O, Jardot P, Giraud-Gatineau A, Houhamdi L, Soumagnac C, Annessi A, Lacoste A, Colson P, Aherfi S, Scola BL. Monitoring the Circulation of SARS-CoV-2 Variants by Genomic Analysis of Wastewater in Marseille, South-East France. Pathogens. 2021 Aug;10(8):1042.
  * NGS libraries prepared using Illumina COVIDSeq protocol
  * single-end sequencing with 36bp read length on NovaSeq 6000
  * reads mapped to reference genome using CLC genomics software v7.5
  * compared nonsynonymous mutations to "classifying mutations that match with 30 SARS-CoV-2 variants circulating in France

  > Fairly crummy bioinformatics, nothing much of note here

- Crits-Christoph A, Kantor RS, Olm MR, Whitney ON, Al-Shayeb B, Lou YC, Flamholz A, Kennedy LC, Greenwald H, Hinkle A, Hetzel J. Genome sequencing of sewage detects regionally prevalent SARS-CoV-2 variants. MBio. 2021 Jan 19;12(1):e02703-20.
  * "Sequencing viral concentrates and RNA extracted from wastewater can identify multiple SARS-CoV-2 genotypes at various abundances konwn to be present in communities, as well as additional genotypic variants not yet observed in local clinical sequencing efforts."
  * Illumina NextSeq (paired 2x75bp reads), total RNA
  * cDNA enrichment with Illumina Respiratory Virus Oligo Panel
  * rRNA depletion with Gut Microbiome probe set
  * mapped with Bowtie2 to all virus genomes in RefSeq
  * removed duplicate reads, only reported genomes with >10% coverage
  * SNV calling with inStrain for read pairs with >90% identity to reference
  * removed PCR duplicates with Sambamba
  * hypergeometric distribution testing
  
  > SARS-CoV-2 genomics done by metagenomics people

- [Guidance for surveillance of SARS-CoV-2 variants: Interim guidance, 9 August 2021](https://www.who.int/publications/i/item/WHO_2019-nCoV_surveillance_variants)

- Wastewater Monitoring of SARS-CoV-2 Variants in England: Demonstration Case Study for Bristol (Dec 2020 - March 2021)
[Summary for SAGE 08/04/21](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/979864/S1193_Wastewater_Monitoring_of_SARS-CoV-2_Variants_in_England_Demonstration_Case_Study_for_Bristol__Dec_2020-March_2021_.pdf)
  * used "multi-mutational approach [...], aimed at identifying 118 discrete mutations [...] which are signature mutations for known Variants of Concern (VoC) and Variants Under Investigation (VUI)"
  * ARTIC v3 protocol (tiled 400bp amplicons)
  * sequenced on Illumina MiSeq (2x250bp reads)
  * raw reads processed with ncov2019-artic-nf v3 pipeline (Tyson et al. 2020)
    * default parameters (align reads to ref genome with minimap v2.17
    * identified SNPs and indels from BAM files with samtools (v0.1.18-r580) and VarScan (v2.3) on 100,000 reads with alignment score >10
    * required mutations to be present in both reads of pair
  * mutations compared against signature mutations defined by Public Health England (https://github.com/phe-genomics/variant_definitions)
  * genome-wide co-occurrence approach for VOC/VUI identification - using CoOccurrence adJusted Analysis and Calling (COJAC; Jahn et al. 2021)
  * "Since a proportion of signature mutations are shared amongst these VOCs/VUIs, amplicons with co-occurring mutations can be listed for multiple variants."

- Jahn K, Dreifuss D, Topolsky I, Kull A, Ganesanandamoorthy P, Fernandez-Cassi X, Bänziger C, Devaux AJ, Stachler E, Caduff L, Cariti F. Detection and surveillance of SARS-CoV-2 genomic variants in wastewater. [medRxiv. 2021 Jan 9:2021-01](https://www.medrxiv.org/content/10.1101/2021.01.08.21249379v2).
  * ARTIC v3 protocol, fragments enriched and barcoded with unique dual indexing
  * libraries sequenced on Illumina NovaSeq 6000 and MiSeq (paired-end 2x250bp)
  * data analyzed using V-pipe
  * individual low-frequency mutations called based on local haplotype reconstruction using ShoRAH
  * developed new tool (called [COJAC](https://github.com/cbg-ethz/cojac/)) for detecting mutational co-occurrence - takes multiple read alignments (BAM files) and counts read pairs with variant-specific mutational patterns
  * "Detecting multiple signature mutations on the same amplicon increases the confidence of mutation calls at very low variant read counts."
  * only used non-synonymous substitutions for quantification
    * frequencies of B.1.1.7 signature substitutions were resampled with replacement, averaged per wastewater sample, smoothed over time by local regression (lowess, Python statsmodels v0.12.1) to construct bootstrap estimates of B.1.1.7 per-day frequency curves
  
  > developed COJAC (co-occurrence method), worth investigating - can probably ignore local haplotype reconstruction step (one of the authors is a big advocate of this approach)

# Alignment with uncertainty

- GNUMAP: Find every possible way to align the read to the reference, then find the quality of each alignmnent.
    - Clement, N. L., Snell, Q., Clement, M. J., Hollenhorst, P. C., Purwar, J., Graves, B. J., Cairns, B. R., & Johnson, W. E. (2010). The GNUMAP algorithm: Unbiased probabilistic mapping of oligonucleotides from next-generation sequencing. Bioinformatics, 26(1), 38–45. https://doi.org/10.1093/bioinformatics/btp614
    
