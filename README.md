# Wastewater Surveillance for SARS-CoV-2 Variants of Concern

A collaboration between:

- The University of Western Ontario (UWO)
    - Dr. Art Poon's bioinformatics lab
    - Dr. Eric Art's molecular biology lab
        - There are too many "Art"s in this project
- The University of Waterloo (UW)
    - Dr. Trevor Charles' molecular biology lab
    - There are too many U's, W's, and O's in this project
- Public Health Agency of Canada's (PHAC) National Microbiology Laboratory (NML)
- Ontario's Ministry of the Environment, Conservation and Parks (MECP)

## Wastewater collection

Wastewater will be collected from Pearson International Airport, the three universities involved, and various OnRoute locations. Some of the sampling will be automated, some will be manual sample collection.

The sequencing methods are based on the methods of [Chrystal Landgraff's preprint](http://medrxiv.org/lookup/doi/10.1101/2021.03.11.21253409)
Basically, this is an implementation of the Artic pipeline that has been modified for wastewater.


## Data Format

Currently, we expect the data to come in FASTQ format containing many, many short reads of SARS-CoV-2 genomes.
The reads that we get will not be aligned, but they will be based on the primers used in the sequencing step (it's unclear if we'll know which primer was used).

## Detecting VoC's

We have a detailed, fool-proof plan to detect the exact number of Variants of Concern that exist in the relevant population, along with the spatio-temporal spread of these variants, all while accurately determining the actual case counts of Covid-19.
This plan is as follows: 

**TODO.**


