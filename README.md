[![DOI](https://zenodo.org/badge/166986710.svg)](https://zenodo.org/badge/latestdoi/166986710)

# Genomic analysis reveals a functional role for myocardial trabeculae in adults

The endocardial surface of the adult human hearts consists of a complex network of muscular trabeculae. Thought to be a remnant of embryonic development, it remains uknown why these complex structures still persistent in the adult organ. Here, we use population genomics, image-based intermediate phenotyping and *in silico* modelling to determine the effect of this complex cardiovascular trait on function.

The full study is available at: BioRxiv link here

The following repository contains all data analysis and processing scripts applied in the study.

## Software requirements
Data processing and analysis is embedded into different [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (v5.1.5) workflows.
Genetics analyses require
plink [v1.9](https://www.cog-genomics.org/plink2) and [v2](https://www.cog-genomics.org/plink/2.0/) as well as [flashpca](https://github.com/gabraham/flashpca).
Data visualisation and processing is done in [R](https://www.r-project.org/), version 3.4.1.
Genetic association testing requires [Bgenie](https://jmarchini.org/bgenie/) v1.3.

## Content
### 1. fractal-analysis
### 2. automated-fractal-analysis
### 3. fractal-analysis-processing
Code for post-processing of fractal-analysis results including interpolation of FD values to a common number of slices across individuals, summary statistics of FD values per individual and a collection of functions for co-registration of myocardial and trabeculation outlines.

### 4. UK-Biobank
Trabeculation GWAS in UKB cohort (discovery cohort) containing analysis pipelines for genotype and phenotype processing,
GWAS, GWAS results processing, functional enrichement and Mendelian randomisation analyses.

### 5. digital-heart
Trabeculation GWAS in Digital-heart project cohort (validation cohort) containing analysis pipelines for genotype and phenotype processing, GWAS and GWAS results processing.






