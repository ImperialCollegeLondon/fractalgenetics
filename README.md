[![DOI](https://zenodo.org/badge/166986710.svg)](https://zenodo.org/badge/latestdoi/166986710)

# Genomic analysis reveals a functional role for myocardial trabeculae in adults

The endocardial surface of the adult human hearts consists of a complex network of muscular trabeculae. Thought to be a remnant of embryonic development, it remains uknown why these complex structures still persistent in the adult organ. Here, we use population genomics, image-based intermediate phenotyping and *in silico* modelling to determine the effect of this complex cardiovascular trait on function.

The full study is available on [BioRxiv](https://www.biorxiv.org/content/10.1101/553651v1).

The following repository contains all data analysis and processing scripts applied in the study.

## Software requirements
Data processing and analysis is embedded into different [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (v5.1.5) workflows.
Genetics analyses require
plink [v1.9](https://www.cog-genomics.org/plink2) and [v2](https://www.cog-genomics.org/plink/2.0/) as well as [flashpca](https://github.com/gabraham/flashpca).
Data visualisation and processing is done in [R](https://www.r-project.org/), version 3.4.1.
Genetic association testing requires [Bgenie](https://jmarchini.org/bgenie/) v1.3.
Fractal analysis requires Matlab, [MATLAB Compiler Runtime (MCR) 2016b](https://uk.mathworks.com/products/compiler/matlab-runtime.html) and [Ghostscript](https://www.ghostscript.com/). Finite element modeling was conducted with [Abaqus/Standard, Simulia](https://www.3ds.com/products-services/simulia/products/abaqus/abaqusstandard/).

## Content
### 1. fractal-analysis
Fractal analysis tool for computation of left ventricular/right ventricular fractal characteristics, using LVSA series from cardiac magnetic resonance images.

### 2. automated-fractal-analysis
Automated fractal analysis of segmented cardiac images based on the [FracAnalyse](fractal-analysis) software. It uses pre-exisiting image segmentations to determine a region of interest within the myocardium for fractal analysis.

### 3. fractal-analysis-processing
Code for post-processing of fractal-analysis results including interpolation of FD values to a common number of slices across individuals, summary statistics of FD values per individual and a collection of functions for co-registration of myocardial and trabeculation outlines.

### 4. UK-Biobank
Trabeculation GWAS in UKB cohort (discovery cohort) containing analysis pipelines for genotype and phenotype processing,
GWAS, GWAS results processing, functional enrichement and Mendelian randomisation analyses.

### 5. digital-heart
Trabeculation GWAS in Digital-heart project cohort (validation cohort) containing analysis pipelines for genotype and phenotype processing, GWAS and GWAS results processing.

### 6. Finite element modelling
Finite element modeling input files (Abaqus Standard, SIMULIA, Dessault Systemes)for the simulation of 5 consecutive cardiac cycles of the left ventricle, starting from the ventricular reference configuration (zero ventricular pressure).







