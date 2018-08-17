# Genetics of natural variation in left ventricular trabeculation phenotypes

In human hearts, trabeculation of the left ventricle varies in pattern and strength. Both, naturally occurring differences
in healthy individuals across ethnicities and clinically relevant hypertrabeculation phenotypes (left-ventricular non-compaction)
observed. The genetics of the natural variation and clinically observed are still poorly understood.

We use cardiac magnetic resonance imaging and genotypes from the UK Biobank to study the genetics of natural variation
in left ventricular trabeculation phenotypes.

Left ventricular trabeculation phenotypes were automatically extracted via fractal analysis 
from cardiac magnetic resonance images ([FD estimation](https://github.com/UK-Digital-Heart-Project/AutoFD)).
This repository contains the work-flows for phenotype and genotype processing and the genetic analysis. 

## Software requirements
Data processing and analysis is embedded into different [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (v5.1.5) workflows.
[genotypes.smk](https://github.com/HannahVMeyer/ukbb-fd/genotypes.smk) and [ancestry.smk](https://github.com/HannahVMeyer/ukbb-fd/ancestry.smk) require
plink [v1.9](https://www.cog-genomics.org/plink2) and [v2](https://www.cog-genomics.org/plink/2.0/) as well as [flashpca](https://github.com/gabraham/flashpca).
Data visualisation and processing is done in [R](https://www.r-project.org/), version 3.4.1.
Genetic association testing requires [Bgenie](https://jmarchini.org/bgenie/) v1.3.


## Analysis steps
### 1. genotypes.smk
Work-flow for converting biobank genotype data in .bgen format to plink format, minor-allele frequency and ld-filtering
of variants. Application-specific (18545) filtering of samples from 500k genotypes in first step to speed up computation.
Parameters for filtering, file names and target directories are supplied in config/config_conversion.yaml.

### 2. ancestry.smk
Work-flow for estimating kinship and ancestry of sample cohort. Takes ld-pruned, maf-filtered files from [genotypes.smk](https://github.com/HannahVMeyer/ukbb-fd/genotypes.smk).
For ancestry estimation, cohort genotypes are fused with HapMap genotypes of known ancestry (processed as described in [HapMap processing](https://www.ncbi.nlm.nih.gov/pubmed/21085122)),
principal components computed and cohort samples within 1.5 times the maximum Euclidean distance of European Hapmap samples
to the centre of the European Hapmap samples are selected as European via [selectPCA.R](https://github.com/HannahVMeyer/ukbb-fd/ancestry/selectPCA.R).
Parameters for filtering, file names and target directories are supplied in config/config_ancestry.yaml. Scripts called by [ancestry.smk](https://github.com/HannahVMeyer/ukbb-fd/ancestry.smk) can be found
in [ancestry](https://github.com/HannahVMeyer/ukbb-fd/ancestry).

### 3. phenotypes.smk
Work-flow for processing ukbb covariates and ukbb-derived FD phenotypes obtained via [FD estimation](https://github.com/UK-Digital-Heart-Project/AutoFD).
UKBB covariates are downloaded, decrypted and converted to R-readable formats via ukbtools. It uses ancestry and relatedness information generated via [ancestry.smk](https://github.com/HannahVMeyer/ukbb-fd/ancestry.smk) 
and filters FD phenotypes and covariates for unrelated samples of European ancestry (via [preparePheno.r](https://github.com/HannahVMeyer/ukbb-fd/phenotypes/preparePheno.r). It tests covariates
for association with the FD phenotypes and saves relevant data in bgenie format for associating testing via [association.smk](https://github.com/HannahVMeyer/ukbb-fd/association.smk).
Scripts called by [phenotypes.smk](https://github.com/HannahVMeyer/ukbb-fd/phenotypes.smk) can be found in [phenotypes](https://github.com/HannahVMeyer/ukbb-fd/phenotypes).

### 4. association.smk
Scripts called by [association.smk](https://github.com/HannahVMeyer/ukbb-fd/association) can be found in [association](https://github.com/HannahVMeyer/ukbb-fd/association).

## config
Config files for snakemake files job and config files for job submission to lsf cluster system with rule-specific requirements ([cluster.json](https://github.com/HannahVMeyer/ukbb-fd/config/cluster.json)).
