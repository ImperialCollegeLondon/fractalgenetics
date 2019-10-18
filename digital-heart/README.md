# Trabeculation and DCM GWAS in Digital-Heart project cohort

### 1. genotypes
Work-flow for genotype formating (left), quality control (middle), phasing and imputation (right) of genotypes called with the gencall algorithm. Additional workflows for check on batch effects via gwas with batch as phenotype in [dag](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/digital-heart/dag). 
Parameters, file names and target directories are supplied in the respective subdirecties config/config_conversion.yaml files.
Scripts called by the snakemake workflows can be found in the respective scripts subdirectories.
Requires plink [v1.9](https://www.cog-genomics.org/plink2), [shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
and [impute2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html).

<p align="center"> 
  <img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/digital-heart/dag/genotypes_formating_dag.png" height="250">
<img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/digital-heart/dag/genotypes_qc_dag.png" height="400">
 <img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/digital-heart/dag/genotypes_impute_dag.png" height="400">
</p>

### 2. phenotypes
Work-flow for processing co-variates and FD phenotypes obtained via [FD estimation](automated-fractal-analysis).
FD phenotypes and covariates are filtered for unrelated samples of European ancestry (via [preparePheno.r](digital-heart/phenotypes/preparePheno.r). It tests covariates
for association with the FD phenotypes and saves relevant data in bgenie format for associating testing via [association.smk](digital-heart/association.smk).
Parameters, file names and target directories are supplied in the [config file](digital-heart/phenotypes/config/config_phenotypes.yaml).

<p align="center"> 
  <img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/digital-heart/dag/phenotypes_dag.png" height="50">
</p>

### 3. FD association
Work-flow for FD GWAS with [Bgenie](https://jmarchini.org/bgenie/).
Scripts called by [association.smk](digital-heart/association/association.smk) can be found in [association](digital-heart/association/scripts).
Parameters, file names and target directories are supplied in the [config file](digital-heart/association/config/config_association.yaml).

<p align="center"> 
  <img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/digital-heart/dag/association_dag.png" height="150">
</p>

### 3. DCM phenotypes and association
Work-flow for DCM GWAS with [plink]([v1.9](https://www.cog-genomics.org/plink2)).
Scripts called by [dcm.smk](digital-heart/DCM/dcm.smk) can be found in [DCM](digital-heart/dcm).
Parameters, file names and target directories are supplied in the [config file](digital-heart/association/config/config_conversion.yaml).

<p align="center"> 
  <img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/digital-heart/dag/dcm_dag.png" height="300">
</p>
