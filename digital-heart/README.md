# Trabeculation GWAS in Digital-Heart project cohort

### 1. genotypes
Work-flow for genotype quality control, phasing and imputation of genotypes called with the gencall algorithm.
Parameters, file names and target directories are supplied in the respective subdirecties config/config_conversion.yaml files.
Scripts called by the snakemake workflows can be found in the respective scripts subdirectories.
Requires plink [v1.9](https://www.cog-genomics.org/plink2), [shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
and [impute2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html).

### 2. phenotypes
Work-flow for processing co-variates and FD phenotypes obtained via [FD estimation](automated-fractal-analysis).
FD phenotypes and covariates are filtered for unrelated samples of European ancestry (via [preparePheno.r](digital-heart/phenotypes/preparePheno.r). It tests covariates
for association with the FD phenotypes and saves relevant data in bgenie format for associating testing via [association.smk](digital-heart/association.smk).
Parameters, file names and target directories are supplied in the [config file](digital-heart/phenotypes/config/config_conversion.yaml).

### 3. association
Work-flow for GWAS with [Bgenie](https://jmarchini.org/bgenie/).
Scripts called by [association.smk](digital-heart/association/association.smk) can be found in [association](digital-heart/association/scripts).
Parameters, file names and target directories are supplied in the [config file](digital-heart/association/config/config_conversion.yaml).

