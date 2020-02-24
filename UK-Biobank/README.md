# Trabeculation GWAS in UKB cohort

### 1. [genotypes.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/genotypes.smk)
Work-flow for converting biobank genotype data in .bgen format to plink format,
for computing minor-allele frequency and ld-filtering of variants.
Application-specific (18545) filtering of samples from 500k genotypes in first
step to speed up computation. Parameters for filtering, file names and target
directories are supplied in config/config_conversion.yaml. One workflow for
discovery and replication samples.
Requires plink [v1.9](https://www.cog-genomics.org/plink2) and
[v2](https://www.cog-genomics.org/plink/2.0/).

<p align="center"> 
<img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/genotypes_dag.png" height="400">
</p>

### 2. [ancestry.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/ancestry.smk)
Work-flow for estimating kinship and ancestry of discovery and replication
cohort. Takes ld-pruned,maf-filtered files from [genotypes.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/genotypes.smk).
For ancestry estimation, cohort genotypes are fused with HapMap genotypes of
known ancestry (processed as described in [HapMap processing](https://www.ncbi.nlm.nih.gov/pubmed/21085122)),
principal components computed and cohort samples within 1.5 times the maximum
Euclidean distance of European Hapmap samples to the centre of the European
Hapmap samples are selected as European via [selectPCA.R](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/ancestry/selectPCA.R).
Parameters for filtering, file names and target directories are supplied in
config/config_ancestry.yaml. Scripts called by [ancestry.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/ancestry).
Requires plink [v1.9](https://www.cog-genomics.org/plink2) and
[flashpca](https://github.com/gabraham/flashpca).

<p align="center"> 
<img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/ancestry_dag.png" height="700">
</p>

### 3. [phenotypes_discovery.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/phenotypes_discovery.smk) and [phenotypes_replication.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/phenotypes_replication.smk)
Work-flow for processing ukbb covariates and ukbb-derived FD phenotypes obtained
via [FD estimation](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/automated-fractal-analysis).
UKBB covariates are downloaded, decrypted and converted to R-readable formats
via ukbtools. It uses ancestry and relatedness information generated via
[ancestry.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/ancestry.smk) 
and filters FD phenotypes and covariates for unrelated samples of European
ancestry (via [phenotypes](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/phenotypes). It tests covariates
for association with the FD phenotypes and saves relevant data in bgenie format
for associating testing via [association](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/association).
Scripts called by [phenotypes_discovery.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/phenotypes_discovery.smk) and [phenotypes_replication.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/phenotypes_replication.smk) can be found in
[phenotypes](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/phenotypes). Two independent workflows for discovery
(left) and replication (right) analysis.

<p align="center"> 
<img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/phenotypes_discovery_dag.png" height="100">
  <img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/phenotypes_replication_dag.png" height="100">
</p>

### 4. [association_discovery.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/association_discovery.smk) and [association_replication.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/association_replication.smk)
Work-flow for GWAS with [Bgenie](https://jmarchini.org/bgenie/),
Mendelian randomisation (MR), functional enrichement and genetic correlation
analysis. MR analysis uses the R package [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR).
Functional enrichement analysis depends on [GARFIELD2](https://www.ebi.ac.uk/birney-srv/GARFIELD)
Genetic correlation analysis relied on
[LDhub](http://ldsc.broadinstitute.org/ldhub/) and requires sign in with user
account. Scripts called by [association_discovery.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/association_discovery.smk) and [association_replication.smk](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/association_replication.smk) can be
found in [association](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/association). Two independent workflows for
discovery (left) and replication (right) analysis.

<p align="center"> 
<img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/association_discovery_dag.png" height="300">
  <img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/association_replication_dag.png" height="300">
</p>

### 5. [pv-loops](https://github.com/ImperialCollegeLondon/fractalgenetics/tree/master/UK-Biobank/pv-loops)
Analysis and example data for visualisation of cardiac mechanics during systole
and diastole by pressure-volume loop analysis describing the relationship between
left ventricular pressure and left ventricular volume at multiple time points
during a complete cardiac cycle.

## [config](https://github.com/ImperialCollegeLondon/fractalgenetics/tree/master/UK-Biobank/config)
Config files for snakemake files job and config files for job submission to lsf cluster system with rule-specific requirements ([cluster.json](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/config/cluster.json)).

