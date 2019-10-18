# Trabeculation GWAS in UKB cohort

### 1. genotypes.smk
Work-flow for converting biobank genotype data in .bgen format to plink format,
minor-allele frequency and ld-filtering of variants.
Application-specific (18545) filtering of samples from 500k genotypes in first
step to speed up computation. Parameters for filtering, file names and target
directories are supplied in config/config_conversion.yaml.
Requires plink [v1.9](https://www.cog-genomics.org/plink2) and
[v2](https://www.cog-genomics.org/plink/2.0/).

<p align="center"> 
<img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/genotypes_dag.png" height="400">
</p>

### 2. ancestry.smk
Work-flow for estimating kinship and ancestry of sample cohort. Takes ld-pruned, maf-filtered files from [genotypes.smk](UK-Biobank/genotypes.smk).
For ancestry estimation, cohort genotypes are fused with HapMap genotypes of known ancestry (processed as described in [HapMap processing](https://www.ncbi.nlm.nih.gov/pubmed/21085122)),
principal components computed and cohort samples within 1.5 times the maximum Euclidean distance of European Hapmap samples
to the centre of the European Hapmap samples are selected as European via [selectPCA.R](UK-Biobank/ancestry/selectPCA.R).
Parameters for filtering, file names and target directories are supplied in config/config_ancestry.yaml. Scripts called by [ancestry.smk](UK-Biobank/ancestry).
Requires plink [v1.9](https://www.cog-genomics.org/plink2) and [flashpca](https://github.com/gabraham/flashpca).

<p align="center"> 
<img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/ancestry_dag.png" height="700">
</p>

### 3. phenotypes.smk
Work-flow for processing ukbb covariates and ukbb-derived FD phenotypes obtained via [FD estimation](automated-fractal-analysis).
UKBB covariates are downloaded, decrypted and converted to R-readable formats via ukbtools. It uses ancestry and relatedness information generated via [ancestry.smk](https://github.com/HannahVMeyer/ukbb-fd/ancestry.smk) 
and filters FD phenotypes and covariates for unrelated samples of European ancestry (via [preparePheno.r](UK-Biobank/phenotypes/preparePheno.r). It tests covariates
for association with the FD phenotypes and saves relevant data in bgenie format for associating testing via [association.smk](UK-Biobank/association.smk).
Scripts called by [phenotypes.smk](UK-Biobank/phenotypes.smk) can be found in [phenotypes](UK-Biobank/phenotypes).

<p align="center"> 
<img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/phenotypes_discovery_dag.png" height="100">
  <img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/phenotypes_replication_dag.png" height="100">
</p>

### 4. association.smk
Work-flow for GWAS with [Bgenie](https://jmarchini.org/bgenie/), Mendelian randomisation (MR) and functional enrichement analysis. MR analysis uses the R package [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR).
Functional enrichement analysis depends on [GARFIELD2](https://www.ebi.ac.uk/birney-srv/GARFIELD)
Scripts called by [association.smk](UK-Biobank/association.smk) can be found in [association](UK-Biobank/association).

<p align="center"> 
<img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/association_discovery_dag.png" height="300">
  <img src="https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/UK-Biobank/dag/association_replication_dag.png" height="300">
</p>

## config
Config files for snakemake files job and config files for job submission to lsf cluster system with rule-specific requirements ([cluster.json](UK-Biobank/config/cluster.json)).

