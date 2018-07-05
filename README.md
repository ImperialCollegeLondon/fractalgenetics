# Genetics of natural variation in left ventricular trabeculation phenotypes 

In human hearts, trabeculation of the left ventricle varies in pattern and strength. Both, naturally occuring differences
in healthy individuals across ethnicities and clinically relevant hypertrabeculation phenotypes (left-ventricular non-compaction)
observed. The genetics of the natural variation and clinically observed are still poorly understood.

We use cardiac magnetic resonance imaging and genotypes from the UK Biobank to study the genetics of natural variation
in left ventricular trabeculation phenotypes.

Left ventricular trabeculation phenotypes were automatically extracted via fractal analysis 
from cardiac magnetic resonance images ([FD estimation](https://github.com/UK-Digital-Heart-Project/AutoFD)).
This repository contains the work-flows for genotype processing and the genetic analysis. 

## conversion.smk and config_conversion.yaml
Work-flow for converting biobank genotype data in .bgen format to plink format, minor-allele frequency and ld-filtering
of variants. Application-specific (18545) filtering of samples from 500k genotypes in first step to speed up computation.
Parameters for filtering, file names and target directories are supplied in config_conversion.yaml.

## ancestry.smk and config_ancestry.yaml
Work-flow for estimating kinship and ancestry of sample cohort. Takes ld-pruned, maf-filtered files from conversion.smk.
For ancestry estimation, cohort genotypes are fused with HapMap genotypes of known ancestry (processed as described in [HapMap processing]()),
principal components computed and cohort samples within 1.5 times the maximum Euclidean distance of European Hapmap samples
to the center of the European Hapmap samples are selected as European via [selectPCA.R](https://github.com/HannahVMeyer/ukbb-fd/blob/master/selectPCA.R).
Parameters for filtering, file names and target directories are supplied in config_ancestry.yaml.

## cluster.json
Config file for snakemake job submission to lsf cluster system with conversion.smk and ancestry.smk rule-specific requirements.
