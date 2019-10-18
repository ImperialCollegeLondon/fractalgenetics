snakemake -s format.smk \
        --rulegraph \
        ~/data/digital-heart/genotype/QC/combined/gencall.combined.bim \
        ~/data/digital-heart/genotype/QC/combined/same_genotype_probes.txt |
            dot -Tpng > ~/analysis/fractalgenetics/digital-heart/dag/genotypes_formating_dag.png


