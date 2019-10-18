snakemake -s impute.smk \
        --rulegraph \
        ~/data/digital-heart/genotype/imputation/combined/formated/gencall.combined.clean.related.genome.qc.bgen \
        ~/data/digital-heart/genotype/imputation/combined/counts/gencall.combined.clean.related.SNPsPerChr.txt |
            dot -Tpng > ~/analysis/fractalgenetics/digital-heart/dag/genotypes_impute_dag.png


