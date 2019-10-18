snakemake -s control_plink.smk \
        --rulegraph \
        ~/data/digital-heart/genotype/control/plink_HVOL.gencall.combined.clean.related_lm_qqplot.pdf |
            dot -Tpng > ~/analysis/fractalgenetics/digital-heart/dag/genotypes_control_genotyped_dag.png

snakemake -s control_impute.smk \
        --rulegraph \
        ~/data/digital-heart/genotype/control/bgenie_HVOL.gencall.combined.clean.related_lm_qqplot.pdf |
            dot -Tpng >  ~/analysis/fractalgenetics/digital-heart/dag/genotypes_control_impute_dag.png


