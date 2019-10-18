snakemake -s association.smk \
        --rulegraph \
        ~/data/digital-heart/gwas/FD/bgenie_slices_lm_pseudomt_qqplot.pdf \
        ~/data/digital-heart/gwas/FD/Pseudomultitrait_slices_sig5e08_genotypes.dosage.gz |
            dot -Tpng > ~/analysis/fractalgenetics/digital-heart/dag/association_dag.png

