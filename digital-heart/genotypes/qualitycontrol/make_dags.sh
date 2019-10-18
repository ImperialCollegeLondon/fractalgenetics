snakemake -s genotypeQC.smk \
        --rulegraph \
        ~/data/digital-heart/genotype/QC/combined/gencall.combined.HapMapIII.eigenvec \
        ~/data/digital-heart/genotype/QC/combined/HVOL.gencall.combined.clean.related.eigenvec \
        ~/data/digital-heart/genotype/QC/combined/HVOL.gencall.combined.clean.related.bed |
            dot -Tpng > ~/analysis/fractalgenetics/digital-heart/dag/genotypes_qc_dag.png

