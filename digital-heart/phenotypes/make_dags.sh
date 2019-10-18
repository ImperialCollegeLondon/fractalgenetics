snakemake -s phenotypes.smk \
        --rulegraph \
        ~/data/digital-heart/phenotypes/FD/FD_slices_bgenie.txt |
            dot -Tpng > ~/analysis/fractalgenetics/digital-heart/dag/phenotypes_dag.png

