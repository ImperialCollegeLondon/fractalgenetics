snakemake -s dcm.smk \
        --rulegraph \
        ~/data/digital-heart/gwas/DCM/HVOL.DCM.FD.Pseudomultitrait_slices_sig5e08.sigSNPs.clean.assoc.logistic.all \
        ~/data/ukbb/ukb-hrt/DCM/FDAlongHeart_DCM_UKB_all_slices9.pdf |
            dot -Tpng > ~/analysis/fractalgenetics/digital-heart/dag/dcm_dag.png


