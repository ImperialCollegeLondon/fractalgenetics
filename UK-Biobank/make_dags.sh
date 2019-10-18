snakemake -s ancestry.smk \
        --rulegraph \
        ~/data/ukbb/ukb-hrt/ancestry/180628_fractal_dimension/ukb_imp_genome_v3_maf0.1.pruned.pca \
        ~/data/ukbb/ukb-hrt/ancestry/180628_fractal_dimension/ukb_imp_genome_v3_maf0.1.pruned.European.pca \
        ~/data/ukbb/ukb-hrt/ancestry/180628_fractal_dimension/HapMapIII_CGRCh37_180628_fractal_dimension_pca.png \
        ~/data/ukbb/ukb-hrt/ancestry/180628_fractal_dimension/ukb_imp_genome_v3_maf0.1.pruned.kinship.rel |
            dot -Tpng > ancestry_dag.png

snakemake -s genotypes.smk \
        --rulegraph \
        ~/data/ukbb/ukb-hrt/maf0.1/180628_fractal_dimension/ukb_imp_genome_v3_maf0.1.pruned.bed \
        ~/data/ukbb/ukb-hrt/maf0.1/190402_fractal_dimension_26k/ukb_imp_genome_v3_maf0.1.pruned.bed |
            dot -Tpng > genotypes_dag.png

snakemake -s phenotypes_discovery.smk \
        --rulegraph \
        ~/data/ukbb/ukb-hrt/phenotypes/180628_fractal_dimension/FD_slices_bgenie.txt |
            dot -Tpdf > phenotypes_discovery_dag.pdf

snakemake -s phenotypes_replication.smk \
        --rulegraph \
        ~/data/ukbb/ukb-hrt/phenotypes/190402_fractal_dimension_26k/FD_slices_bgenie.txt |
            dot -Tpdf > phenotypes_replication_dag.pdf

snakemake -s association_discovery.smk \
        --rulegraph \
        ~/data/ukbb/ukb-hrt/gwas/180628_fractal_dimension/Pseudomultitrait_slices_sig5e08_ldFiltered.txt \
        ~/data/ukbb/ukb-hrt/gwas/180628_fractal_dimension/Pseudomultitrait_slices_sig5e08_genotypes.dosage.gz \
        ~/data/ukbb/ukb-hrt/gwas/180628_fractal_dimension/Distribution_slices_beta.pdf \
        ~/data/ukbb/ukb-hrt/MR/180628_fractal_dimension/MRbase_summary.rds \
        ~/data/ukbb/ukb-hrt/ldhub/180628_fractal_dimension/Rg_summary_all.csv \
        ~/data/ukbb/ukb-hrt/GTEX/180628_fractal_dimension/gTEX_geneexpression.pdf \
        ~/data/ukbb/ukb-hrt/annotation/180628_fractal_dimension/Functional_enrichment_summary.pdf |
            dot -Tpdf > association_discovery_dag.pdf

snakemake -s association_replication.smk \
        --rulegraph \
        ~/data/ukbb/ukb-hrt/gwas/190402_fractal_dimension_26k/Replication_slices_association_significant_in_discovery.txt |
            dot -Tpdf > association_replication_dag.pdf
