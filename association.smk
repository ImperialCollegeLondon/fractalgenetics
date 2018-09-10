# how to run fro ./:
# snakemake -s association.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete

configfile: "config/config_association.yaml"

rule all:
    input:
        expand("{ukb}/gwas/bgenie_{name}_lm_pseudomt_qqplot.pdf",
            ukb=config["ukbdir"],
            name=['summary_EDV', 'summary', 'slices_EDV', 'slices'])

rule generateSNPfiles:
    input:
        rsids="{dir}/maf0.001/ukb_imp_chr{chr}_v3_maf0.001.bim"
    output:
        "{dir}/maf0.001/ukb_imp_chr{chr}_v3_maf0.001.rsid"
    shell:
        "cut -f 2 {input.rsids} > {output}"

rule summary_EDV:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/FD_phenotypes_bgenie.txt",
            ukb=config["ukbdir"]),
        covs=expand("{ukb}/phenotypes/FD_covariates_EDV_bgenie.txt",
            ukb=config["ukbdir"]),
        rsids=expand("{ukb}/maf0.001/ukb_imp_chr{{chr}}_v3_maf0.001.rsid",
                    ukb=config["ukbdir"])
    params:
        n=config["n"]
    output:
        "{dir}/gwas/bgenie_summary_EDV_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --include_rsids  {input.rsids} \
            --pvals --exc_missing_inds \
            --scale_phenotypes \
            --thread {params.n} \
            --out {wildcards.dir}/gwas/bgenie_summary_EDV_lm_st_chr{wildcards.chr}"

rule summary:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/FD_phenotypes_bgenie.txt",
            ukb=config["ukbdir"]),
        covs=expand("{ukb}/phenotypes/FD_covariates_bgenie.txt",
            ukb=config["ukbdir"]),
        rsids=expand("{ukb}/maf0.001/ukb_imp_chr{{chr}}_v3_maf0.001.rsid",
                    ukb=config["ukbdir"])
    params:
        n=config["n"]
    output:
        "{dir}/gwas/bgenie_summary_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --include_rsids  {input.rsids} \
            --pvals --exc_missing_inds \
            --scale_phenotypes \
            --thread {params.n} \
            --out {wildcards.dir}/gwas/bgenie_summary_lm_st_chr{wildcards.chr}"

rule slices_EDV:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/FD_slices_bgenie.txt",
            ukb=config["ukbdir"]),
        covs=expand("{ukb}/phenotypes/FD_covariates_EDV_bgenie.txt",
            ukb=config["ukbdir"]),
        rsids=expand("{ukb}/maf0.001/ukb_imp_chr{{chr}}_v3_maf0.001.rsid",
                    ukb=config["ukbdir"])
    params:
        n=config["n"]
    output:
        "{dir}/gwas/bgenie_slices_EDV_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --include_rsids  {input.rsids} \
            --pvals --exc_missing_inds \
            --scale_phenotypes \
            --thread {params.n} \
            --out {wildcards.dir}/gwas/bgenie_slices_EDV_lm_st_chr{wildcards.chr}"

rule slices:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/FD_phenotypes_bgenie.txt",
            ukb=config["ukbdir"]),
        covs=expand("{ukb}/phenotypes/FD_covariates_bgenie.txt",
            ukb=config["ukbdir"]),
        rsids=expand("{ukb}/maf0.001/ukb_imp_chr{{chr}}_v3_maf0.001.rsid",
                    ukb=config["ukbdir"])
    params:
        n=config["n"]
    output:
        "{dir}/gwas/bgenie_slices_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --include_rsids  {input.rsids} \
            --pvals --exc_missing_inds \
            --scale_phenotypes \
            --thread {params.n} \
            --out {wildcards.dir}/gwas/bgenie_slices_lm_st_chr{wildcards.chr}"

rule summaryResults:
    input:
        "{dir}/phenotypes/FD_phenotypes_EUnorel.csv"
    output:
        "{dir}/gwas/bgenie_{name}_lm_pseudomt_qqplot.pdf",
        "{dir}/gwas/bgenie_{name}_lm_pseudomt_manhattanplot.pdf"
    shell:
        "Rscript association/association_results.R \
            --pheno {input} \
            --name {wildcards.name} \
            --directory {wildcards.dir}/gwas \
            --showProgress "
