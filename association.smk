# how to run fro ./:
# snakemake -s association.smk --jobs 5000 --latency-wait 30 --cluster-config cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete

configfile: "config_association.yaml"

rule all:
    input:
        expand("{ukb}/gwas/bgenie_lm_st_chr{chr}.gz",
            ukb=config["ukbdir"],
            chr=range(1,23))

rule generateSNPfiles:
    input:
        rsids="{dir}/maf0.001/ukb_imp_chr{chr}_v3_maf0.001.bim"
    output:
        "{dir}/maf0.001/ukb_imp_chr{chr}_v3_maf0.001.rsid"
    shell:
        "cut -f 2 {input.rsids} > {output}"

rule bgenie:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/FD_phenotypes_bgenie.txt",
            ukb=config["ukbdir"]),
        covs=expand("{ukb}/phenotypes/FD_covariates_bgenie.txt",
            ukb=config["ukbdir"]),
        rsids=expand("{ukb}/maf0.001/ukb_imp_chr{{chr}}_v3_maf0.001.rsid",
                    ukb=config["ukbdir"])
    output:
        "{dir}/gwas/bgenie_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --include_rsids  {input.rsids} \
            --pvals --exc_missing_inds \
            --scale_phenotypes \
            --out {dir}/gwas/bgenie_lm_st_chr{chr}"
