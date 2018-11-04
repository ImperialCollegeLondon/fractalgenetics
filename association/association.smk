# how to run fro ./:
# snakemake -s association.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_association.yaml"

def filesLDSC(wildcards):
    if wildcards.type == "summary":
        fdlist = ['MeanGlobalFD', 'MeanBasalFD', 'MeanApicalFD', 'MeanMidFD',
        'MaxMidFD','MaxBasalFD', 'MaxApicalFD']
        ldfiles = ['{}/gwas/ldsc_summary_{}.sumstats.gz'.format(
            wildcards.dir, fd) for fd in fdlist]
    if wildcards.type == "slices":
        ldfiles = ['{}/gwas/ldsc_slices_Slice_{}.sumstats.gz'.format(
            wildcards.dir, n) for n in range(1,11)]
    return ldfiles


rule all:
    input:
        expand("{dir}/bgenie_{name}_lm_st_chr{chr}.gz",
            dir=config["gwasdir"],
            name=['summary', 'slices'],
            chr=range(1,23)),
        expand("{dir}/sumstats_{type}_pseudomt.txt",
            type=['summary', 'slices'],
            dir=config["gwasdir"]),
        #expand("{ukb}/gwas/ldsc_{type}.sumstats.gz",
        #    type=['summary_{}'.format(x) for x in ['MeanGlobalFD',
        #        'MeanBasalFD', 'MeanApicalFD', 'MeanMidFD', 'MaxMidFD',
        #        'MaxBasalFD', 'MaxApicalFD']],
        #    ukb=config["ukbdir"]),
        #expand("{ukb}/gwas/ldsc_{type}_fd.log",
        #    type=['summary', 'slices', 'summary_pseudomt'],
        #    ukb=config["ukbdir"])



rule summary:
    input:
        geno=expand("{geno}/{name}.chr{{chr}}.qc.bgen",
            name=config['name'],
            geno=config["genodir"]),
        pheno=expand("{pheno}/FD_phenotypes_bgenie.txt",
            pheno=config["phenodir"]),
        covs=expand("{pheno}/FD_covariates_bgenie.txt",
            pheno=config["phenodir"])
    params:
        n=config["n"]
    output:
        "{dir}/bgenie_summary_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --pvals --exc_missing_inds \
            --scale_phenotypes \
            --thread {params.n} \
            --out {wildcards.dir}/bgenie_summary_lm_st_chr{wildcards.chr}"


rule slices:
    input:
        geno=expand("{geno}/{name}.chr{{chr}}.qc.bgen",
            name=config['name'],
            geno=config["genodir"]),
        pheno=expand("{pheno}/FD_slices_bgenie.txt",
            pheno=config["phenodir"]),
        covs=expand("{pheno}/FD_covariates_bgenie.txt",
            pheno=config["phenodir"])
    params:
        n=config["n"]
    output:
        "{dir}/bgenie_slices_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --pvals --exc_missing_inds \
            --scale_phenotypes \
            --thread {params.n} \
            --out {wildcards.dir}/bgenie_slices_lm_st_chr{wildcards.chr}"

rule summaryResults:
    input:
        pheno=expand("{pheno}/FD_phenotypes.csv",
            pheno=config["phenodir"]),
        bgenie=expand("{{dir}}/bgenie_summary_lm_st_chr{chr}.gz",
            chr=range(1,23))
    output:
        "{dir}/bgenie_summary_lm_pseudomt_qqplot.pdf",
        "{dir}/bgenie_summary_lm_pseudomt_manhattanplot.pdf",
        "{dir}/sumstats_summary_pseudomt.txt"
    shell:
        "Rscript scripts/association_results.R \
            --pheno {input.pheno} \
            --name summary \
            --directory {wildcards.dir} \
            --showProgress "

rule slicesResults:
    input:
        pheno=expand("{pheno}/FD_slices.csv",
            pheno=config["phenodir"]),
        bgenie=expand("{{dir}}/bgenie_slices_lm_st_chr{chr}.gz",
            chr=range(1,23))
    output:
        "{dir}/bgenie_slices_lm_pseudomt_qqplot.pdf",
        "{dir}/bgenie_slices_lm_pseudomt_manhattanplot.pdf",
        "{dir}/sumstats_slices_pseudomt.txt"
    shell:
        "Rscript scripts/association_results.R \
            --pheno {input.pheno} \
            --name slices \
            --directory {wildcards.dir} \
            --showProgress "

