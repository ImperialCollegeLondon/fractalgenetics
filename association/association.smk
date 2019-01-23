# how to run from ./:
# snakemake -s association.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_association.yaml"

rule all:
    input:
        expand("{dir}/bgenie_{analysis}_lm_st_chr{chr}.gz",
            dir=config["gwasdir"],
            analysis=['slices'],
            chr=range(1,23)),
        expand("{dir}/{analysis}_{suffix}",
            dir=config["gwasdir"],
            analysis=['slices'],
            suffix=["concordance.txt", "concordance.pdf",
                "concordance_summary.txt"]),
        expand("{dir}/{analysis}_{suffix}",
            dir=config["gwasdir"],
            analysis=['slices'],
            suffix=["lm_pseudomt_qqplot.pdf", "lm_pseudomt_manhattanplot.pdf")

rule slices:
    input:
        geno=expand("{geno}/{name}.chr{{chr}}.qc.bgen",
            name=config['name'],
            geno=config["genodir"]),
        pheno=expand("{pheno}/FD_{{analysis}}_bgenie.txt",
            pheno=config["phenodir"]),
        covs=expand("{pheno}/FD_covariates_bgenie.txt",
            pheno=config["phenodir"])
    params:
        n=config["n"]
    output:
        "{dir}/bgenie_{analysis}_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --pvals --exc_missing_inds \
            --thread {params.n} \
            --out {wildcards.dir}/bgenie_{wildcards.analysis}_lm_st_chr{wildcards.chr}"

rule results:
    input:
        gwas=expand("{{dir}}/bgenie_slices_lm_st_chr{chr}.gz",
            chr=range(1,23))
        pheno=expand("{pheno}/FD_{{analysis}}_bgenie.txt",
            pheno=config["phenodir"]),
    output:
        "{dir}/gwas/bgenie_{analysis}_lm_pseudomt_qqplot.pdf",
        "{dir}/gwas/bgenie_{analysis}_lm_pseudomt_manhattanplot.pdf",
    shell:
        "Rscript association/association_results.R \
            --pheno {input.pheno} \
            --name {wildcards.analysis} \
            --directory {wildcards.dir}/association \
            --showProgress "

rule concordance:
    input:
        bgenie=expand("{{dir}}/bgenie_slices_lm_st_chr{chr}.gz",
            chr=range(1,23))
    output:
        "{dir}/{analysis}_concordance.txt",
        "{dir}/{analysis}_concordance.pdf",
        "{dir}/{analysis}_concordance_summary.txt"
    params:
        ukbdir=config['ukbdir']
    wildcard_constraints:
        analysis="\w+"
    shell:
        "Rscript scripts/concordance.R \
            --name {wildcards.analysis} \
            --directory {wildcards.dir} \
            --ukbdir {params.ukbdir} \
            --showProgress "


