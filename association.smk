# how to run fro ./:
# snakemake -s association.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output}
# -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_association.yaml"

def filesLDSC(wildcards):
    if wildcards.type == "summary":
        fdlist = ['globalFD', 'meanBasalFD', 'meanApicalFD', 'maxBasalFD', 'maxApicalFD']
        ldfiles = ['{}/gwas/ldsc_summary_{}.sumstats.gz'.format(
            wildcards.dir, fd) for fd in fdlist]
    if wildcards.type == "slices":
        ldfiles = ['{}/gwas/ldsc_slices_Slice_{}.sumstats.gz'.format(
            wildcards.dir, n) for n in range(1,11)]
    return ldfiles


rule all:
    input:
        expand("{ukb}/gwas/bgenie_{name}_lm_st_chr{chr}.gz",
            ukb=config["ukbdir"],
            name=['summary', 'slices'],
            chr=range(1,23)),
        expand("{ukb}/gwas/ldsc_{type}.sumstats.gz",
            type=['summary_{}'.format(x) for x in ['globalFD', 'meanBasalFD',
                'meanApicalFD', 'maxBasalFD', 'maxApicalFD']],
            ukb=config["ukbdir"]),
        expand("{ukb}/gwas/ldsc_{type}_fd.log",
            type=['summary', 'slices', 'summary_pseudomt'],
            ukb=config["ukbdir"])

rule generateSNPfiles:
    input:
        rsids="{dir}/maf0.001/ukb_imp_chr{chr}_v3_maf0.001.bim"
    output:
        "{dir}/maf0.001/ukb_imp_chr{chr}_v3_maf0.001.rsid"
    shell:
        "cut -f 2 {input.rsids} > {output}"


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


rule slices:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/FD_slices_bgenie.txt",
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
        "{dir}/gwas/bgenie_summary_lm_pseudomt_qqplot.pdf",
        "{dir}/gwas/bgenie_summary_lm_pseudomt_manhattanplot.pdf"
    shell:
        "Rscript association/association_results.R \
            --pheno {input} \
            --name summary \
            --directory {wildcards.dir}/gwas \
            --showProgress "

rule slicesResults:
    input:
        "{dir}/phenotypes/FD_slices_EUnorel.csv"
    output:
        "{dir}/gwas/bgenie_slices_lm_pseudomt_qqplot.pdf",
        "{dir}/gwas/bgenie_slices_lm_pseudomt_manhattanplot.pdf"
    shell:
        "Rscript association/association_results.R \
            --pheno {input} \
            --name slices \
            --directory {wildcards.dir}/gwas \
            --showProgress "


rule ldScoreFormating:
    input:
        snplist="{}/eur_w_ld_chr/w_hm3.snplist".format(config['LDdir']),
        bgenie="{dir}/gwas/sumstats_{type}.txt"
    output:
        munge="{dir}/gwas/ldsc_{type}.sumstats.gz"
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        munge_sumstats.py \
            --sumstats {input.bgenie} \
            --maf-min 0.05 \
            --out {wildcards.dir}/gwas/ldsc_{wildcards.type} \
            --merge-alleles {input.snplist}
        """

rule ldScoreRegression:
    input:
        filesLDSC
    output:
        ldsc="{dir}/gwas/ldsc_{type}_fd.log"
    params:
        euroLD="{}/eur_w_ld_chr/".format(config['LDdir']),
        files=lambda wildcards, input: ",".join(map(str, input))
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        ldsc.py \
            --rg {params.files} \
            --ref-ld-chr {params.euroLD} \
            --w-ld-chr {params.euroLD} \
            --out {wildcards.dir}/gwas/ldsc_{wildcards.type}_fd
        python association/ldsc_parser.py -f {output} -o {wildcards.dir}/gwas \
            -n {wildcards.type}
        """

rule ldScoreRegressionMultitrait:
    input:
        "{dir}/gwas/ldsc_{type}.sumstats.gz"
    output:
        ldsc="{dir}/gwas/ldsc_{type}_fd.log"
    params:
        euroLD="{}/eur_w_ld_chr/".format(config['LDdir'])
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        ldsc.py \
            --h2 {input} \
            --ref-ld-chr {params.euroLD} \
            --w-ld-chr {params.euroLD} \
            --out {wildcards.dir}/gwas/ldsc_{wildcards.type}_fd
        python association/ldsc_parser.py -f {output} -o {wildcards.dir}/gwas \
            -n {wildcards.type}
        """
