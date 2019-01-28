# how to run fro ./:
# snakemake -s association.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_association.yaml"

rule all:
    input:
        expand("{ukb}/gwas/bgenie_{analysis}_lm_st_chr{chr}.gz",
            ukb=config["ukbdir"],
            analysis=['slices'],
            chr=range(1,23)),
        expand("{ukb}/gwas/bgenie_{analysis}_lm_pseudomt_{plot}.pdf",
           analysis=['slices'],
           plot=['qqplot', 'manhattanplot'],
           ukb=config["ukbdir"]),
        expand("{ukb}/tags/European_ukb_imp_chr{chr}_v3_maf{maf}_{kb}kb_r{r2}.tags.list",
            ukb=config["ukbdir"],
            maf=config['maf'],
            kb=config['kbwindow'],
            r2=config['r2'],
            chr=range(1,23)),
        expand("{ukb}/gwas/Pseudomultitrait_slices_sig5e08_genotypes.dosage.gz",
            ukb=config["ukbdir"])

rule generateSNPfiles:
    input:
        rsids="{dir}/maf0.001/ukb_imp_chr{chr}_v3_maf0.001.bim"
    output:
        "{dir}/maf0.001/ukb_imp_chr{chr}_v3_maf0.001.rsid"
    shell:
        "cut -f 2 {input.rsids} |sort|uniq > {output}"

rule generateLDtags:
    input:
        keep="{dir}/ancestry/European_FID_IID.keep",
        rsid="{dir}/maf0.001/ukb_imp_chr{chr}_v3_maf0.001.rsid",
        bed="{dir}/plink/ukb_imp_chr{chr}_v3.bed",
        bim="{dir}/plink/ukb_imp_chr{chr}_v3.bim",
        fam="{dir}/plink/ukb_imp_chr{chr}_v3.fam",
    params:
        r2=config['r2'],
        kb=config['kbwindow']
    output:
        "{dir}/tags/European_ukb_imp_chr{chr}_v3_maf{maf}_{kb}kb_r{r2}.tags.list"
    shell:
        "plink --bed {input.bed} \
            --fam {input.fam} \
            --bim {input.bim} \
            --keep {input.keep} \
            --extract {input.rsid} \
            --show-tags all \
            --tag-r2 {params.r2} --tag-kb {params.kb} \
            --out {wildcards.dir}/tags/European_ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}_{params.kb}kb_r{params.r2}"

rule association:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/FD_{{analysis}}_bgenie.txt",
            ukb=config["ukbdir"]),
        covs=expand("{ukb}/phenotypes/FD_covariates_bgenie.txt",
            ukb=config["ukbdir"]),
        rsids=expand("{ukb}/maf0.001/ukb_imp_chr{{chr}}_v3_maf0.001.rsid",
                    ukb=config["ukbdir"])
    params:
        n=config["n"]
    output:
        "{dir}/gwas/bgenie_{analysis}_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --include_rsids  {input.rsids} \
            --pvals --exc_missing_inds \
            --thread {params.n} \
            --out {wildcards.dir}/gwas/bgenie_{wildcards.analysis}_lm_st_chr{wildcards.chr}"

rule results:
    input:
        tags=expand("{{dir}}/tags/European_ukb_imp_chr{chr}_v3_maf{maf}_{kb}kb_r{r2}.tags.list",
            maf=config['maf'],
            r2=config['r2'],
            kb=config['kbwindow'],
            chr=range(1,23)),
        gwas=expand("{{dir}}/gwas/bgenie_{{analysis}}_lm_st_chr{chr}.gz",
            chr=range(1,23)),
        pheno="{dir}/phenotypes/FD_{analysis}_EUnorel.csv"
    output:
        "{dir}/gwas/bgenie_{analysis}_lm_pseudomt_qqplot.pdf",
        "{dir}/gwas/bgenie_{analysis}_lm_pseudomt_manhattanplot.pdf",
    shell:
        "Rscript association/association_results.R \
            --pheno {input.pheno} \
            --name summary \
            --directory {wildcards.dir}/gwas \
            --showProgress "


rule extractSigGenotypes:
    input:
        geno=expand("{geno}/ukb_imp_chr{chr}_v3.bgen",
            chr=range(1,23),
            geno=config["genodir"]),
        samples="{dir}/rawdata/ukb18545_imp_chr1_v3_s487378.sample",
        fdsamples="{dir}/phenotypes/FD_slices_EUnorel.csv",
        sig="{dir}/gwas/Pseudomultitrait_slices_sig5e08.txt"
    output:
        bimbam="{dir}/gwas/Pseudomultitrait_slices_sig5e08_genotypes.dosage.gz",
        rsids="{dir}/gwas/Pseudomultitrait_slices_sig5e08_qctool.IDs",
        toupdate="{dir}/gwas/Pseudomultitrait_slices_sig5e08_qctool.toUpdate",
        sampleids="{dir}/phenotypes/FD_samples_qctool.IDs",
    params:
        genodir=config["genodir"]
    shell:
        """
        cut -d " " -f 2 {input.sig} | tail -n +2 |tr '\n' ' ' > {output.rsids}
        cut -d " " -f 1,2 {input.sig} > {output.toupdate}
        cut -d "," -f 1 {input.fdsamples} | tr '\n' ' ' > {output.sampleids}
        qctool -g {params.genodir}/ukb_imp_chr#_v3.bgen \
            -incl-rsids {output.rsids} -incl-samples {output.sampleids} \
            -og {output.bimbam} -s {input.samples}
        """

