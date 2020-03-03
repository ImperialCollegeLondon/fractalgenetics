# how to run fro ./:
# snakemake -s association_discovery.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} \
# -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_association.yaml"

rule all:
    input:
        expand("{ukb}/gwas/{pheno}/bgenie_{analysis}_STRAIN_lm_st_chr{chr}.gz",
            ukb=config["ukbdir"],
            pheno=config['discovery'],
            analysis=['slices'],
            chr=range(1,23)),
        expand("{ukb}/gwas/{pheno}/Pseudomultitrait_{analysis}_STRAIN_sig5e08_ldFiltered.txt",
            ukb=config["ukbdir"],
            analysis=['slices'],
            pheno=config['discovery']),

rule generateSNPfiles:
    input:
        rsids="{dir}/maf0.001/{pheno}/ukb_imp_chr{chr}_v3_maf0.001.bim"
    output:
        "{dir}/maf0.001/{pheno}/ukb_imp_chr{chr}_v3_maf0.001.rsid"
    shell:
        "cut -f 2 {input.rsids} |sort|uniq > {output}"

rule generateLDtags:
    input:
        keep="{dir}/ancestry/{pheno}/European_FID_IID.keep",
        rsid="{dir}/maf0.001/{pheno}/ukb_imp_chr{chr}_v3_maf0.001.rsid",
        bed="{dir}/plink/{pheno}/ukb_imp_chr{chr}_v3.bed",
        bim="{dir}/plink/{pheno}/ukb_imp_chr{chr}_v3.bim",
        fam="{dir}/plink/{pheno}/ukb_imp_chr{chr}_v3.fam",
    params:
        r2=config['r2'],
        kb=config['kbwindow']
    output:
        "{dir}/tags/{pheno}/European_ukb_imp_chr{chr}_v3_maf{maf}_{kb}kb_r{r2}.tags.list"
    shell:
        "plink --bed {input.bed} \
            --fam {input.fam} \
            --bim {input.bim} \
            --keep {input.keep} \
            --extract {input.rsid} \
            --show-tags all \
            --tag-r2 {params.r2} --tag-kb {params.kb} \
            --out {wildcards.dir}/tags/{wildcards.pheno}/European_ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}_{params.kb}kb_r{params.r2}"

rule association:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/{{pheno}}/FD_{{analysis}}_bgenie.txt",
            ukb=config["ukbdir"]),
        covs=expand("{ukb}/phenotypes/{{pheno}}/FD_covariates_bgenie.txt",
            ukb=config["ukbdir"]),
        rsids=expand("{ukb}/maf0.001/{{pheno}}/ukb_imp_chr{{chr}}_v3_maf0.001.rsid",
                    ukb=config["ukbdir"])
    params:
        n=config["n"]
    output:
        "{dir}/gwas/{pheno}/bgenie_{analysis}_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --include_rsids  {input.rsids} \
            --pvals --exc_missing_inds \
            --thread {params.n} \
            --out {wildcards.dir}/gwas/{wildcards.pheno}/bgenie_{wildcards.analysis}_lm_st_chr{wildcards.chr}"

rule association_strain:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/{{pheno}}/FD_{{analysis}}_bgenie.txt",
            ukb=config["ukbdir"]),
        covs=expand("{ukb}/phenotypes/{{pheno}}/FD_covariates_STRAIN_bgenie.txt",
            ukb=config["ukbdir"]),
        rsids=expand("{ukb}/maf0.001/{{pheno}}/ukb_imp_chr{{chr}}_v3_maf0.001.rsid",
                    ukb=config["ukbdir"])
    params:
        n=config["n"]
    output:
        "{dir}/gwas/{pheno}/bgenie_{analysis}_STRAIN_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --include_rsids  {input.rsids} \
            --pvals --exc_missing_inds \
            --thread {params.n} \
            --out {wildcards.dir}/gwas/{wildcards.pheno}/bgenie_{wildcards.analysis}_STRAIN_lm_st_chr{wildcards.chr}"

rule results_strain:
    input:
        tags=expand("{{dir}}/tags/{{pheno}}/European_ukb_imp_chr{chr}_v3_maf{maf}_{kb}kb_r{r2}.tags.list",
            maf=config['maf'],
            r2=config['r2'],
            kb=config['kbwindow'],
            chr=range(1,23)),
        gwas=expand("{{dir}}/gwas/{{pheno}}/bgenie_{{analysis}}_STRAIN_lm_st_chr{chr}.gz",
            chr=range(1,23)),
        pheno="{dir}/phenotypes/{pheno}/FD_{analysis}_EUnorel.csv"
    output:
        "{dir}/gwas/{pheno}/bgenie_{analysis}_STRAIN_lm_pseudomt_qqplot.pdf",
        "{dir}/gwas/{pheno}/bgenie_{analysis}_STRAIN_lm_pseudomt_manhattanplot.pdf",
        "{dir}/gwas/{pheno}/bgenie_{analysis}_STRAIN_lm_st_genomewide.csv",
        "{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_STRAIN_sig5e08.txt",
        "{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_STRAIN_sig5e08_ldFiltered.txt",
    shell:
        "Rscript association/association-results.R \
            --pheno {input.pheno} \
            --cohort {wildcards.pheno} \
            --name {wildcards.analysis}_STRAIN \
            --directory {wildcards.dir} \
            --showProgress "

rule results:
    input:
        tags=expand("{{dir}}/tags/{{pheno}}/European_ukb_imp_chr{chr}_v3_maf{maf}_{kb}kb_r{r2}.tags.list",
            maf=config['maf'],
            r2=config['r2'],
            kb=config['kbwindow'],
            chr=range(1,23)),
        gwas=expand("{{dir}}/gwas/{{pheno}}/bgenie_{{analysis}}_lm_st_chr{chr}.gz",
            chr=range(1,23)),
        pheno="{dir}/phenotypes/{pheno}/FD_{analysis}_EUnorel.csv"
    output:
        "{dir}/gwas/{pheno}/bgenie_{analysis}_lm_pseudomt_qqplot.pdf",
        "{dir}/gwas/{pheno}/bgenie_{analysis}_lm_pseudomt_manhattanplot.pdf",
        "{dir}/gwas/{pheno}/bgenie_{analysis}_lm_st_genomewide.csv",
        "{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_sig5e08.txt",
        "{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_sig5e08_ldFiltered.txt",
    shell:
        "Rscript association/association-results.R \
            --pheno {input.pheno} \
            --cohort {wildcards.pheno} \
            --name {wildcards.analysis} \
            --directory {wildcards.dir} \
            --showProgress "

rule functional_enrichment:
    input:
        "{dir}/gwas/{pheno}/bgenie_{analysis}_lm_st_genomewide.csv",
    output:
        "{dir}/annotation/{pheno}/Functional_enrichment_{analysis}.pdf",
        "{dir}/annotation/{pheno}/Functional_enrichment_{analysis}_all.pdf",
    shell:
        "Rscript association/functional-enrichment.R \
            --name {wildcards.analysis} \
            --directory {wildcards.dir}/gwas/{wildcards.pheno} \
            --showProgress "

rule extract_genotypes:
    input:
        geno=expand("{geno}/ukb_imp_chr{chr}_v3.bgen",
            chr=range(1,23),
            geno=config["genodir"]),
        samples="{dir}/rawdata/ukb18545_imp_chr1_v3_s487378.sample",
        fdsamples="{dir}/phenotypes/{pheno}/FD_{analysis}_EUnorel.csv",
        sig="{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_sig5e08.txt"
    output:
        bimbam="{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_sig5e08_genotypes.dosage.gz",
        rsids="{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_sig5e08_qctool.IDs",
        toupdate="{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_sig5e08_qctool.toUpdate",
        sampleids="{dir}/phenotypes/{pheno}/FD_{analysis}_samples_qctool.IDs",
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

rule effect:
    input:
        "{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_sig5e08_ldFiltered.txt"
    output:
        "{dir}/gwas/{pheno}/Distribution_{analysis}_beta.pdf",
    shell:
        "Rscript association/effect-size-distribution.R \
            --directory {wildcards.dir}/gwas/{wildcards.pheno} \
            --showProgress "
rule mr:
    input:
        instruments="{dir}/gwas/{pheno}/Pseudomultitrait_slices_sig5e08_ldFiltered.txt",
        bimbam="{dir}/gwas/{pheno}/Pseudomultitrait_slices_sig5e08_genotypes.dosage.gz",
        oauth="{dir}/MR/mrbase.oauth"
    output:
        "{dir}/MR/{pheno}/MRbase_{analysis}.rds"
    shell:
        "Rscript association/mendelian-randomisation.R \
            --gwasdir {wildcards.dir}/gwas/{wildcards.pheno} \
            --mrdir {wildcards.dir}/MR/{wildcards.pheno} \
            --showProgress "

rule ldhub:
    input:
        "{dir}/gwas/{pheno}/bgenie_{analysis}_lm_st_genomewide.csv",
        "{dir}/ldhub/{pheno}/sumstats_{analysis}_MeanBasalFD.txt",
        "{dir}/ldhub/{pheno}/sumstats_{analysis}_MeanMidFD.txt",
        "{dir}/ldhub/{pheno}/sumstats_{analysis}_MeanApicalFD.txt",
    output:
        "{dir}/ldhub/{pheno}/Rg_{analysis}_all.csv"
    shell:
        "Rscript association/LDhub-results.R \
            --directory {wildcards.dir}/ldhub/{wildcards.pheno} \
            --showProgress "

rule gtex:
    input:
        "{dir}/gwas/{pheno}/bgenie_slices_lm_st_genomewide.csv",
        expand("{{dir}}/GTEX/{{pheno}}/expression_atlas-homo_sapiens_{loci}.tsv",
            loci=["C4A", "GOSR2", "MTSS1", "PDZRN3", "PKP1", "SAP30L",
                  "SSXP10", "TDH", "ZNF358"])
    output:
        "{dir}/GTEX/{pheno}/gTEX_geneexpression.pdf"
    shell:
        "Rscript association/gtex-expression.R \
            --prefix expression_atlas-homo_sapiens \
            --directory {wildcards.dir}/GTEX/{wildcards.pheno} \
            --showProgress "
