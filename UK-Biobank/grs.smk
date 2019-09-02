# how to run from ./:
# snakemake -s grs.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name}
# -q {cluster.queue} -n {cluster.n} \
# -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_grs.yaml"

def assign_linker(wildcards):
    if wildcards.type == "fd" or wildcards.type == "ukb_replication":
        return "{}/rawdata/ukb18545_imp_chr1_v3_s487378.sample".format(wildcards.dir)
    else:
        return "{}/rawdata/ukb40616_imp_chr1_v3_s487317.sample".format(wildcards.dir)

def assign_replication(wildcards):
    if wildcards.cohort == "ukb":
        samples="{}/rawdata/ukb18545_imp_chr1_v3_s487378.sample".format(wildcards.dir)
        pheno="{}/phenotypes/190402_fractal_dimension_26k/FD_summary_bgenie.txt".format(wildcards.dir)
        cov="{}/phenotypes/190402_fractal_dimension_26k/FD_covariates_bgenie.txt".format(wildcards.dir)
    if wildcards.cohort == "dh":
        dd="/homes/hannah/data/digital-heart"
        samples="{}/genotype/imputation/combined/genotypes/gencall.combined.clean.related.chr1.sample".format(dd)
        pheno="{}/phenotype/FD/FD_summary_bgenie.txt".format(dd)
        cov="{}/phenotype/FD/FD_covariates_bgenie.txt".format(dd)
    return({"samples": samples, "pheno": pheno, "cov": cov)


rule all:
    input:
        expand("{ukb}/genotypes/ukb_cal_genome.bim",
             ukb=config["ukbdir"])
        expand("{ukb}/grs/{region}_grs_genotyped.summary",
             ukb=config["ukbdir"],
             region=['MeanBasalFD', 'MeanApicalFD', 'MeanMidFD']),
        expand("{ukb}/grs/prsice_pheno_replication_{cohort}.txt",
             ukb=config["ukbdir"],
             cohort=['ukb', 'dh'])

rule select_ukb_hf_phenotypes:
    input:
        ukbdata="{dir}/rawdata/ukb40616.txt",
    output:
        hf="{dir}/heart_failure_phenotypes/hf_eid.csv",
        nicm="{dir}/heart_failure_phenotypes/nicm_eid.csv",
        icm="{dir}/heart_failure_phenotypes/icm_eid.csv",
        cad="{dir}/heart_failure_phenotypes/cad_eid.csv",
        aragam_nicm="{dir}/heart_failure_phenotypes/aragam_nicm_eid.csv",
        sz_nicm="{dir}/heart_failure_phenotypes/sz_nicm_eid.csv",
        ids="{dir}/heart_failure_phenotypes/heart_failure_samples_id.rds"
    shell:
        """
        Rscript grs/select-hf-phenotypes.R \
            --outdir {wildcards.dir}/heart_failure_phenotypes \
            --data {input.ukbdata} \
        """

rule select_ukb_hf_covariates:
    input:
        ukbdata="{dir}/rawdata/ukb40616.txt",
        samples="{dir}/rawdata/ukb40616_imp_chr1_v3_s487317.sample",
        relatedness="{dir}/rawdata/ukb40616_rel_s488282.dat",
    output:
        covs="{dir}/heart_failure_phenotypes/heart_failure_covariates.csv",
        overview="{dir}/heart_failure_phenotypes/heart_failure_samples_overview.csv"
    shell:
        """
        Rscript grs/select-hf-covariates.R \
            --outdir {wildcards.dir}/heart_failure_phenotypes \
            --data {input.ukbdata} \
            --samples {input.samples} \
            --relatedness {input.relatedness}
        """

rule collect_bim_rsids:
    input:
        rsids=expand("{{dir}}/genotypes/ukb_cal_chr{chr}_v2.bim.gz",
                    chr=range(1,23))
    output:
        rsid="{dir}/genotypes/ukb_cal_genome.bim",
    shell:
        """
        zcat {input.rsids} > {output.rsid}
        """

rule format_discovery_prsice:
    input:
        samples="{dir}/rawdata/ukb18545_imp_chr1_v3_s487378.sample",
        gwas="{dir}/gwas/180628_fractal_dimension/bgenie_summary_lm_st_genomewide.csv",
        pheno="{dir}/phenotypes/180628_fractal_dimension/FD_summary_bgenie.txt",
        cov="{dir}/phenotypes/180628_fractal_dimension/FD_covariates_bgenie.txt",
        rsid="{dir}/genotypes/ukb_cal_genome.bim",
    params:
    output:
        cov="{dir}/grs/prsice_covariates_discovery_ukb.txt",
        basal="{dir}/grs/prsice_pheno_MeanBasalFD_discovery_ukb.txt",
        mid="{dir}/grs/prsice_pheno_MeanMidFD_discovery_ukb.txt",
        apical="{dir}/grs/prsice_pheno_MeanApicalFD_discovery_ukb.txt",
    shell:
        """
        Rscript grs/format-discovery-prsice.R \
            --directory {wildcards.dir} \
            --samples {input.samples} \
            --pheno {input.pheno} \
            --covariates {input.cov} \
            --gwas {input.gwas} \
            --genotypedids {input.rsid} \
            --showProgress
        """

rule construct_grs:
    input:
        bim=expand("{{dir}}/genotypes/ukb_cal_chr{chr}_v2.bim",
            chr=range(1,23)),
        bed=expand("{{dir}}/genotypes/ukb_cal_chr{chr}_v2.bed",
            chr=range(1,23)),
        fam="{dir}/rawdata/ukb18545_cal_chr1_v2_s488282.fam",
        gwas="{dir}/grs/prsice_association_summary_{region}.txt",
        pheno="{dir}/grs/prsice_pheno_summary_{region}.txt",
        covs="{dir}/grs/prsice_covariates_summary.txt",
        heartsamples="{dir}/heart_phenotypes/fd_eid.csv",
    params:
        geno=config["genodir"],
    output:
        results="{dir}/grs/{region}_grs_genotyped.summary"
    shell:
        """
        Rscript ~/software/PRSice/PRSice.R --dir {wildcards.dir} \
            --prsice PRSice \
            --base {input.gwas} \
            --A1 a_0 \
            --A2 a_1 \
            --bp pos \
            --chr chr \
            --pvalue p \
            --snp rsid \
            --stat beta \
            --beta \
            --target {wildcards.dir}/genotypes/ukb_cal_chr#_v2,{input.fam} \
            --keep {input.heartsamples} \
            --type bed \
            --pheno {input.pheno} \
            --ignore-fid \
            --binary-target F \
            --allow-inter \
            --perm 10000 \
            --cov {input.covs} \
            --cov-factor sex \
            --out {wildcards.dir}/grs/{wildcards.region}_grs_genotyped \
            --print-snp \
            --thread 8 \
            --memory 15Gb \
            --seed 101 \
            """

rule format_replication_prsice:
    input:
        unpack(assign_replication)
    output:
        pheno="{dir}/grs/prsice_pheno_replication_{cohort}.txt",
        cov="{dir}/grs/prsice_covariates_replication_{cohort}.txt",
    shell:
        """
        Rscript grs/format-replication-prsice.R \
            --outdir {wildcards.dir}/grs \
            --cohort {wildcards.cohort} \
            --samples {input.samples} \
            --pheno {input.pheno} \
            --cov {input.cov} \
            --showProgress
        """

rule replication_grs:
    input:
        bim=expand("{{dir}}/genotypes/ukb_cal_chr{chr}_v2.bim",
            chr=range(1,23)),
        bed=expand("{{dir}}/genotypes/ukb_cal_chr{chr}_v2.bed",
            chr=range(1,23)),
        fam="{dir}/rawdata/ukb18545_cal_chr1_v2_s488282.fam",
        gwas="{dir}/grs/prsice_association_summary_{region}.txt",
        pheno="{dir}/grs/prsice_pheno_summary_{region}.txt",
        covs="{dir}/grs/prsice_covariates_summary.txt",
        heartsamples="{dir}/heart_phenotypes/fd_eid.csv",
    params:
        geno=config["genodir"],
    output:
        results="{dir}/grs/{region}_grs_genotyped.summary"
    shell:
        """
        Rscript ~/software/PRSice/PRSice.R --dir {wildcards.dir} \
            --prsice PRSice \
            --base {input.gwas} \
            --A1 a_0 \
            --A2 a_1 \
            --bp pos \
            --chr chr \
            --pvalue p \
            --snp rsid \
            --stat beta \
            --beta \
            --target {wildcards.dir}/genotypes/ukb_cal_chr#_v2,{input.fam} \
            --keep {input.heartsamples} \
            --type bed \
            --pheno {input.pheno} \
            --ignore-fid \
            --binary-target F \
            --allow-inter \
            --perm 10000 \
            --cov {input.covs} \
            --cov-factor sex \
            --out {wildcards.dir}/grs/{wildcards.region}_grs_genotyped \
            --print-snp \
            --thread 8 \
            --memory 15Gb \
            --seed 101 \
            """

rule format_hf_prsice:
    input:
        hf="{dir}/heart_failure_phenotypes/hf_eid.csv",
        nicm="{dir}/heart_failure_phenotypes/nicm_eid.csv",
        icm="{dir}/heart_failure_phenotypes/icm_eid.csv",
        cad="{dir}/heart_failure_phenotypes/cad_eid.csv",
        aragam_nicm="{dir}/heart_failure_phenotypes/aragam_nicm_eid.csv",
        sz_nicm="{dir}/heart_failure_phenotypes/sz_nicm_eid.csv",
        ids="{dir}/heart_failure_phenotypes/heart_failure_samples_id.rds"
    output:
        pheno="{dir}/grs/prsice_heart_failures_ukb.txt",
        cov="{dir}/grs/prsice_covariates_heart_failures_ukb.txt",
        overview="{dir}/grs/heart_failures_ukb_overview.csv",
    shell:
        """
        Rscript grs/format-hf-prsice.R \
            --outdir {wildcards.dir}/grs \
            --samples {input.samples} \
            --hf {input.hf} \
            --icm {input.icm} \
            --cad {input.cad} \
            --nicm {input.nicm} \
            --sznicm {input.sznicm} \
            --aragamnicm {input.aragamnicm} \
            --idoverview {input.ids} \
            --cov {input.cov} \
            --showProgress
        """

rule hf_grs:
    input:
        bim=expand("{{dir}}/genotypes/ukb_cal_chr{chr}_v2.bim",
            chr=range(1,23)),
        bed=expand("{{dir}}/genotypes/ukb_cal_chr{chr}_v2.bed",
            chr=range(1,23)),
        fam="{dir}/rawdata/ukb18545_cal_chr1_v2_s488282.fam",
        gwas="{dir}/grs/prsice_association_summary_{region}.txt",
        pheno="{dir}/grs/prsice_pheno_summary_{region}.txt",
        covs="{dir}/grs/prsice_covariates_summary.txt",
        heartsamples="{dir}/heart_phenotypes/fd_eid.csv",
    params:
        geno=config["genodir"],
    output:
        results="{dir}/grs/{region}_grs_genotyped.summary"
    shell:
        """
        Rscript ~/software/PRSice/PRSice.R --dir {wildcards.dir} \
            --prsice PRSice \
            --base {input.gwas} \
            --A1 a_0 \
            --A2 a_1 \
            --bp pos \
            --chr chr \
            --pvalue p \
            --snp rsid \
            --stat beta \
            --beta \
            --target {wildcards.dir}/genotypes/ukb_cal_chr#_v2,{input.fam} \
            --keep {input.heartsamples} \
            --type bed \
            --pheno {input.pheno} \
            --ignore-fid \
            --binary-target F \
            --allow-inter \
            --perm 10000 \
            --cov {input.covs} \
            --cov-factor sex \
            --out {wildcards.dir}/grs/{wildcards.region}_grs_genotyped \
            --print-snp \
            --thread 8 \
            --memory 15Gb \
            --seed 101 \
            """
