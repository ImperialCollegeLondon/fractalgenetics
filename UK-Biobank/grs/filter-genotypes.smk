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

rule all:
    input:
        expand("{ukb}/maf0.001/{pheno}/ukb_imp_genome_v3_maf0.001.rsid",
            ukb=config["ukbdir"],
            pheno=config["pheno"]),
        #expand("{ukb}/grs/genotypes/genotypes_{type}.{fileformat}",
        #    fileformat=["bed", "dosage.gz", "bgen"],
        #    ukb=config["ukbdir"],
        #    type=['fd', 'hf', 'icm', 'sz_nicm', 'nicm', 'cad',
        #          'ukb_replication'])
        expand("{ukb}/grs/{region}_grs_genotyped.summary",
             ukb=config["ukbdir"],
             region=['MeanBasalFD', 'MeanApicalFD', 'MeanMidFD']),
        expand("{ukb}/genotypes/ukb_cal_genome.bim",
             ukb=config["ukbdir"])

rule select_ukb_hf_phenotypes:
    input:
        ukbdata="{dir}/rawdata/ukb40616.txt",
    output:
        hf="{dir}/heart_failure_phenotypes/hf_eid.csv",
        nicm="{dir}/heart_failure_phenotypes/nicm_eid.csv",
        icm="{dir}/heart_failure_phenotypes/icm_eid.csv",
        cad="{dir}/heart_failure_phenotypes/cad_eid.csv",
        aragam_nicm="{dir}/heart_failure_phenotypes/aragam_nicm_eid.csv",
        sz_nicm="{dir}/heart_failure_phenotypes/sz_nicm_eid.csv"
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
        Rscript grs/select-hf-phenotypes.R \
            --outdir {wildcards.dir}/heart_failure_phenotypes \
            --data {input.ukbdata} \
            --samples {input.samples} \
            --relatedness {input.relatedness}
        """

rule format_data_prsice:
    input:
        hf="{dir}/heart_failure_phenotypes/hf_eid.csv",
        nicm="{dir}/heart_failure_phenotypes/nicm_eid.csv",
        icm="{dir}/heart_failure_phenotypes/icm_eid.csv",
        cad="{dir}/heart_failure_phenotypes/cad_eid.csv",
        aragam_nicm="{dir}/heart_failure_phenotypes/aragam_nicm_eid.csv",
        sz_nicm="{dir}/heart_failure_phenotypes/sz_nicm_eid.csv",
        hf_covs="{dir}/heart_failure_phenotypes/heart_failure_covariates.csv",
        samples18545="{dir}/rawdata/ukb18545_imp_chr1_v3_s487378.sample",
        samples40616="{dir}/rawdata/ukb40616_imp_chr1_v3_s487317.sample",
        samplesDH=expand("{dir}/genotype/imputation/combined/genotypes/gencall.combined.clean.related.chr1.sample",
            dir=config['dhdir']),
        discoveryfile="{dir}/phenotypes/180628_fractal_dimension/FD_summary_bgenie.txt",
        cov_discoveryfile="{dir}/phenotypes/180628_fractal_dimension/FD_covariates_bgenie.txt",
        replicationfile="{dir}/phenotypes/190402_fractal_dimension_26k/FD_summary_bgenie.txt",
        cov_replicationfile="{dir}/phenotypes/190402_fractal_dimension_26k/FD_covariates_bgenie.txt",
        dhfile=expand("{dir}/phenotype/FD/FD_summary_bgenie.txt",
            dir=config['dhdir']),
        cov_dhfile=expand("{dir}/phenotype/FD/FD_covariates_bgenie.txt",
            dir=config['dhdir'])
    params:
        cohort=config["pheno"]
    output:
        pheno="{dir}/grs/prsice_hffailures_ukb.txt",
        pheno="{dir}/grs/hffailures_ukb_overview.txt",
        pheno="{dir}/grs/hffailures_ukb_qc_overview.txt",
        pheno="{dir}/grs/prsice_covariates_replication_dh.txt",
        pheno="{dir}/grs/prsice_covariates_replication_ukb.txt",
        pheno="{dir}/grs/prsice_covariates_discovery_ukb.txt",
    shell:
        """
        Rscript grs/format-data-prsice.R \
            --directory {wildcards.dir} \
            --cohort {params.cohort} \
            --sample18545 {input.sample18545} \
            --samples40616 {input.samples40616} \
            --samplesDH {input.samplesDH} \
            --discoveryfile {input.discoveryfile} \
            --cov_discoveryfile {input.cov_discoveryfile} \
            --replicationfile {input.replicationfile} \
            --cov_replicationfile {input.covreplicationfile} \
            --dhfile {input.dhfile} \
            --cov_dhfile {input.cov_dhfile} \
            --showProgress
        """

rule format_rsids:
    input:
        rsids=expand("{ukb}/maf0.001/{{pheno}}/ukb_imp_chr{chr}_v3_maf0.001.rsid",
                    ukb=config["ukbdir"],
                    chr=range(1,23))
    output:
        rsid="{dir}/maf0.001/{pheno}/ukb_imp_genome_v3_maf0.001.rsid",
        qctool="{dir}/maf0.001/{pheno}/ukb_imp_genome_v3_maf0.001_qctool.rsid"
    shell:
        """
        cat {input.rsids} > {output.rsid}
        tr '\n' ' ' < {output.rsid} > {output.qctool}
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

rule collect_sampleids:
    input:
        ukb_discovery="{dir}/phenotypes/180628_fractal_dimension/FD_slices_EUnorel.csv",
        ukb_replication="{dir}/phenotypes/190402_fractal_dimension_26k/FD_slices_EUnorel.csv",
        dh_replication=expand("{dhdir}/phenotype/FD/FD_slices.csv",
            dhdir=config['dhdir'])
    output:
        ukb_discovery="{dir}/heart_phenotypes/fd_eid.csv",
        ukb_replication="{dir}/heart_phenotypes/ukb_replication_eid.csv",
        dh_replication="{dir}/heart_phenotypes/dh_replication_eid.csv",
    shell:
        """
        cut -d "," -f 1 {input.ukb_discovery} | tail -n+2 > {output.ukb_discovery}
        cut -d "," -f 1 {input.ukb_replication} | tail -n+2 > {output.ukb_replication}
        cut -d "," -f 1 {input.dh_replication} | tail -n+2 > {output.dh_replication}
        """

rule format_sampleids:
    input:
        heartsamples="{dir}/heart_phenotypes/{type}_eid.csv",
    output:
        sampleids="{dir}/grs/genotypes/samples_{type}_qctool.IDs",
    shell:
        """
        tr '\n' ' ' < {input.heartsamples} > {output.sampleids}
        """

rule extract_genotypes_bgen:
    input:
        geno=expand("{geno}/ukb_imp_chr{chr}_v3.bgen",
            chr=range(1,23),
            geno=config["genodir"]),
        sampleids="{dir}/grs/genotypes/samples_{type}_qctool.IDs",
        rsids=expand("{{dir}}/maf0.001/{pheno}/ukb_imp_genome_v3_maf0.001_qctool.rsid",
            pheno=config['pheno']),
        samples=assign_linker
    output:
        plink="{dir}/grs/genotypes/genotypes_{type}.bgen",
    params:
        genodir=config["genodir"]
    shell:
        """
        qctool -g {params.genodir}/ukb_imp_chr#_v3.bgen \
            -incl-samples {input.sampleids} \
            -incl-rsids {input.rsids} \
            -ofiletype bgen \
            -og {wildcards.dir}/grs/genotypes/genotypes_{wildcards.type} \
            -s {input.samples}
        """
rule extract_genotypes_plink:
    input:
        geno=expand("{geno}/ukb_imp_chr{chr}_v3.bgen",
            chr=range(1,23),
            geno=config["genodir"]),
        sampleids="{dir}/grs/genotypes/samples_{type}_qctool.IDs",
        rsids=expand("{{dir}}/maf0.001/{pheno}/ukb_imp_genome_v3_maf0.001_qctool.rsid",
            pheno=config['pheno']),
        samples=assign_linker
    output:
        plink="{dir}/grs/genotypes/genotypes_{type}.bed",
    params:
        genodir=config["genodir"]
    shell:
        """
        qctool -g {params.genodir}/ukb_imp_chr#_v3.bgen \
            -incl-samples {input.sampleids} \
            -incl-rsids {input.rsids} \
            -ofiletype binary_ped \
            -og {wildcards.dir}/grs/genotypes/genotypes_{wildcards.type} \
            -s {input.samples}
        """

rule extract_genotypes_bimbam:
    input:
        geno=expand("{geno}/ukb_imp_chr{chr}_v3.bgen",
            chr=range(1,23),
            geno=config["genodir"]),
        rsids=expand("{{dir}}/maf0.001/{pheno}/ukb_imp_genome_v3_maf0.001_qctool.rsid",
            pheno=config['pheno']),
        sampleids="{dir}/grs/genotypes/samples_{type}_qctool.IDs",
        samples=assign_linker
    output:
        bimbam="{dir}/grs/genotypes/genotypes_{type}.dosage.gz",
    params:
        genodir=config["genodir"]
    shell:
        """
        qctool -g {params.genodir}/ukb_imp_chr#_v3.bgen \
            -incl-samples {input.sampleids} \
            -incl-rsids {input.rsids} \
            -ofiletype dosage \
            -og {wildcards.dir}/grs/genotypes/genotypes_{wildcards.type} \
            -s {input.samples}
        """

rule grs:
    input:
        #geno=expand("{geno}/ukb_imp_chr{chr}_v3.bgen",
        #    chr=range(1,23),
        #    geno=config["genodir"]),
        fam="{dir}/rawdata/ukb18545_cal_chr1_v2_s488282.fam",
        bim=expand("{{dir}}/genotypes/ukb_cal_chr{chr}_v2.bim",
            chr=range(1,23)),
        bed=expand("{{dir}}/genotypes/ukb_cal_chr{chr}_v2.bed",
            chr=range(1,23)),
        gwas="{dir}/grs/prsice_association_summary_{region}.txt",
        pheno="{dir}/grs/prsice_pheno_summary_{region}.txt",
        covs="{dir}/grs/prsice_covariates_summary.txt",
        rsids="{dir}/grs/MeanBasalFD_grs.valid",
        heartsamples="{dir}/heart_phenotypes/fd_eid.csv",
        #rsids=expand("{{dir}}/maf0.001/{pheno}/ukb_imp_genome_v3_maf0.001.rsid",
        #    pheno=config['pheno']),
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
