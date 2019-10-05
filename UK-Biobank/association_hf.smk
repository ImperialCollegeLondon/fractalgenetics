#### example cluster call #####
# snakemake -s association_hf.smk --jobs 5000 --latency-wait 120 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

###### load config file #####
configfile: "config/config_hf_association.yaml"


##### target rules #####
rule all:
    input:
        # rule/combine.smk
        expand("{dir}/grs/all_samples_heart_failures_controls.txt",
            dir=config["dir"]),
        expand("{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.bed",
            dir=config["dir"]),
        expand("{dir}/grs/genotypes/genotypes_{type}_sigSNPs.bed",
            dir=config["dir"],
            type=["hf", "cad", "icm", "nicm", "sz_nicm", "controls"]),
        expand("{dir}/grs/genotypes/controls_{type}_sigSNPs.clean.bed",
            dir=config["dir"],
            type=["hf", "cad", "icm", "nicm", "sz_nicm"]),
        expand("{dir}/grs/genotypes/controls_{type}_covariates.txt",
            dir=config["dir"],
            type=["hf", "cad", "icm", "nicm", "sz_nicm"]),
        expand("{dir}/grs/genotypes/controls_{type}_sigSNPs.clean.assoc.logistic.all.sig",
            dir=config["dir"],
            type=["hf", "cad", "icm", "nicm", "sz_nicm"]),


rule all_samples:
    input:
        cov="{dir}/grs/prsice_covariates_heart_failures_ukb.txt",
    output:
        all="{dir}/grs/covariates_heart_failures_controls.txt",
        id="{dir}/grs/all_samples_heart_failures_controls.txt",
    shell:
        """
        awk '$2 != "NA"' {input.cov} > {output.all}
        cut -d ' '  -f 1 {output.all} |tail -n +2 | tr '\n' ' ' > {output.id}
        """


rule all_samples_genotypes:
    input:
        geno=expand("{geno}/ukb_imp_chr{chr}_v3.bgen",
            chr=range(1,23),
            geno=config["genodir"]),
        sampleids="{dir}/grs/all_samples_heart_failures_controls.txt",
        samples="{dir}/rawdata/ukb40616_imp_chr1_v3_s487317.sample",
        rsids=expand("{dir}/Pseudomultitrait_slices_sig5e08_qctool.IDs",
            dir=config['sigdir'],
            pheno=config['pheno']),
    params:
        genodir=config['genodir']
    output:
        bed="{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.bed",
        bim="{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.bim",
        fam="{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.fam",
    shell:
        """
        qctool -g {params.genodir}/ukb_imp_chr#_v3.bgen \
            -incl-rsids {input.rsids} \
            -ofiletype binary_ped \
            -s {input.samples} \
            -og {wildcards.dir}/grs/genotypes/all_samples_heart_failures_controls_tmp
        """

rule test_merge:
    input:
        bed="{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.bed",
        bim="{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.bim",
        fam="{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.fam",
    output:
        mis="{dir}/grs/genotypes/merge_test_sigSNPs.missnp",
    params:
    shell:
        """
        set +e
        plink --bfile {wildcards.dir}/grs/genotypes/all_samples_heart_failures_controls_tmp \
            --bmerge {wildcards.dir}/grs/genotypes/all_samples_heart_failures_controls_tmp \
            --out {wildcards.dir}/grs/genotypes/merge_test_sigSNPs
        """

rule extract_sig_genotypes_plink:
    input:
        bed="{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.bed",
        bim="{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.bim",
        fam="{dir}/grs/genotypes/all_samples_heart_failures_controls_tmp.fam",
        mis="{dir}/grs/genotypes/merge_test_sigSNPs.missnp",
        samples="{dir}/heart_failure_phenotypes/{type}_eid.csv",
    output:
        bed="{dir}/grs/genotypes/genotypes_{type}_sigSNPs.bed",
        bim="{dir}/grs/genotypes/genotypes_{type}_sigSNPs.bim",
        fam="{dir}/grs/genotypes/genotypes_{type}_sigSNPs.fam",
        samples="{dir}/grs/genotypes/genotypes_{type}_qctool.id",
    params:
    shell:
        """
        awk 'NR==FNR {{a[$2]=$0; next}} $1 in a {{print a[$1]}} ' {input.fam} \
            {input.samples} > {output.samples}
        plink --bfile {wildcards.dir}/grs/genotypes/all_samples_heart_failures_controls_tmp \
            --exclude {input.mis} \
            --keep {output.samples} \
            --make-bed \
            --allow-extra-chr \
            --allow-no-sex \
            --out {wildcards.dir}/grs/genotypes/genotypes_{wildcards.type}_sigSNPs
        """

rule merge:
    input:
        hf_bed="{dir}/genotypes_{type}_sigSNPs.bed",
        hf_fam="{dir}/genotypes_{type}_sigSNPs.fam",
        hf_bim="{dir}/genotypes_{type}_sigSNPs.bim",
        controls_bed="{dir}/genotypes_controls_sigSNPs.bed",
        controls_fam="{dir}/genotypes_controls_sigSNPs.fam",
        controls_bim="{dir}/genotypes_controls_sigSNPs.bim",
    output:
        both_bed="{dir}/controls_{type}_sigSNPs.bed",
        both_bim="{dir}/controls_{type}_sigSNPs.bim",
        both_fam="{dir}/controls_{type}_sigSNPs.fam"
    shell:
        """
        plink --bfile {wildcards.dir}/genotypes_{wildcards.type}_sigSNPs \
            --bmerge {wildcards.dir}/genotypes_controls_sigSNPs \
            --make-bed \
            --allow-extra-chr \
            --allow-no-sex \
            --make-pheno {input.hf_fam} '*' \
            --out {wildcards.dir}/controls_{wildcards.type}_sigSNPs
        """

rule covariates:
    input:
        fam="{dir}/grs/genotypes/controls_{type}_sigSNPs.fam",
        cov="{dir}/grs/covariates_heart_failures_controls.txt",
    output:
        "{dir}/grs/genotypes/controls_{type}_covariates.txt"
    shell:
        """
        awk 'NR==FNR {{first=$1; $1=""; a[first]=$0; next}} $2 in a {{print $1" "$2""a[$2]}}' \
            <(tr "," ' '< {input.cov}) {input.fam} > {output}
        """

rule missing:
    input:
        both_bed="{dir}/controls_{type}_sigSNPs.bed",
        both_bim="{dir}/controls_{type}_sigSNPs.bim",
        both_fam="{dir}/controls_{type}_sigSNPs.fam"
    output:
        missing="{dir}/controls_{type}_sigSNPs.missing",
        ids="{dir}/controls_{type}_sigSNPs.missing.ids",
        both_bed="{dir}/controls_{type}_sigSNPs.clean.bed",
        both_bim="{dir}/controls_{type}_sigSNPs.clean.bim",
        both_fam="{dir}/controls_{type}_sigSNPs.clean.fam"
    shell:
        """
        plink --bfile {wildcards.dir}/controls_{wildcards.type}_sigSNPs \
            --test-missing \
            --allow-extra-chr \
            --allow-no-sex \
            --pfilter 0.01 \
            --out {wildcards.dir}/controls_{wildcards.type}_sigSNPs
        tr -s ' ' < {output.missing} | sed 's/ //' | cut -d " " -f 2 | \
            tail -n +2 > {output.ids}
        plink --bfile {wildcards.dir}/controls_{wildcards.type}_sigSNPs \
            --allow-extra-chr \
            --allow-no-sex \
            --exclude {output.ids} \
            --make-bed \
            --out {wildcards.dir}/controls_{wildcards.type}_sigSNPs.clean
        """

rule association_cc:
    input:
        both_bed="{dir}/controls_{type}_sigSNPs.clean.bed",
        both_bim="{dir}/controls_{type}_sigSNPs.clean.bim",
        both_fam="{dir}/controls_{type}_sigSNPs.clean.fam",
        cov="{dir}/controls_{type}_covariates.txt",
    output:
        logist="{dir}/controls_{type}_sigSNPs.clean.assoc.logistic",
        perm="{dir}/controls_{type}_sigSNPs.clean.assoc.logistic.perm",
        combined="{dir}/controls_{type}_sigSNPs.clean.assoc.logistic.all"
    shell:
        """
        plink --bfile {wildcards.dir}/controls_{wildcards.type}_sigSNPs.clean \
            --logistic hide-covar perm \
            --covar {input.cov} \
            --allow-extra-chr \
            --allow-no-sex \
            --out {wildcards.dir}/controls_{wildcards.type}_sigSNPs.clean
        tr -s ' ' < {output.perm} | sed 's/ //' | tr ' ' '\t' | cut -f 3,4 | \
            paste {output.logist} -  > {output.combined}
        """



rule significant:
    input:
        cc="{dir}/controls_{type}_sigSNPs.clean.assoc.logistic.all",
        genotyped=expand('{dir}/ukb_cal_genome.bim',
            dir=config['genotyped'])
    params:
        sigloci=16
    output:
        sig="{dir}/controls_{type}_sigSNPs.clean.assoc.logistic.all.sig",
        geno="{dir}/controls_{type}_sigSNPs.clean.assoc.logistic.all.sig.genotyped"
    shell:
        """
        head -n 1 {input.cc} > {output.geno}
        awk 'FNR==NR && ($10*{params.sigloci}) < 0.05' {input.cc} > {output.sig}
        awk 'FNR==NR {{a[$2]=$0; next}} $2 in a {{print a[$2]}}' {output.sig} \
            {input.genotyped} >> {output.geno}
        """

