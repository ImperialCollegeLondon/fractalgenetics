from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

###### load config file #####
configfile: "config/config_dcm.yaml"


##### target rules #####
rule all:
    input:
        # rule/combine.smk
        expand("{dir}/{diagnosis}.FD.{type}.sigSNPs.dosage",
            diagnosis=['HVOL', 'DCM'],
            type=['Pseudomultitrait_Slices_sig5e08'],
            #type=['Pseudomultitrait_Slices_sig5e08', 'GCC', 'random'],
            dir=config["dir"]),
        expand("{dir}/HVOL.DCM.FD.{type}.sigSNPs.{suffix}",
            suffix=['bed', 'bim', 'fam'],
            type=['Pseudomultitrait_Slices_sig5e08'],
            #type=['Pseudomultitrait_Slices_sig5e08', 'GCC', 'random'],
            dir=config["dir"]),
        expand("{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.model.perm.{suffix}",
            suffix=['sig', 'sig.genotyped'],
            type=['Pseudomultitrait_Slices_sig5e08'],
            #type=['Pseudomultitrait_Slices_sig5e08', 'GCC', 'random'],
            dir=config["dir"])

rule extract:
    input:
        bgen=expand("{dir}/formated/{name}.genome.qc.bgen",
            name=config['name'],
            dir=config['imputedir']),
        sample=expand("{dir}/genotypes/{name}.chr1.sample",
            name=config['name'],
            dir=config['imputedir']),
        rsids=expand('{dir}/{{type}}_qctool.IDs',
            dir=config['sigdir']),
        update=expand('{dir}/{{type}}_qctool.toUpdate',
            dir=config['sigdir']),
        diagnosis=expand('{dir}/European.{{diagnosis}}.gencall.combined.qctool.IDs',
            dir=config['genodir'])
    output:
        bed="{dir}/{diagnosis}.FD.{type}.sigSNPs.bed",
        bim="{dir}/{diagnosis}.FD.{type}.sigSNPs.bim",
        fam="{dir}/{diagnosis}.FD.{type}.sigSNPs.fam",
        bimbam="{dir}/{diagnosis}.FD.{type}.sigSNPs.dosage",
        stats="{dir}/{diagnosis}.FD.{type}.sigSNPs.stats.txt"
    wildcard_constraints:
        diagnosis="\w+"
    shell:
        """
        qctool -g {input.bgen} -og {output.bimbam} -s {input.sample} \
            -incl-samples {input.diagnosis} -incl-rsids {input.rsids}
        qctool -g {input.bgen} -og {output.bed} -s {input.sample} \
            -incl-samples {input.diagnosis} -incl-rsids {input.rsids} \
            -ofiletype binary_ped
        plink --bfile {wildcards.dir}/{wildcards.diagnosis}.FD.{wildcards.type}.sigSNPs \
            --allow-extra-chr \
            --update-chr {input.update} 1 2 \
            --make-bed \
            --out {wildcards.dir}/{wildcards.diagnosis}.FD.{wildcards.type}.sigSNPs
        qctool -g {input.bgen} -s {input.sample} \
            -incl-samples {input.diagnosis} -incl-rsids {input.rsids} \
            -snp-stats -osnp {output.stats}
        """

rule merge:
    input:
        dcm_bed="{dir}/DCM.FD.{type}.sigSNPs.bed",
        dcm_bim="{dir}/DCM.FD.{type}.sigSNPs.bim",
        dcm_fam="{dir}/DCM.FD.{type}.sigSNPs.fam",
        hvol_bed="{dir}/HVOL.FD.{type}.sigSNPs.bed",
        hvol_bim="{dir}/HVOL.FD.{type}.sigSNPs.bim",
        hvol_fam="{dir}/HVOL.FD.{type}.sigSNPs.fam"
    output:
        both_bed="{dir}/HVOL.DCM.FD.{type}.sigSNPs.bed",
        both_bim="{dir}/HVOL.DCM.FD.{type}.sigSNPs.bim",
        both_fam="{dir}/HVOL.DCM.FD.{type}.sigSNPs.fam"
    shell:
        """
        plink --bfile {wildcards.dir}/HVOL.FD.{wildcards.type}.sigSNPs \
            --bmerge {wildcards.dir}/DCM.FD.{wildcards.type}.sigSNPs \
            --make-bed \
            --allow-extra-chr \
            --allow-no-sex \
            --make-pheno {input.dcm_fam} '*' \
            --out {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs
        """

rule covariates:
    input:
        fam="{dir}/HVOL.DCM.FD.{type}.sigSNPs.fam",
        cov=expand("{dir}/20160412_All_BRU_format.txt",
            dir=config['phenodir']
            )
    output:
        "{dir}/HVOL.DCM.FD.{type}.covariates.txt"
    shell:
        """
        awk 'NR==FNR {{if ($3=="F") $3=2; if ($3=="M") $3=1; a[$1]=$3" "$5; \
            next}} $2 in a {{print $1" "$2" "a[$2]}}' \
            {input.cov} {input.fam} > {output}
        """

rule missing:
    input:
        both_bed="{dir}/HVOL.DCM.FD.{type}.sigSNPs.bed",
        both_bim="{dir}/HVOL.DCM.FD.{type}.sigSNPs.bim",
        both_fam="{dir}/HVOL.DCM.FD.{type}.sigSNPs.fam",
        related="{dir}/HVOL.DCM.related.IDs"
    output:
        missing="{dir}/HVOL.DCM.FD.{type}.sigSNPs.missing",
        ids="{dir}/HVOL.DCM.FD.{type}.sigSNPs.missing.ids",
        both_bed="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bed",
        both_bim="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bim",
        both_fam="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.fam",
    shell:
        """
        plink --bfile {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs \
            --test-missing \
            --allow-extra-chr \
            --allow-no-sex \
            --pfilter 0.01 \
            --out {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs
        tr -s ' ' < {output.missing} | sed 's/ //' | cut -d " " -f 2 | \
            tail -n +2 > {output.ids}
        plink --bfile {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs \
            --allow-extra-chr \
            --allow-no-sex \
            --remove {input.related} \
            --exclude {output.ids} \
            --make-bed \
            --out {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs.clean
        """

rule association:
    input:
        both_bed="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bed",
        both_bim="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bim",
        both_fam="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.fam",
        cov="{dir}/HVOL.DCM.FD.{type}.covariates.txt"
    output:
        logist="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic",
        perm="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic.perm",
        combined="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic.all"
    shell:
        """
        plink --bfile {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs.clean \
            --logistic hide-covar perm \
            --covar {input.cov} \
            --allow-extra-chr \
            --allow-no-sex \
            --out {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs.clean
        tr -s ' ' < {output.perm} | sed 's/ //' | tr ' ' '\t' | cut -f 3,4 | \
            paste {output.logist} -  > {output.combined}
        """

rule casecontrol:
    input:
        both_bed="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bed",
        both_bim="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bim",
        both_fam="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.fam",
        cov="{dir}/HVOL.DCM.FD.{type}.covariates.txt"
    output:
        cc="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.model",
        perm="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.model.best.perm",
        combined="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.model.all"
    shell:
        """
        plink --bfile {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs.clean \
            --model mperm=50000 \
            --allow-extra-chr \
            --allow-no-sex \
            --ci 0.95 \
            --out {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs.clean
        awk 'FNR==NR {{a[$2]=$3;next}} $2 in a {{print $0,a[$2]}}' \
            {output.perm} {output.cc} > {output.combined}
        """


rule significant:
    input:
        cc="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.model.all",
        genotyped=expand('{dir}/gencall.combined.bim',
            dir=config['genodir'])
    params:
        sigloci=17
    output:
        sig="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.model.perm.sig",
        geno="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.model.perm.sig.genotyped"
    shell:
        """
        awk 'FNR==NR && ($11*{params.sigloci}) < 0.05' {input.cc} > {output.sig}
        awk 'FNR==NR {{a[$2]=$0; next}} $2 in a {{print a[$2]}}' {output.sig} \
            {input.genotyped} > {output.geno}
        """
#### example cluster call #####
# snakemake -s impute.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

