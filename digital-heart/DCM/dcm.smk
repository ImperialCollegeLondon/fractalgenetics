from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

###### load config file #####
configfile: "config/config_dcm.yaml"


##### target rules #####
rule all:
    input:
        expand("{dir}/gwas/DCM/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic.all.{suffix}",
            suffix=['sig', 'sig.genotyped'],
            type=['Pseudomultitrait_Slices_sig5e08'],
            dir=config["dir"]),
        #expand("{ukb}/DCM/FDAlongHeart_DCM_UKB_all_slices{interpolate}.pdf",
        #    ukb=config['ukbdir'],
        #    interpolate=config["interpolate"])

rule process_FD:
    input:
        pheno="{dir}/phenotype/FD/20191002_DCM_FD_all.csv",
        covs="{dir}/phenotype/2Dphenotype/20160412_All_BRU_format.txt",
        samples="{dir}/genotype/imputation/combined/genotypes/gencall.combined.clean.related.chr1.sample",
        europeans="{dir}/genotype/QC/combined/DCM.gencall.combined.clean.related.fam",
    output:
        "{dir}/phenotype/FD/FD_summary_EUnorel.csv",
        "{dir}/phenotype/FD/FD_slices_EUnorel.csv",
    params:
        interpolate=config['interpolate'],
        genodir=config['genodir'],
        plink=config['plink']
    shell:
        "Rscript 'preparePheno.r' \
            --outdir {wildcards.dir}/phenotype/FD \
            --cov {input.covs} \
            --pheno {input.pheno} \
            --samples {input.samples} \
            --europeans {input.europeans} \
            --interpolate {params.interpolate}"

rule fd_ukb_dcm:
    input:
        ukb=expand("{{ukb}}/phenotypes/{pheno}/FD_slices_EUnorel.csv",
            pheno=config['discovery']),
        dcm=expand("{dh}/FD/FD_slices_EUnorel.csv",
            dh=config['phenodir'])
    output:
        "{ukb}/DCM/FDAlongHeart_DCM_UKB_all_slices{interpolate}.pdf"
    shell:
        "Rscript 'compare-FD-dcm-ukb.R' \
            --dir {wildcards.ukb}/DCM \
            --dcm {input.dcm} \
            --interpolate {wildcards.interpolate} \
            --ukb {input.ukb}"


rule extract:
    input:
        bgen=expand("{dir}/formated/{name}.genome.qc.bgen",
            name=config['name'],
            dir=config['imputedir']),
        sample=expand("{dir}/genotypes/{name}.chr1.sample",
            name=config['name'],
            dir=config['imputedir']),
        rsids=expand('{dir}/gwas/180628_fractal_dimension/Pseudomultitrait_slices_sig5e08_qctool.IDs',
            dir=config['ukbdir']),
        update=expand('{dir}/gwas/180628_fractal_dimension/Pseudomultitrait_slices_sig5e08_qctool.toUpdate',
            dir=config['ukbdir']),
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
        cov=expand("{dir}/2Dphenotype/20160412_All_BRU_format.txt",
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

rule freq:
    input:
        both_bed="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bed",
        both_bim="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bim",
        both_fam="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.fam",
    output:
        freq="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.frq"
    shell:
        """
        plink --bfile {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs.clean \
            --freq \
            --allow-extra-chr \
            --allow-no-sex \
            --out {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs.clean
        """

rule association:
    input:
        both_bed="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bed",
        both_bim="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.bim",
        both_fam="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.fam",
        cov="{dir}/HVOL.DCM.FD.{type}.covariates.txt",
        freq="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.frq"
    output:
        logist="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic",
        perm="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic.perm",
        combined="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic.all"
    shell:
        """
        plink --bfile {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs.clean \
            --logistic beta hide-covar perm \
            --covar {input.cov} \
            --allow-extra-chr \
            --allow-no-sex \
            --ci 0.95 \
            --seed 45908 \
            --out {wildcards.dir}/HVOL.DCM.FD.{wildcards.type}.sigSNPs.clean
        tr -s ' ' < {output.perm} | sed 's/ //' | cut -d " " -f 3,4 \
            > {output.perm}.tmp
        tr -s ' ' < {output.logist} | sed 's/ //' | \
            cut -d " " --complement -f 4  > {output.logist}.tmp
        tr -s ' ' < {input.freq} | sed 's/ //' | \
            cut -d " " -f 3,4,5  > {input.freq}.tmp
        paste -d " "  {output.logist}.tmp {output.perm}.tmp {input.freq}.tmp \
            > {output.combined}
        rm {output.perm}.tmp {output.logist}.tmp {input.freq}.tmp
        """

rule significant:
    input:
        combined="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic.all",
        genotyped=expand('{dir}/gencall.combined.bim',
            dir=config['genodir'])
    params:
        sigloci=16
    output:
        sig="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic.all.sig",
        geno="{dir}/HVOL.DCM.FD.{type}.sigSNPs.clean.assoc.logistic.all.sig.genotyped"
    shell:
        """
        awk 'FNR==NR && $12 < 0.05' {input.combined} > {output.sig}
        awk 'FNR==NR {{a[$2]=$0; next}} $2 in a {{print a[$2]}}' {output.sig} \
            {input.genotyped} > {output.geno}
        """


#### example cluster call #####
# snakemake -s impute.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

