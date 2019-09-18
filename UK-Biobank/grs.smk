# how to run from ./:
# snakemake -s grs.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} \
# -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_grs.yaml"

def assign_linker(wildcards):
    if wildcards.type == "fd":
        return "{}/rawdata/ukb18545_imp_chr1_v3_s487378.sample".format(wildcards.dir)
    else:
        return "{}/rawdata/ukb40616_imp_chr1_v3_s487317.sample".format(wildcards.dir)

rule all:
    input:
        expand("{ukb}/maf0.001/{pheno}/ukb_imp_genome_v3_maf0.001.rsid",
            ukb=config["ukbdir"],
            pheno=config["pheno"]),
        expand("{ukb}/grs/genotypes/genotypes_{type}.{fileformat}",
            fileformat=["bed", "dosage.gz"],
            ukb=config["ukbdir"],
            type=['fd', 'hf', 'icm', 'sz_nicm', 'nicm', 'cad']),

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

rule format_sampleids:
    input:
        heartsamples="{dir}/heart_phenotypes/{type}_eid.csv",
    output:
        sampleids="{dir}/grs/genotypes/samples_{type}_qctool.IDs",
    shell:
        """
        tr '\n' ' ' < {input.heartsamples} > {output.sampleids}
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
            -og {output.plink} -s {input.samples}
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
            -og {output.bimbam} -s {input.samples}
        """

#rule extract_genotypes_hf:
#    input:
#        geno=expand("{geno}/ukb_imp_chr{chr}_v3.bgen",
#            chr=range(1,23),
#            geno=config["genodir"]),
#        samples="{dir}/Heart_failure_phenotypes/all_eid.sample",
#        hfsamples="{dir}/Heart_failure_phenotypes/{type}_eid.csv"
#    output:
#        bimbam="{dir}/grs/genotypes_{type}.dosage.gz",
#        bed="{dir}/grs/genotypes_{type}.bed",
#        sampleids="{dir}/grs/samples_{type}_qctool.IDs",
#    params:
#        genodir=config["genodir"]
#    shell:
#       """
#        cut -d "," -f 1 {input.hfsamples} | tr '\n' ' ' > {output.sampleids}
#        qctool -g {params.genodir}/ukb_imp_chr#_v3.bgen \
#            -incl-samples {output.sampleids} \
#            -ofiletype dosage \
#            -og {output.bimbam} -s {input.samples}
#        qctool -g {params.genodir}/ukb_imp_chr#_v3.bgen \
#            -incl-samples {output.sampleids} \
#            -ofiletype binary_ped \
#            -og {output.bed} -s {input.samples}
#        """

rule grs:
    input:
        samples="{dir}/grs/genotypes/genotypes_{type}.bed"
    output:
        results="{dir}/grs/{type}_grs.csv"
    shell:
        "Rscript PRSice.R --dir {wildcards.dir} \
            --prsice PRSice \
            --base TOY_BASE_GWAS.assoc \
            --target TOY_TARGET_DATA \
            --thread 1 \
            --stat BETA \
            --beta \
            --binary-target F"
