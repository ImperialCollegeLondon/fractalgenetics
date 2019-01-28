# how to run fro ./:
# snakemake -s control_plink.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_control.yaml"

rule all:
    input:
        expand("{dir}/control/plink_HVOL.{name}.combined.{suffix}_batch.assoc.linear",
            dir=config["controldir"],
            name=config["name"],
            suffix=config["suffix"]),
        expand("{dir}/control/HVOL.{name}.combined.{suffix}.batchPhenotypes.fam",
            dir=config["controldir"],
            name=config["name"],
            suffix=config["suffix"])

rule prepare_files:
    input:
        bgen="{dir}/imputation/combined/genotypes/gencall.combined.clean.related.chr1.sample",
        HVOL="{dir}/QC/combined/European.HVOL.gencall.combined.txt"
    params:
        batch1="sanger12",
        batch2="singapore12",
        batch3="singapore3"
    output:
        "{dir}/control/HVOL.{name}.combined.{suffix}.batchPhenotypes.fam"
    shell:
        "Rscript scripts/prepare_files.R \
            --directory={wildcards.dir} \
            --name={wildcards.name} \
            --batch1={params.batch1} \
            --batch2={params.batch2} \
            --batch3={params.batch3} \
            --suffix={wildcards.suffix} \
            --bgen={input.bgen} \
            --HVOl={input.HVOL} \
            --showProgress"

rule association_genotyped:
    input:
        bed=expand("{geno}/{{name}}.combined.{{suffix}}.bed",
            geno=config["genodir"]),
        bim=expand("{geno}/{{name}}.combined.{{suffix}}.bim",
            geno=config["genodir"]),
        fam="{dir}/control/HVOL.{name}.combined.{suffix}.batchPhenotypes.fam"
    output:
        "{dir}/control/plink_HVOL.{name}.combined.{suffix}_batch.assoc.linear"
    shell:
        """
        plink --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --prune --mac 20 \
            --out {wildcards.dir}/plink_HVOL.{wildcards.name}.combined.{wildcards.suffix}_batch
        """

rule association_results:
    input:
        "{dir}/plink_HVOL.{name}_batch.assoc.linear"
    output:
        "{dir}/plink_HVOL.{name}_lm_qqplot.pdf",
        "{dir}/plink_HVOL.{name}_lm_manhattanplot.pdf"
    shell:
        "Rscript scripts/association_plink.R \
            --name HVOL.{wildcards.name}_batch \
            --directory {wildcards.dir} \
            --showProgress"

