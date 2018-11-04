# how to run fro ./:
# snakemake -s control_plink.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_control.yaml"

rule all:
    input:
        expand("{dir}/plink_{name}_batch.qassoc",
            dir=config["controldir"],
            name=config["name"])

rule association_genotyped:
    input:
        bed=expand("{geno}/{{name}}.bed", geno=config["genodir"]),
        bim=expand("{geno}/{{name}}.bim", geno=config["genodir"]),
        fam="{dir}/HVOL.{name}.batchPhenotypes.fam"
    output:
        "{dir}/plink_HVOL.{name}_batch.assoc.linear"
    shell:
        """
        plink --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --prune --mac 20 \
            --out {wildcards.dir}/plink_HVOL.{wildcards.name}_batch
        """

rule associationResults:
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
