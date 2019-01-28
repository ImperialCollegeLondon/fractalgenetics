# how to run fro ./:
# snakemake -s control.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_control.yaml"

rule all:
    input:
        expand("{dir}/bgenie_HVOL.{name}_batch_lm_st_chr{chr}.gz",
            dir=config["controldir"],
            name=config["name"],
            chr=range(1,23)),

rule association_impute:
    input:
        geno=expand("{geno}/{{name}}.chr{{chr}}.qc.bgen",
            geno=config["imputedir"]),
        pheno="{dir}/Batch_phenotypes_HVOL.{name}_bgenie.txt"
    params:
        n=config["n"]
    output:
        "{dir}/bgenie_HVOL.{name}_batch_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --scale_phenotypes \
            --pvals --exc_missing_inds \
            --thread {params.n} \
            --out {wildcards.dir}/bgenie_HVOL.{wildcards.name}_batch_lm_st_chr{wildcards.chr}"

rule results_impute:
    input:
        expand("{{dir}}/bgenie_HVOL.{{name}}_batch_lm_st_chr{chr}.gz",
            chr=range(1,23))
    output:
        "{dir}/bgenie_HVOL.{name}_lm_qqplot.pdf",
        "{dir}/bgenie_HVOL.{name}_lm_manhattanplot.pdf"
    shell:
        "Rscript scripts/association_impute.R \
            --name=HVOL.{wildcards.name}_batch \
            --directory={wildcards.dir} \
            --showProgress "
