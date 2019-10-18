# how to run fro ./:
# snakemake -s phenotypes.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete

configfile: "config/config_phenotypes.yaml"

rule all:
    input:
        expand("{dir}/phenotypes/FD_slices_bgenie.txt",
            dir=config["dir"]),

rule formatPhenotypes:
    input:
        pheno="{dir}/phenotype/FD/20181116_HVOLSegmentations_FD.csv",
        cov="{dir}/phenotype/2Dphenotype/20160705_GenScan.txt",
        samples="{dir}/genotype/imputation/combined/genotypes/gencall.combined.clean.related.chr1.sample",
        europeans="{dir}/genotype/QC/combined/European.HVOL.gencall.combined.IDs",
        pcs="{dir}/genotype/QC/combined/HVOL.gencall.combined.clean.related.eigenvec",
        path2plink="/homes/hannah/bin/plink",
        genodir="{dir}/genotype/QC/combined"
    output:
        "{dir}/phenotypes/FD/FD_slices_bgenie.txt",
    shell:
        "Rscript 'phenotypes/preparePheno.r' -gd={input.genodir} \
            -o={wildcards.dir}/phenotype/FD \
            -c={input.cov} \
            -p={input.pheno} \
            -s={input.samples} \
            -e={input.europeans} \
            -pcs={input.pcs} \
            -path2plink={input.path2plink}"

