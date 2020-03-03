# how to run from ./:
# snakemake -s mr.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} \
# -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_mr.yaml"

rule all:
    input:
        expand("{ukb}/MR/{pheno}/MR_hermes_dcm.pdf",
            pheno=config['discovery'],
            ukb=config["ukbdir"]),

rule mr:
    input:
        instruments="{dir}/gwas/{pheno}/Pseudomultitrait_slices_sig5e08_ldFiltered.txt",
        bimbam="{dir}/gwas/{pheno}/Pseudomultitrait_slices_sig5e08_genotypes.dosage.gz",
        dcm=expand("{dir}/HVOL.DCM.FD.Pseudomultitrait_Slices_sig5e08.sigSNPs.clean.assoc.logistic.all",
            dir=config["dcmdir"]),
        hermes=expand("{dir}/lookup_results.tsv",
            dir=config['hermesdir'])
    output:
        "{dir}/MR/{pheno}/MR_hermes_dcm.pdf"
    shell:
        "Rscript association/mendelian-randomisation.R \
            --gwasdir {wildcards.dir}/gwas/{wildcards.pheno} \
            --mrdir {wildcards.dir}/MR/{wildcards.pheno} \
            --dcm {input.dcm} \
            --hermes {input.hermes} \
            --showProgress "

