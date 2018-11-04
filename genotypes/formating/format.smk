# how to run fro ./:
# snakemake -s format.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_format.yaml"

rule all:
    input:
        expand("{dir}/{alg}/{call}.{alg}.{suffix}",
            dir=config["dir"],
            call='gencall',
            alg=['sanger12', 'singapore12', 'singapore3'],
            suffix=['bim', 'bed', 'fam']),
        expand("{dir}/combined/{suffix}",
            dir=config["dir"],
            suffix=['overlap_genotypeProbes.pdf', 'same_genotype_probes.txt']),
        expand("{dir}/{alg}/{call}.{alg}.{suffix}",
            dir=config["dir"],
            alg=['combined'],
            call='gencall',
            suffix=['bim', 'bed', 'fam']),
        #expand("{dir}/European.HVOL.sanger12.singapore123.txt",
         #   dir=config["dir"])

rule format_raw_genotypes:
    input:
        preQCfail="{qcdir}/{call}.{alg}.preQCfail.IDs",
    params:
        center=lambda wildcards: config['center'][wildcards.alg],
        rawdir=lambda wildcards: config['rawdir'][wildcards.alg],
        rawdata=lambda wildcards: config['rawdata'][wildcards.alg],
        sample=lambda wildcards: config['sample'][wildcards.alg]
    output:
        "{qcdir}/{call}.{alg}.raw.bim",
        "{qcdir}/{call}.{alg}.raw.fam",
        "{qcdir}/{call}.{alg}.raw.bed"
    shell:
        """
        bash scripts/formatRawGenotypes.sh {wildcards.call}.{wildcards.alg} \
            {wildcards.qcdir} {params.rawdir} {params.center} \
            {params.rawdata} {params.sample}
        """

rule match_to_reference:
    input:
        "{qcdir}/{call}.{alg}.raw.bim",
        "{qcdir}/{call}.{alg}.raw.fam",
        "{qcdir}/{call}.{alg}.raw.bed"
    params:
        UK10K1KGdir=config['UK10K1KGdir'],
        ensemblNotInUK10K1KG=config['ensemblNotInUK10K1KG']
    output:
        "{qcdir}/{call}.{alg}.bim",
        "{qcdir}/{call}.{alg}.fam",
        "{qcdir}/{call}.{alg}.bed"
    wildcard_constraints:
        call="\w+",
        alg="[\w\d]+"
    shell:
        """
        bash scripts/match2reference.sh {wildcards.call}.{wildcards.alg} \
            {wildcards.qcdir} {params.UK10K1KGdir} {params.ensemblNotInUK10K1KG}
        """

rule check_probes:
    input:
        sanger12=config['probes']['sanger12'],
        singapore12=config['probes']['singapore12'],
        singapore3=config['probes']['singapore3']
    output:
        "{dir}/overlap_genotypeProbes.pdf",
        "{dir}/same_genotype_probes.txt"
    shell:
        """
        Rscript scripts/compare_probeID_across_batches.R \
            --sanger12={input.sanger12} \
            --singapore12={input.singapore12} \
            --singapore3={input.singapore3} \
            --directory={wildcards.dir}
        """

rule check_chips:
    input:
        expand("{{dir}}/{alg}/{{call}}.{alg}.bim",
            alg=['singapore12', 'singapore3', 'sanger12']),
        expand("{{dir}}/{alg}/{{call}}.{alg}.bed",
            alg=['singapore12', 'singapore3', 'sanger12']),
        expand("{{dir}}/{alg}/{{call}}.{alg}.fam",
            alg=['singapore12', 'singapore3', 'sanger12']),
        commonprobes="{dir}/combined/same_genotype_probes.txt"
    params:
        lmissTh=config['lmissTh']
    output:
        "{dir}/combined/{call}.combined.bim",
        "{dir}/combined/{call}.combined.bed",
        "{dir}/combined/{call}.combined.fam"
    shell:
        """
        bash scripts/create_combined_dataset.sh \
            {wildcards.dir} \
            'sanger12 singapore12 singapore3' \
            {wildcards.call} \
            {params.lmissTh} \
            {input.commonprobes} \
        """

