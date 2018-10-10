# how to run fro ./:
# snakemake -s genotypeQC.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_genotypeQC.yaml"

rule all:
    input:
        expand("{dir}/{alg}/{call}.{alg}.{suffix}",
            dir=config["dir"],
            call='gencall',
            alg=['sanger12', 'singapore12', 'singapore3'],
            suffix=['failsex', 'lmiss', 'genome', 'het', 'HapMapIII.eigenvec'])
            #suffix=['bim', 'bed', 'fam'])

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
        bash format/formatRawGenotypes.sh {wildcards.call}.{wildcards.alg} \
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
    shell:
        """
        bash format/match2reference.sh {wildcards.call}.{wildcards.alg} \
            {wildcards.qcdir} {params.UK10K1KGdir} {params.ensemblNotInUK10K1KG}
        """


rule genotypeQC:
    input:
        "{qcdir}/{alg}.bim",
        "{qcdir}/{alg}.fam",
        "{qcdir}/{alg}.bed"
    params:
        highld=config['highld'],
    output:
        "{qcdir}/{alg}.failsex",
        "{qcdir}/{alg}.lmiss",
        "{qcdir}/{alg}.het",
        "{qcdir}/{alg}.genome"
    shell:
        """
        bash QC/genotypeQC.sh {wildcards.alg} {wildcards.qcdir} \
            {params.highld}
        """

rule ancestryCheck:
    input:
        "{qcdir}/{alg}.bim",
        "{qcdir}/{alg}.fam",
        "{qcdir}/{alg}.bed"
    params:
        highld=config['highld'],
        hapmap=config['hapmap'],
        hapmapdir=config['hapmapdir']
    output:
        "{qcdir}/{alg}.{ref}.eigenvec",
        "{qcdir}/{alg}.{ref}.eigenval"
    shell:
        """
        bash QC/ancestryCheck.sh {wildcards.alg} {wildcards.qcdir} \
            {params.highld} {wildcards.ref} {params.hapmap} {params.hapmapdir}
        """

rule filtering_and_plots:
    input:
        sample="",
    output:
        "{qcdir}/"
    shell:
        """
        Rscript QC/genotypeQC.R --vanilla \
            --default-packages=R.utils \
            --alg={wildcards.alg} \
            --qcdir={wildcars.qcdir} \
            --maleTh={params.maleTh} \
            --femaleTh={params.femaleTh} \
            --imissTh={params.imissTh} \
            --hetTh={params.hetTh}  \
            --highIBDTh={params.highIBDTh} \
            --lmissTh={params.lmissTh} \
            --hweTh={params.hweTh} \
            --mafTh={params.mafTh} \
            --sample={input.sample}
        """



