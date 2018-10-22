# how to run fro ./:
# snakemake -s genotypeQC.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_genotypeQC.yaml"

rule all:
    input:
        expand("{dir}/{alg}/{call}.{alg}.{suffix}",
            dir=config["dir"],
            call='gencall',
            alg='combined',
            suffix=['clean.related.bim', 'clean.related.bed', 'clean.related.fam']),
        expand("{dir}/{alg}/{call}.{alg}.{suffix}",
            dir=config["dir"],
            call='gencall',
            alg='combined',
            suffix=['HapMapIII.eigenvec', 'HapMapIII.eigenval']),
        expand("{dir}/{alg}/European.HVOL.{call}.{alg}.txt",
            call='gencall',
            alg='combined',
            dir=config["dir"])

rule genotypeQC:
    input:
        "{qcdir}/{alg}/{call}.{alg}.bim",
        "{qcdir}/{alg}/{call}.{alg}.fam",
        "{qcdir}/{alg}/{call}.{alg}.bed",
    params:
        highld=config['highld'],
    output:
        "{qcdir}/{alg}/{call}.{alg}.sexcheck",
        "{qcdir}/{alg}/{call}.{alg}.lmiss",
        "{qcdir}/{alg}/{call}.{alg}.het",
        "{qcdir}/{alg}/{call}.{alg}.genome",
        "{qcdir}/{alg}/{call}.{alg}.prune.in",
        "{qcdir}/{alg}/{call}.{alg}.pruned.bim"
    wildcard_constraints:
        alg="\w*"
    shell:
        """
        bash scripts/genotypeQC.sh {wildcards.call}.{wildcards.alg} \
            {wildcards.qcdir}/{wildcards.alg} {params.highld}
        """

rule ancestryCheck:
    input:
        "{qcdir}/{alg}/{call}.{alg}.bim",
        "{qcdir}/{alg}/{call}.{alg}.fam",
        "{qcdir}/{alg}/{call}.{alg}.bed",
        "{qcdir}/{alg}/{call}.{alg}.prune.in",
        "{qcdir}/{alg}/{call}.{alg}.pruned.bim"
    params:
        highld=config['highld'],
        hapmap=config['hapmap'],
        hapmapdir=config['hapmapdir']
    output:
        "{qcdir}/{alg}/{call}.{alg}.{ref}.eigenvec",
        "{qcdir}/{alg}/{call}.{alg}.{ref}.eigenval"
    wildcard_constraints:
        alg="\w*"
    shell:
        """
        bash scripts/ancestryCheck.sh {wildcards.call}.{wildcards.alg} \
            {wildcards.qcdir}/{wildcards.alg} {params.highld} {wildcards.ref} \
            {params.hapmap} {params.hapmapdir}
        """

rule filtering_and_plots:
    input:
        "{qcdir}/{alg}/{call}.{alg}.bim",
        "{qcdir}/{alg}/{call}.{alg}.fam",
        "{qcdir}/{alg}/{call}.{alg}.bed",
        "{qcdir}/{alg}/{call}.{alg}.sexcheck",
        "{qcdir}/{alg}/{call}.{alg}.lmiss",
        "{qcdir}/{alg}/{call}.{alg}.het",
        "{qcdir}/{alg}/{call}.{alg}.genome",
        "{qcdir}/{alg}/{call}.{alg}.{ref}.eigenvec",
        "{qcdir}/{alg}/{call}.{alg}.{ref}.eigenval"
    output:
        "{qcdir}/{alg}/{call}.{alg}.clean.related.bim",
        "{qcdir}/{alg}/{call}.{alg}.clean.related.fam",
        "{qcdir}/{alg}/{call}.{alg}.clean.related.bed"
    params:
        maleTh=config['maleTh'],
        femaleTh=config['femaleTh'],
        imissTh=config['imissTh'],
        hetTh=config['hetTh'],
        highIBDTh=config['highIBDTh'],
        lmissTh=config['lmissTh'],
        hweTh=config['hweTh'],
        macTh=config['macTh'],
        externalSex=config['externalSex'],
        externalSexID=config['externalSexID'],
        externalSexSex=config['externalSexSex'],
        refSamples=config['refSamples'],
        refColors=config['refColors'],
        path2plink=config['path2plink'],
        ref=config['reference']
    shell:
        """
        Rscript QC/genotypeQC.R \
            --name={wildcards.call}.{wildcards.alg} \
            --directory={wildcards.qcdir} \
            --check_sex \
            --fixMixup \
            --maleTh={params.maleTh} \
            --femaleTh={params.femaleTh} \
            --externalSex={params.externalSex} \
            --externalSexSex={params.externalSexSex} \
            --externalSexID={params.externalSexID} \
            --check_het_and_miss \
            --imissTh={params.imissTh} \
            --hetTh={params.hetTh}  \
            --check_relatedness \
            --relatednessTh={params.highIBDTh} \
            --check_ancestry \
            --prefixMerge={wildcards.call}.{wildcards.alg}.{params.ref} \
            --refSamples={params.refSamples} \
            --refColors={params.refColors} \
            --lmissTh={params.lmissTh} \
            --hweTh={params.hweTh} \
            --macTh={params.macTh} \
            --plink={params.path2plink} \
            --plot \
            --showProgress
        """

rule HVOLsamples:
    input:
        "{dir}/{alg}/{call}.{alg}.clean.related.fam",
    params:
        diagnosis=config['diagnosis']
    output:
        HVOL="{dir}/{alg}/European.HVOL.{call}.{alg}.txt",
    shell:
        """
        awk 'FNR==NR {{a[$1]; next}} ($2 in a && $6 == "HVOL")' {input} \
            {params.diagnosis} > {output.HVOL}
        """
