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
            suffix=['failsex','clean.related.bim', 'clean.related.bed', 'clean.related.fam']),
        expand("{dir}/{alg}/{call}.{alg}.{suffix}",
            dir=config["dir"],
            call='gencall',
            alg=['sanger12', 'singapore12', 'singapore3'],
            suffix=['HapMapIII.eigenvec', 'HapMapIII.eigenval']),
        expand("{dir}/{alg}/{file}",
            dir=config["dir"],
            alg=['combined'],
            file=['gencall.combined.clean.related.bim',
                    'gencall.combined.clean.related.bed',
                    'gencall.combined.clean.related.fam']),
        expand("{dir}/European.HVOL.sanger12.singapore123.txt",
            dir=config["dir"])

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
    wildcard_constraints:
        call="\w+",
        alg="[\w\d]+"
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
        "{qcdir}/{alg}.bim",
        "{qcdir}/{alg}.fam",
        "{qcdir}/{alg}.bed",
        "{qcdir}/{alg}.failsex",
        "{qcdir}/{alg}.lmiss",
        "{qcdir}/{alg}.het",
        "{qcdir}/{alg}.genome"
    output:
        "{qcdir}/{alg}.clean.related.bim",
        "{qcdir}/{alg}.clean.related.fam",
        "{qcdir}/{alg}.clean.related.bed"
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
            --name={wildcards.alg} \
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
            --prefixMerge={wildcards.alg}.{params.ref} \
            --refSamples={params.refSamples} \
            --refColors={params.refColors} \
            --lmissTh={params.lmissTh} \
            --hweTh={params.hweTh} \
            --macTh={params.macTh} \
            --plink={params.path2plink} \
            --plot \
            --showProgress
        """

rule europeanSamples:
    input:
        sanger12="{dir}/sanger12/gencall.sanger12.clean.related.fam",
        singapore12="{dir}/singapore12/gencall.singapore12.clean.related.fam",
        singapore3="{dir}/singapore3/gencall.singapore3.clean.related.fam",
    params:
        diagnosis=config['diagnosis']
    output:
        all="{dir}/European.sanger12.singapore123.txt",
        HVOL="{dir}/European.HVOL.sanger12.singapore123.txt"
    shell:
        """
        cat {input.sanger12} {input.singapore12} {input.singapore3} > {output.all}
        awk 'FNR==NR {{a[$1]; next}} ($2 in a && $6 == "HVOL")' {output.all} \
            {params.diagnosis} > {output.HVOL}
        """

rule mergeDatasets:
    input:
        sanger12_bim="{dir}/sanger12/gencall.sanger12.clean.related.bim",
        sanger12_bed="{dir}/sanger12/gencall.sanger12.clean.related.bed",
        sanger12_fam="{dir}/sanger12/gencall.sanger12.clean.related.fam",
        singapore12_bim="{dir}/singapore12/gencall.singapore12.clean.related.bim",
        singapore12_bed="{dir}/singapore12/gencall.singapore12.clean.related.bed",
        singapore12_fam="{dir}/singapore12/gencall.singapore12.clean.related.fam",
        singapore3_bim="{dir}/singapore3/gencall.singapore3.clean.related.bim",
        singapore3_bed="{dir}/singapore3/gencall.singapore3.clean.related.bed",
        singapore3_fam="{dir}/singapore3/gencall.singapore3.clean.related.fam"
    output:
        all_bim="{dir}/combined/gencall.combined.clean.related.bim",
        all_bed="{dir}/combined/gencall.combined.clean.related.bed",
        all_fam="{dir}/combined/gencall.combined.clean.related.fam",
        mergelist="{dir}/combined/gencall.combined.mergelist"
    shell:
        """
        echo -e '{input.sanger12_bed} {input.sanger12_bim} {input.sanger12_fam}\n\
                 {input.singapore12_bed} {input.singapore12_bim} {input.singapore12_fam}\n\
                 {input.singapore3_bed} {input.singapore3_bim} {input.singapore3_fam}\n' \
                 > {output.mergelist}
        plink --merge-list {output.mergelist} \
            --make-bed \
            --out {wildcards.dir}/combined/gencall.combined.clean.related
        """
