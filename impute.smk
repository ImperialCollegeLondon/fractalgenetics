# how to run fro ./:
# snakemake -s impute.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

import numpy as np
from functools import reduce
from operator import add

# rule input functions
def createStringPlink(wildcards):
    if wildcards.chr == "X_PAR1":
        plink_str="--chr X --from-bp 60001 --to-bp 2699520"
    elif wildcards.chr == "X_PAR2":
        plink_str="--chr X --from-bp 154931044 --to-bp 999999999"
    elif wildcards.chr == "X":
        plink_str="--chr X --from-bp 2699521 --to-bp 154931043"
    else:
        plink_str="--chr {}".format(wildcards.chr)
    return plink_str

def createStringImpute(wildcards):
    if wildcards.chr == "X_PAR1":
        impute_str="-chrX -Xpar"
    elif wildcards.chr == "X_PAR2":
        impute_str="-chrX -Xpar"
    elif wildcards.chr == "X":
        impute_str="-chrX"
    else:
        impute_str=""
    return impute_str

def chunkSettings(wildcards):
    dir="/homes/hannah/data/hmeyer/UK10K1000Genomes"
    with open ("{}/chunkBoundariesChr{}.txt".format(dir, wildcards.chr)) as file:
        content = file.readlines()
        content = [x.strip() for x in content]
        for line in content:
            parts = line.split()
            if parts[1] == str(wildcards.chunk):
                if len(parts) > 2:
                    return "{} {}".format(parts[2], parts[3])
        return False

def chunks2concatenate(wildcards):
    return expand("{{dir}}/imputed/chr{{chr}}/{{name}}.chr{{chr}}.{chunk}.gen",
                chunk=chunksPerChromosome(chromosome=[wildcards.chr],
                    dir=config['referencedir'],
                    noSnpInInterval=noSnpInInterval)['chunks'])

def chunks2check(wildcards):
    warnings=expand("{{dir}}/imputed/chr{{chr}}/{{name}}.chr{{chr}}.{chunk}.gen_warnings",
                chunk=chunksPerChromosome(chromosome=[wildcards.chr],
                    dir=config['referencedir'],
                    noSnpInInterval=noSnpInInterval)['chunks'])
    return warnings

# global parameters functions
def chunksPerChromosome(chromosome, dir, noSnpInInterval=None):
    if noSnpInInterval is not None:
        exclude = dict(item.split(".") for item in noSnpInInterval)
    else:
        exclude = None

    chunks_list = []
    chr_list = []

    for c in chromosome:
        with open ("{}/chunkBoundariesChr{}.txt".format(dir, c)) as file:
            content = file.readlines()
            content = [x.strip() for x in content]
            for line in content:
                parts = line.split()
                if len(parts) > 2:
                    chr = parts[0]
                    chunk = parts[1]
                    if exclude is not None:
                        if not(chr in exclude and exclude[chr] == chunk):
                            chunks_list.append(chunk)
                            chr_list.append(chr)
                    else:
                        chunks_list.append(chunk)
                        chr_list.append(chr)
    return {'chunks': chunks_list, 'chr': chr_list}


# config file
configfile: "config/config_impute.yaml"

# global parameters
## after first run manually extraacted from *no_chunks by combineCheckChunks
noSnpInInterval = ['9.22', '1.48', '21.1']
CHR = list(range(1,23)) + ["X"]
AUTOSOMES = list(range(1,23))
perChr = chunksPerChromosome(chromosome=CHR, dir=config['referencedir'],
    noSnpInInterval=noSnpInInterval)
CHRlist = perChr['chr']
CHUNKlist = perChr['chunks']

## QC columns in qctool -snp-stats output
qctool_missing_column = 19
qctool_info_column = 17
qctool_maf_column = 14
qctool_hwe_column = 8
qctool_rsid_column = 2

rule all:
    input:
        # rule splitChromosomes
        expand("{dir}/unphased/{name}.chr{chr}.{suffix}",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23),
            suffix=['bim', 'bed', 'fam']),
        # rule phasing
        expand("{dir}/phased/{name}.chr{chr}.{suffix}",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23),
            suffix=['hap.gz', 'sample']),
        # rule imputation
        expand(expand("{dir}/imputed/chr{{chr}}/{name}.chr{{chr}}.{{chunk}}.{suffix}",
                dir=config["dir"],
                name='gencall.combined.clean.related',
                suffix=['gen']),
            zip, chunk=CHUNKlist, chr=CHRlist),
        # rule checkChunks
        expand("{dir}/imputed/chr{chr}/{name}.chr{chr}.warnings_overview",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23)),
        # rule combineCheckChunks
        expand("{dir}/imputed/{name}.warnings_overview",
            dir=config["dir"],
            name='gencall.combined.clean.related'),
        # rule concatenateChunks
        expand("{dir}/genotypes/{name}.chr{chr}.{suffix}",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23),
            suffix=['sample', 'gen']),
        # rule snpQC
        expand("{dir}/genotypes/{name}.chr{chr}.qc.gen.gz",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23)),
        # rule combine_and_convert
        expand(expand("{dir}/{{format}}/{name}.genome.qc.{{suffix}}",
                dir=config["dir"],
                name='gencall.combined.clean.related'),
            zip, format=['bgen', 'bimbam'], suffix=['bgen', 'dosage.gz'])


rule splitChromosomes:
    input:
        bim=expand("{qcdir}/combined/{{name}}.bim",
            qcdir=config['qcdir']),
        bed=expand("{qcdir}/combined/{{name}}.bed",
            qcdir=config['qcdir']),
        fam=expand("{qcdir}/combined/{{name}}.fam",
            qcdir=config['qcdir'])
    params:
        stringPlink=lambda wildcards: createStringPlink(wildcards)
    output:
        bim="{dir}/unphased/{name}.chr{chr}.bim",
        bed="{dir}/unphased/{name}.chr{chr}.bed",
        fam="{dir}/unphased/{name}.chr{chr}.fam",
    shell:
        """
        plink --bed {input.bed} \
              --bim {input.bim} \
              --fam {input.fam} \
              {params.stringPlink} \
              --make-bed \
              --out {wildcards.dir}/unphased/{wildcards.name}.chr{wildcards.chr}
        """

rule phasing:
    input:
        bim="{dir}/unphased/{name}.chr{chr}.bim",
        bed="{dir}/unphased/{name}.chr{chr}.bed",
        fam="{dir}/unphased/{name}.chr{chr}.fam"
    output:
        hap="{dir}/phased/{name}.chr{chr}.hap.gz",
        sample="{dir}/phased/{name}.chr{chr}.sample"
    params:
        phasethreads=config['threads'],
        window=config['window'],
        states=config['states'],
        map=lambda wildcards: "{}/genetic_map_chr{}_combined_b37.txt".format(
            config['referencedir'], wildcards.chr),
        stringPhasing=lambda wildcards: "--chrX" if wildcards.chr =="X" else ""
    shell:
        """
        shapeit --thread {params.phasethreads} \
            --window {params.window} \
            --states {params.states} \
            --effective-size 11418 \
            -B {wildcards.dir}/unphased/{wildcards.name}.chr{wildcards.chr} \
            --input-map {params.map} \
            --output-log {wildcards.dir}/phased/chr{wildcards.chr}.shapeit \
            --output-max {output.hap} {output.sample} {params.stringPhasing}
        """

rule imputation:
    input:
        hap="{dir}/phased/{name}.chr{chr}.hap.gz",
        sample="{dir}/phased/{name}.chr{chr}.sample"
    output:
        gen="{dir}/imputed/chr{chr}/{name}.chr{chr}.{chunk}.gen"
    params:
        buffer=config['buffer_size'],
        k_hap=config['k_hap'],
        states=config['states'],
        Ne=config['Ne'],
        ref=config['referencedir'],
        map=lambda wildcards: "{}/genetic_map_chr{}_combined_b37.txt".format(
            config['referencedir'], wildcards.chr),
        reflegend=lambda wildcards: "{}/chr{}.{:02d}.shapeit.legend.gz".format(
            config['referencedir'], wildcards.chr, int(wildcards.chunk)),
        refhap=lambda wildcards: "{}/chr{}.{:02d}.shapeit.hap.gz".format(
            config['referencedir'], wildcards.chr, int(wildcards.chunk)),
        boundaries=lambda wildcards: chunkSettings(wildcards),
        stringImpute=lambda wildcards: createStringImpute(wildcards),
    shell:
         """
         impute2 -allow_large_regions \
                    -m {params.map} \
                    -h {params.refhap} -l {params.reflegend} -use_prephased_g \
                    -known_haps_g {input.hap} \
                    -sample_g {input.sample} \
                    -k_hap {params.k_hap} \
                    -Ne {params.Ne} -buffer {params.buffer} -verbose \
                    -int {params.boundaries} \
                    -o {output.gen} {params.stringImpute}
         """

rule checkChunks:
    input:
        chunks2check
    params:
        chunks=lambda wildcards: chunksPerChromosome(chromosome=[wildcards.chr],
            dir=config['referencedir'],
            noSnpInInterval=noSnpInInterval)['chunks']
    output:
        warnings="{dir}/imputed/chr{chr}/{name}.chr{chr}.warnings_overview"
    shell:
        """
        dir={wildcards.dir}/imputed/chr{wildcards.chr}
        for chunk in {params.chunks}; do
            summary={wildcards.name}.chr{wildcards.chr}.$chunk.gen_summary
            warnings={wildcards.name}.chr{wildcards.chr}.$chunk.gen_warnings
            if grep -q 'no SNPs in the imputation' $dir/$summary; then \
                echo "{wildcards.chr}\t$chunk" >> \
                    $dir/{wildcards.name}.chr{wildcards.chr}.nochunks
            elif grep -q 'ERROR' $dir/$summary; then
                echo "{wildcards.chr}\t$chunk" >> \
                    $dir/{wildcards.name}.chr{wildcards.chr}.errors_overview
            fi
            echo "{wildcards.chr}\t$chunk" >> {output.warnings}
            cat $dir/$warnings >> {output.warnings}
        done
        """

rule combineCheckChunks:
    input:
        warnings=expand("{{dir}}/imputed/chr{chr}/{{name}}.chr{chr}.warnings_overview",
            chr=AUTOSOMES)
    output:
        warnings="{dir}/imputed/{name}.warnings_overview"
    wildcard_constraints:
        name="[\.\w]*"
    shell:
        """
        for chr in `seq 1 22`; do
            dir={wildcards.dir}/imputed/chr$chr
            if [ -f $dir/{wildcards.name}.chr$chr.nochunks ]; then
                cat $dir/{wildcards.name}.chr$chr.nochunks >> \
                    {wildcards.dir}/imputed/{wildcards.name}.nochunks
            fi
            if [ -f $dir/{wildcards.name}.chr$chr.errors_overview ] ; then
                cat  $dir/{wildcards.name}.chr$chr.errors_overview  >> \
                    {wildcards.dir}/imputed/{wildcards.name}.errors_overview
            fi
            cat {input.warnings} >> {output.warnings}
        done
        """

rule concatenateChunks:
    input:
        chunks2concatenate
    output:
        chr="{dir}/genotypes/{name}.chr{chr}.gen",
        sample="{dir}/genotypes/{name}.chr{chr}.sample"
    shell:
        """
        impute={wildcards.dir}/imputed/chr{wildcards.chr}
        phase={wildcards.dir}/phased
        cat $impute/{wildcards.name}.chr{wildcards.chr}.*.gen > {output.chr}
        awk 'NR==2 {$7="P"} {print $0}' \
            $phased/{wildcards.name}.chr{wildcards.chr}.sample > {output.sample}
        """

rule snpQC:
    input:
        gen="{dir}/genotypes/{name}.chr{chr}.gen",
    output:
        stats="{dir}/genotypes/{name}.chr{chr}.gen.stats",
        infofail="{dir}/genotypes/{name}.chr{chr}.gen.infofail",
        maffail="{dir}/genotypes/{name}.chr{chr}.gen.maffail",
        hwefail="{dir}/genotypes/{name}.chr{chr}.gen.hwefail",
        missingfail="{dir}/genotypes/{name}.chr{chr}.gen.missingfail",
        fail="{dir}/genotypes/{name}.chr{chr}.gen.fail",
        genqc="{dir}/genotypes/{name}.chr{chr}.qc.gen.gz",
    params:
        rsid=qctool_rsid_column,
        info=qctool_info_column,
        infoTh=config['infoTh'],
        maf=qctool_maf_column,
        mafTh=config['mafTh'],
        hwe=qctool_hwe_column,
        hweTh=config['hweTh'],
        missing=qctool_missing_column,
        missingTh=config['missingTh']
    shell:
        """
        qctool -g {input.gen} -snp-stats -osnp {output.stats}
        awk -F' ' '!/^#/ && (${params.info} < {params.infoTh}) \
            {{print ${params.rsid}}}' {output.stats} > {output.infofail}
        awk -F' ' '!/^#/ && (${params.maf} < {params.mafTh}) \
            {{print ${params.rsid}}}' {output.stats} > {output.maffail}
        awk -F' ' '!/^#/ && (${params.hwe} < {params.hweTh}) \
            {{print ${params.rsid}}}' {output.stats} > {output.hwefail}
        awk -F' ' '!/^#/ && (${params.missing} > {params.missingTh}) \
            {{print ${params.rsid}}}' {output.stats} > {output.missingfail}
        cat {output.infofail} {output.maffail} {output.missingfail} \
            {output.hwefail} | sort | uniq > {output.fail}
        qctool -g {input.gen} -og {output.genqc} -excl-rsids {output.fail}
        """

rule combine_and_convert:
    input:
        genqc=expand("{{dir}}/genotypes/{{name}}.chr{chr}.qc.gen.gz",
            chr=AUTOSOMES),
        sample=expand("{{dir}}/genotypes/{{name}}.chr{chr}.sample",
            chr=AUTOSOMES)
    output:
        bgen="{dir}/bgen/{name}.genome.qc.bgen",
        bimbam="{dir}/bimbam/{name}.genome.qc.dosage.gz"
    shell:
        """
        qctool -g {wildcards.dir}/genotypes/{wildcards.name}.chr#.gen \
            -og {output.bgen} -s {input.sample}
        qctool -g {wildcards.dir}/genotypes/{wildcards.name}.chr#.gen \
            -og {output.bimbam}
        """
