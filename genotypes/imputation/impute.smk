from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

###### load config file #####
configfile: "config/config_impute.yaml"


###### target functions  #####
def chunksPerChromosome(chromosome, dir, noSnpInInterval=None):
    if noSnpInInterval is not None:
        nosnp = noSnpInInterval.split(",")
        exclude = dict(item.split(".") for item in nosnp)
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

# global parameters
#chromosomes = list(range(1,23)) + ["X"]
chromosomes = list(range(1,23))
perChr = chunksPerChromosome(chromosome=chromosomes, dir=config['referencedir'],
    noSnpInInterval=config['noSnpInInterval'])

##### target rules #####
rule all:
    input:
        # rules/prepare.smk
        expand("{dir}/unphased/{name}.chr{chr}.{suffix}",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23),
            suffix=['bim', 'bed', 'fam']),
        # rules/phase.smk
        expand("{dir}/phased/{name}.chr{chr}.{suffix}",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23),
            suffix=['hap.gz', 'sample']),
        # rules/impute.smk
        expand(expand("{dir}/imputed/chr{{chr}}/{name}.chr{{chr}}.{{chunk}}.{suffix}",
                dir=config["dir"],
                name='gencall.combined.clean.related',
                suffix=['gen']),
            zip, chunk=perChr['chunks'], chr=perChr['chr']),
        # rules/check.smk
        expand("{dir}/imputed/chr{chr}/{name}.chr{chr}.warnings_overview",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23)),
        expand("{dir}/imputed/{name}.warnings_overview",
            dir=config["dir"],
            name='gencall.combined.clean.related'),
        expand("{dir}/genotypes/{name}.chr{chr}.{suffix}",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23),
            suffix=['sample', 'gen']),
        # rules/qc.smk
        expand("{dir}/genotypes/{name}.chr{chr}.qc.gen.gz",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23)),
        # rule/combine.smk
        expand("{dir}/formated/{name}.chr{chr}.qc.{suffix}",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            chr=range(1,23),
            suffix=['bgen', 'dosage.gz']),
        expand("{dir}/formated/{name}.genome.qc.{suffix}",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            suffix=['bgen', 'dosage.gz']),
        # rules/counts.smk
        expand("{dir}/counts/{name}.SNPsPerChr.{suffix}",
            dir=config["dir"],
            name='gencall.combined.clean.related',
            suffix=['txt', 'pdf'])

##### load rules #####
include: "rules/prepare.smk"
include: "rules/phase.smk"
include: "rules/impute.smk"
include: "rules/check.smk"
include: "rules/qc.smk"
include: "rules/convert.smk"
include: "rules/count.smk"

##### example cluster call #####
# snakemake -s impute.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

