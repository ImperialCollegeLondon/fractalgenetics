def chunks2concatenate(wildcards):
    return expand("{{dir}}/imputed/chr{{chr}}/{{name}}.chr{{chr}}.{chunk}.gen",
                chunk=chunksPerChromosome(chromosome=[wildcards.chr],
                    dir=config['referencedir'],
                    noSnpInInterval=config['noSnpInInterval'])['chunks'])

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

rule genotypeCounts:
    input:
        chunksPerChr=chunks2concatenate,
        chr="{dir}/genotypes/{name}.chr{chr}.gen",
        chrQC="{dir}/genotypes/{name}.chr{chr}.qc.gen.gz",
    output:
        chr="{dir}/counts/chr{chr}SNPsPerChunk.txt"
    params:
        refdir=config['referencedir']
    shell:
        """
        perl PhasingAndImputation/imputeQC.pl --indir {wildcards.dir} \
            --outdir {wildcards.dir}/counts \
            --refdir {params.refdir} \
            --noSnps noSnpInIntervalString \
            --verbose
        """
