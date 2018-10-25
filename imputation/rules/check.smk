def chunks2concatenate(wildcards):
    return expand("{{dir}}/imputed/chr{{chr}}/{{name}}.chr{{chr}}.{chunk}.gen",
                chunk=chunksPerChromosome(chromosome=[wildcards.chr],
                    dir=config['referencedir'],
                    noSnpInInterval=config['noSnpInInterval'])['chunks'])

def chunks2check(wildcards):
    warnings=expand("{{dir}}/imputed/chr{{chr}}/{{name}}.chr{{chr}}.{chunk}.gen_warnings",
                chunk=chunksPerChromosome(chromosome=[wildcards.chr],
                    dir=config['referencedir'],
                    noSnpInInterval=config['noSnpInInterval'])['chunks'])
    return warnings

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

rule checkChunks:
    input:
        chunks2check
    params:
        chunks=lambda wildcards: chunksPerChromosome(chromosome=[wildcards.chr],
            dir=config['referencedir'],
            noSnpInInterval=config['noSnpInInterval'])['chunks']
    output:
        warnings="{dir}/imputed/chr{chr}/{name}.chr{chr}.warnings_overview"
    shell:
        """
        dir={wildcards.dir}/imputed/chr{wildcards.chr}
        for chunk in {params.chunks}; do
            summary={wildcards.name}.chr{wildcards.chr}.$chunk.gen_summary
            warnings={wildcards.name}.chr{wildcards.chr}.$chunk.gen_warnings
            if grep -q 'ERROR: There are no' $dir/$summary; then \
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
            chr=chromosomes)
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
        phased={wildcards.dir}/phased
        cat $impute/{wildcards.name}.chr{wildcards.chr}.*.gen > {output.chr}
        awk 'NR==2 {{$7="P"}} {{print $0}}' \
            $phased/{wildcards.name}.chr{wildcards.chr}.sample > {output.sample}
        """
