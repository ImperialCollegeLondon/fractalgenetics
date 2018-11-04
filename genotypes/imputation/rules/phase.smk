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
