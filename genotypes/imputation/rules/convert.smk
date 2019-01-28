AUTOSOMES = list(range(1,23))

rule convert:
    input:
        genqc="{dir}/genotypes/{name}.chr{chr}.qc.gen.gz",
        sample="{dir}/genotypes/{name}.chr1.sample"
    output:
        bgen="{dir}/formated/{name}.chr{chr}.qc.bgen",
        bimbam="{dir}/formated/{name}.chr{chr}.qc.dosage.gz"
    shell:
        """
        qctool -g {input.genqc} -og {output.bgen} -s {input.sample}
        qctool -g {input.genqc} -og {output.bimbam} -s {input.sample}
        """

rule combine_and_convert:
    input:
        genqc=expand("{{dir}}/genotypes/{{name}}.chr{chr}.qc.gen.gz",
            chr=AUTOSOMES),
        sample=expand("{{dir}}/genotypes/{{name}}.chr{chr}.sample",
            chr=1)
    output:
        bgen="{dir}/formated/{name}.genome.qc.bgen",
        bimbam="{dir}/formated/{name}.genome.qc.dosage.gz"
    shell:
        """
        qctool -g {wildcards.dir}/genotypes/{wildcards.name}.chr#.qc.gen.gz \
            -og {output.bgen} -s {input.sample}
        qctool -g {wildcards.dir}/genotypes/{wildcards.name}.chr#.qc.gen.gz \
            -og {output.bimbam} -s {input.sample}
        """

