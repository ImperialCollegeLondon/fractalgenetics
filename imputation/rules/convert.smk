AUTOSOMES = list(range(1,23))

rule combine_and_convert:
    input:
        genqc=expand("{{dir}}/genotypes/{{name}}.chr{chr}.qc.gen.gz",
            chr=AUTOSOMES),
        sample=expand("{{dir}}/genotypes/{{name}}.chr{chr}.sample",
            chr=1)
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
