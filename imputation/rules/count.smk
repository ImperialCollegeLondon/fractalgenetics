rule genotypeCounts:
    input:
        chr=expand("{{dir}}/genotypes/{{name}}.chr{chr}.gen",
            chr=chromosomes),
        chrQC=expand("{{dir}}/genotypes/{{name}}.chr{chr}.qc.gen.gz",
            chr=chromosomes)
    output:
        txt="{dir}/counts/{name}.SNPsPerChr.txt",
        pdf="{dir}/counts/{name}.SNPsPerChr.pdf"
    params:
        refdir=config['referencedir'],
        nosnp=config['noSnpInInterval']
    shell:
        """
        perl scripts/imputeQC.pl --indir {wildcards.dir} \
            --outdir {wildcards.dir}/counts \
            --refdir {params.refdir} \
            --nosnp {params.nosnp} \
            --name {wildcards.name} \
            --verbose
        Rscript scripts/imputeQC.R \
            --directory {wildcards.dir}/counts \
            --name {wildcards.name}
        """
