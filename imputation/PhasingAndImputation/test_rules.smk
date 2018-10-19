
        expand(expand("{dir}/imputed/chr{{chr}}/{name}.chr{{chr}}.{{chunk}}.{suffix}",
                dir=config["dir"],
                name='gencall.combined.clean.related',
                suffix=['nochunks', 'warnings_overview', 'errors_overview']),
            zip,
            chr=np.repeat(list(range(1,23)) + ["X"], CHUNKS),
            chunk=reduce(add, [list(range(1,x+1)) for x in CHUNKS])),

rule checkChunks:
    input:
        summary=expand("{{dir}}/{{name}}.chr{{chr}}.{chunk}.gen_summary",
            chunk=1),
            #chunk=glob_wildcards("{dir}/{name}.chr{chr}.{chunk}.gen_summary").chr),
        warnings=expand("{{dir}}/{{name}}.chr{{chr}}.{chunk}.gen_warnings",
            chunk=glob_wildcards("{dir}/{name}.chr{chr}.{chunk}.gen_summary").chunk)
    params:
        chunks = glob_wildcards("{dir}/{name}.chr{chr}.{chunk}.gen_summary").chunk
    output:
        nochunks="{dir}/{name}.chr{chr}_nochunks",
        warnings="{dir}/{name}.chr{chr}_warnings_overview",
        errors="{dir}/{name}.chr{chr}_errors_overview",
    shell:
        """
        for c in {params.chunks}:
            if grep -q 'no SNPs in the imputation interval' {input.summary}; then \
                echo "{wildcards.chr}\t$c" >> {output.nochunks}
            elif grep -q 'ERROR' {input.summary}; then
                echo "{wildcards.chr}\t$c" >> {output.errors}
            fi
            echo "{wildcards.chr}\t$c" >> {output.warnings}
            cat {input.warnings} >> {output.warnings}
        """
