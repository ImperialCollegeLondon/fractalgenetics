## QC columns in qctool -snp-stats output
qctool_missing_column = 19
qctool_info_column = 17
qctool_maf_column = 14
qctool_hwe_column = 8
qctool_rsid_column = 2

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
