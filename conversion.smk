# snakemake -s conversion.smk --jobs 5000 --latency-wait 30 --cluster-config cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete

import pdb
configfile: "config_conversion.yaml"

rule all:
    input:
        expand("{ukb}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.bed",
            ukb=config["ukbdir"],
            maf=config["maf"],
            chr=range(1,23)),
        expand("{ukb}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicateIDs",
            ukb=config["ukbdir"],
            maf=config["maf"],
            chr=range(1,23)),
        expand("{ukb}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicatePos",
            ukb=config["ukbdir"],
            maf=config["maf"],
            chr=range(1,23)),
        expand("{ukb}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicatePosIDs",
            ukb=config["ukbdir"],
            maf=config["maf"],
            chr=range(1,23)),
        expand("{ukb}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.allDuplicates",
            ukb=config["ukbdir"],
            maf=config["maf"],
            chr=range(1,23)),
        expand("{ukb}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned_no_duplicates.bed",
            ukb=config["ukbdir"],
            maf=config["maf"],
            chr=range(1,23)),
        expand("{ukb}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
            ukb=config["ukbdir"],
            maf=config["maf"]),
        expand("{ukb}/popstructure/ukb_imp_v3_kinship.rel",
            ukb=config["ukbdir"]),
        expand("{ukb}/popstructure/ukb_imp_v3_pca",
            ukb=config["ukbdir"])

rule getSamples:
    input:
        pheno=expand("{ukb}/rawdata/FD.csv", ukb=config["ukbdir"]),
    output:
        "{dir}/FD_samples.csv"
    log:
        "{dir}/log/FD_samples.log"
    shell:
        "(cut -d ',' -f 1 {input.pheno} | "
        "tail -n +2 > {output}) 2> {log}"

rule convertPlink:
    input:
        gen=expand("{genodir}/ukb_imp_chr{{chr}}_v3.bgen",
            genodir=config["genodir"]),
        samples=expand("{ukb}/rawdata/ukb18545_imp_chr1_v3_s487378.sample",
            ukb=config["ukbdir"]),
        filter=expand("{ukb}/rawdata/FD_samples.csv",
            ukb=config["ukbdir"])
    output:
        "{dir}/ukb_imp_chr{chr}_v3.bed"
    log:
        "{dir}/log/chr{chr}_convert.log"
    shell:
        "(plink2 --bgen {input.gen} \
            --make-bed \
            --sample {input.samples} \
            --keep {input.filter} \
            --out {wildcards.dir}/ukb_imp_chr{wildcards.chr}_v3) 2> {log}"

rule filterMaf:
    input:
        bed=expand("{ukb}/plink/ukb_imp_chr{{chr}}_v3.bed",
            ukb=config["ukbdir"]),
        bim=expand("{ukb}/plink/ukb_imp_chr{{chr}}_v3.bim",
            ukb=config["ukbdir"]),
        fam=expand("{ukb}/plink/ukb_imp_chr{{chr}}_v3.fam",
            ukb=config["ukbdir"])
    output:
        "{dir}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.bed"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --make-bed \
            --maf {wildcards.maf} \
            --out {wildcards.dir}/maf{wildcards.maf}/ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}"


rule computeLD:
    input:
        bed=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.bed",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        bim=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.bim",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        fam=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.fam",
            maf=config["maf"],
            ukb=config["ukbdir"])
    output:
        "{dir}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.prune.in"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --indep-pairwise 50kb 1 0.8 \
            --out {wildcards.dir}/maf{wildcards.maf}/ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}"

rule pruneLD:
    input:
        bed=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.bed",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        bim=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.bim",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        fam=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.fam",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        prunein=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.prune.in",
            maf=config["maf"],
            ukb=config["ukbdir"])
    output:
        "{dir}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.bed"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --extract {input.prunein} \
            --make-bed \
            --out {wildcards.dir}/maf{wildcards.maf}/ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}.pruned"
#echo -e "$outdir/ukb_imp_chr${chr}_v3" >> $outdir/file_list

rule findDuplicates:
    input:
        bim=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.pruned.bim",
            maf=config["maf"],
            ukb=config["ukbdir"])
    output:
        duplicateIDs="{dir}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicateIDs",
        duplicatePos="{dir}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicatePos",
        duplicatePosIDs="{dir}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicatePosIDs",
        allDuplicates="{dir}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.allDuplicates"
    shell:
        """
        awk -F "\\t" 'id[$2]++ >= 1 {{print $2}}' {input.bim} > {output.duplicateIDs}
        awk -F "\\t" 'pos[$4]++ >= 1 {{print $4}}' {input.bim} > {output.duplicatePos}
        grep -w -f {output.duplicatePos} {input.bim} | cut -f 2 > {output.duplicatePosIDs}
        cat {output.duplicatePosIDs} {output.duplicateIDs} | sort | uniq > {output.allDuplicates}
        """

rule filterDuplicates:
    input:
        duplicates=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.pruned.allDuplicates",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        bed=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.pruned.bed",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        bim=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.pruned.bim",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        fam=expand("{ukb}/maf{maf}/ukb_imp_chr{{chr}}_v3_maf{maf}.pruned.fam",
            maf=config["maf"],
            ukb=config["ukbdir"])
    output:
        "{dir}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned_no_duplicates.bed"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --exclude {input.duplicates} \
            --make-bed \
            --out {wildcards.dir}/maf{wildcards.maf}/ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}.pruned_no_duplicates"

rule createFilelist:
    output:
        "'{dir}/maf{maf}/ukb_imp_v3_maf{maf}.pruned_no_duplicates.file_list'.format(dir=config['ukbdir'], maf=config['maf'])"
    run:
        files2merge = ["{dir}/maf{maf}/ukb_imp_chr{chr}_v3_maf{maf}.pruned_no_duplicates".format(
            dir=config['ukbdir'], maf=config['maf'], chr=x) for x in range(1, 23)]

        with open({output}, "w") as handle:
            for f in files2merge:
                handle.write("{}\n".format(f))

rule mergeFiles:
    input:
        list=expand("{ukb}/maf{maf}/ukb_imp_v3_maf{maf}.pruned_no_duplicates.file_list",
            ukb=config["ukbdir"],
            maf=config["maf"])
    output:
        "{dir}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.bed"
    shell:
        "plink --merge-list {input.list} \
            --make-bed \
            --out {wildcards.dir}/maf{wildcards.maf}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned"

rule kinship:
    input:
        bed=expand("{ukb}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        bim=expand("{ukb}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        fam=expand("{ukb}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
            maf=config["maf"],
            ukb=config["ukbdir"])
    output:
        "{dir}/ukb_imp_v3_kinship.rel"
    shell:
        "plink2 --bed {input.bed} \
                --fam {input.fam} \
                --bim {input.bim} \
                --maf 0.1 \
                --make-rel square \
                --out {wildcards.dir}/ukb_imp_v3_kinship.rel"

rule pca:
    input:
        bed=expand("{ukb}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        bim=expand("{ukb}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
            maf=config["maf"],
            ukb=config["ukbdir"]),
        fam=expand("{ukb}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
            maf=config["maf"],
            ukb=config["ukbdir"])
    output:
        "{dir}/ukb_imp_v3_pca"
    shell:
        "flashpca --bed {input.bed} \
                --fam {input.fam} \
                --bim {input.bim} \
                --ndim 50 \
                --outpc {wildcards.dir}/ukb_imp_v3_pca \
                --suffix _ukb_imp_v3.txt"

