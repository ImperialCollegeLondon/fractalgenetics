# snakemake -s genotypes.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue}
# -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete


configfile: "config/config_genotypes.yaml"

rule all:
    input:
        expand("{ukb}/rawdata/{pheno}_samples.csv",
            ukb=config["ukbdir"],
            pheno=config["pheno"]),
        expand("{ukb}/plink/{pheno}/ukb_imp_chr{chr}_v3.bed",
             ukb=config["ukbdir"],
             maf=config["maf"],
             pheno=config["pheno"],
             chr=range(1,23)),
        expand("{ukb}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.bim",
             ukb=config["ukbdir"],
             maf=config["maf"],
             pheno=config["pheno"],
             chr=range(1,23)),
        expand("{ukb}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicateIDs",
             ukb=config["ukbdir"],
             maf=config["maf"],
             pheno=config["pheno"],
             chr=range(1,23)),
        expand("{ukb}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned_no_duplicates.bed",
             ukb=config["ukbdir"],
             maf=config["maf"],
             pheno=config["pheno"],
             chr=range(1,23)),
         expand("{ukb}/maf{maf}/{pheno}/ukb_imp_v3_kinship.rel",
             pheno=config["pheno"],
             maf=config["maf"],
             ukb=config["ukbdir"]),
         expand("{ukb}/maf{maf}/{pheno}/ukb_imp_v3_pca",
             pheno=config["pheno"],
             maf=config["maf"],
             ukb=config["ukbdir"])

rule getSamples:
    input:
        pheno="{dir}/rawdata/{pheno}.csv",
    output:
        "{dir}/rawdata/{pheno}_samples.csv"
    log:
        "{dir}/log/{pheno}_samples.log"
    shell:
        "(cut -d ',' -f 1 {input.pheno} | "
        "tail -n +2 > {output}) 2> {log}"

rule convertPlink:
    input:
        gen=expand("{genodir}/ukb_imp_chr{{chr}}_v3.bgen",
            genodir=config["genodir"]),
        samples="{dir}/rawdata/ukb18545_imp_chr1_v3_s487378.sample",
        filter="{dir}/rawdata/{pheno}_samples.csv"
    output:
        "{dir}/plink/{pheno}/ukb_imp_chr{chr}_v3.bed"
        "{dir}/plink/{pheno}/ukb_imp_chr{chr}_v3.bim"
        "{dir}/plink/{pheno}/ukb_imp_chr{chr}_v3.fam"
    log:
        "{dir}/plink/{pheno}/log/chr{chr}_convert.log"
    shell:
        "(plink2 --bgen {input.gen} \
            --make-bed \
            --sample {input.samples} \
            --keep {input.filter} \
            --out {wildcards.dir}/plink/{wildcards.pheno}/ukb_imp_chr{wildcards.chr}_v3) 2> {log}"

rule filterMaf:
    input:
        bed="{dir}/plink/{pheno}/ukb_imp_chr{chr}_v3.bed",
        bim="{dir}/plink/{pheno}/ukb_imp_chr{chr}_v3.bim",
        fam="{dir}/plink/{pheno}/ukb_imp_chr{chr}_v3.fam",
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.bed"
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.bim"
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.fam"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --make-bed \
            --maf {wildcards.maf} \
            --out {wildcards.dir}/maf{wildcards.maf}/{wildcards.pheno}/ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}"


rule computeLD:
    input:
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.fam",
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.prune.in"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --indep-pairwise 50kb 1 0.8 \
            --out {wildcards.dir}/maf{wildcards.maf}/{wildcards.pheno}/ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}"

rule pruneLD:
    input:
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.fam",
        prunein="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.prune.in",
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.bed",
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.bim",
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.fam"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --extract {input.prunein} \
            --make-bed \
            --out {wildcards.dir}/maf{wildcards.maf}/{wildcards.pheno}/ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}.pruned"

rule findDuplicates:
    input:
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.bim",
    output:
        duplicateIDs="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicateIDs",
        duplicatePos="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicatePos",
        duplicatePosIDs="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.duplicatePosIDs",
        allDuplicates="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.allDuplicates"
    shell:
        """
        awk -F "\\t" 'id[$2]++ >= 1 {{print $2}}' {input.bim} > {output.duplicateIDs}
        awk -F "\\t" 'pos[$4]++ >= 1 {{print $4}}' {input.bim} > {output.duplicatePos}
        grep -w -f {output.duplicatePos} {input.bim} | cut -f 2 > {output.duplicatePosIDs}
        cat {output.duplicatePosIDs} {output.duplicateIDs} | sort | uniq > {output.allDuplicates}
        """

rule filterDuplicates:
    input:
        duplicates="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.allDuplicates",
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned.fam"
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned_no_duplicates.bed",
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned_no_duplicates.bim",
        "{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned_no_duplicates.fam"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --exclude {input.duplicates} \
            --make-bed \
            --out {wildcards.dir}/maf{wildcards.maf}/{wildcards.pheno}/ukb_imp_chr{wildcards.chr}_v3_maf{wildcards.maf}.pruned_no_duplicates"

rule createFilelist:
    input:
        bed=expand("{{dir}}/maf{{maf}}/{{pheno}}/ukb_imp_chr{chr}_v3_maf{{maf}}.pruned_no_duplicates.bed",
             chr=range(1,23)),
        bim=expand("{{dir}}/maf{{maf}}/{{pheno}}/ukb_imp_chr{chr}_v3_maf{{maf}}.pruned_no_duplicates.bim",
             chr=range(1,23)),
        fam=expand("{{dir}}/maf{{maf}}/{{pheno}}/ukb_imp_chr{chr}_v3_maf{{maf}}.pruned_no_duplicates.fam",
             chr=range(1,23)),
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_v3_maf{maf}.pruned_no_duplicates.file_list"
    run:
        files2merge = ["{dir}/maf{maf}/{pheno}/ukb_imp_chr{chr}_v3_maf{maf}.pruned_no_duplicates".format(
            dir=config['ukbdir'], maf=config['maf'], pheno=config['pheno'], chr=x) for x in range(1, 23)]

        with open(output[0], "w") as handle:
            for f in files2merge:
                handle.write("{}\n".format(f))

rule mergeFiles:
    input:
        list="{dir}/maf{maf}/{pheno}/ukb_imp_v3_maf{maf}.pruned_no_duplicates.file_list",
        bed=expand("{{dir}}/maf{{maf}}/{{pheno}}/ukb_imp_chr{chr}_v3_maf{{maf}}.pruned_no_duplicates.bed",
             chr=range(1,23)),
        bim=expand("{{dir}}/maf{{maf}}/{{pheno}}/ukb_imp_chr{chr}_v3_maf{{maf}}.pruned_no_duplicates.bim",
             chr=range(1,23)),
        fam=expand("{{dir}}/maf{{maf}}/{{pheno}}/ukb_imp_chr{chr}_v3_maf{{maf}}.pruned_no_duplicates.fam",
             chr=range(1,23)),
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        "{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        "{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.fam"
    shell:
        "plink --merge-list {input.list} \
            --make-bed \
            --out {wildcards.dir}/maf{wildcards.maf}/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned"

rule kinship:
    input:
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_v3_kinship.rel"
    shell:
        "plink2 --bed {input.bed} \
                --fam {input.fam} \
                --bim {input.bim} \
                --maf 0.1 \
                --make-rel square \
                --out {wildcards.dir}/{wildcards.pheno}/ukb_imp_v3_kinship.rel"

rule pca:
    input:
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_v3_pca"
    shell:
        "flashpca --bed {input.bed} \
                --fam {input.fam} \
                --bim {input.bim} \
                --ndim 50 \
                --outpc {wildcards.dir}/{wildcards.pheno}/ukb_imp_v3_pca \
                --suffix _ukb_imp_v3.txt"
