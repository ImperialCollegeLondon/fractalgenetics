# how to run from ./:
# snakemake -s ancestry.smk --jobs 5000 --latency-wait 200 --cluster-config config/cluster.json 
# --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete

configfile: "config/config_ancestry.yaml"

rule all:
    input:
        expand("{ukb}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
            ukb=config["ukbdir"],
            pheno=config["pheno"],
            maf=config["maf"]),
        expand("{ukb}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.kinship.rel",
            ukb=config["ukbdir"],
            pheno=config["pheno"],
            maf=config["maf"]),
        expand("{ukb}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.pca",
            ukb=config["ukbdir"],
            pheno=config["pheno"],
            maf=config["maf"]),
        expand("{ukb}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.pca",
            ukb=config["ukbdir"],
            pheno=config["pheno"],
            maf=config["maf"],
            hapmap=config["hapmap"]),
        expand("{ukb}/ancestry/{pheno}/{hapmap}_{pheno}_pca.png",
            pheno=config["pheno"],
            ukb=config["ukbdir"],
            hapmap=config["hapmap"]),
        expand("{ukb}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.European.pca",
            pheno=config["pheno"],
            ukb=config["ukbdir"],
            maf=config["maf"])

# maf in
rule ukbSNPs:
    input:
        bed=expand("{ukb}/maf{maf}/{{pheno}}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
            ukb=config["ukbdir"],
            maf=config['maf-in']),
        bim=expand("{ukb}/maf{maf}/{{pheno}}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
            ukb=config["ukbdir"],
            maf=config['maf-in']),
        fam=expand("{ukb}/maf{maf}/{{pheno}}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
            ukb=config["ukbdir"],
            maf=config['maf-in'])
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        "{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        "{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.snplist"
    wildcard_constraints:
        maf=0.1
    shell:
        """
        plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --make-bed \
            --maf {wildcards.maf} \
            --write-snplist \
            --out {wildcards.dir}/maf{wildcards.maf}/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned
        cp {wildcards.dir}/maf{wildcards.maf}/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.snplist \
           {wildcards.dir}/ancestry/{wildcards.pheno}
        """

rule kinshipUKB:
    input:
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.fam"
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.kinship.rel"
    shell:
        "plink2 --bed {input.bed} \
                --fam {input.fam} \
                --bim {input.bim} \
                --make-rel square \
                --out {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.kinship"

rule pcaUKB:
    input:
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.fam"
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.pca"
    shell:
        "flashpca --bed {input.bed} \
                --fam {input.fam} \
                --bim {input.bim} \
                --ndim 50 \
                --outpc {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.pca \
                --outvec  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.eigenvectors \
                --outload  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.SNPloadings \
                --outval  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.eigenvalues \
                --outpve  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.variance_explained"

rule extractUKBfromHapmap:
    input:
        bed=expand("{hapdir}/{hapmap}.bed",
            hapmap=config["hapmap"],
            hapdir=config["hapdir"]),
        bim=expand("{hapdir}/{hapmap}.bim",
            hapmap=config["hapmap"],
            hapdir=config["hapdir"]),
        fam=expand("{hapdir}/{hapmap}.fam",
            hapmap=config["hapmap"],
            hapdir=config["hapdir"]),
        snplist=expand("{{dir}}/ancestry/{{pheno}}/ukb_imp_genome_v3_maf{maf}.pruned.snplist",
            maf=config['maf'])
    output:
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.bed",
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.bim",
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.fam",
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.snplist"
    wildcard_constraints:
        hapmap=config['hapmap']
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --make-bed \
            --extract {input.snplist} \
            --write-snplist \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/{wildcards.hapmap}.intersection"

rule extractHapmapFromUKB:
    input:
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
        snplist=expand("{{dir}}/ancestry/{{pheno}}/{hapmap}.intersection.snplist",
            hapmap=config["hapmap"])
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.bed",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.bim",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.fam"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --make-bed \
            --extract {input.snplist} \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.intersection"

rule compareSNPs:
    input:
        ukb_bed="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.bed",
        ukb_bim="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.bim",
        ukb_fam="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.fam",
        hapmap_bed="{dir}/ancestry/{pheno}/{hapmap}.intersection.bed",
        hapmap_bim="{dir}/ancestry/{pheno}/{hapmap}.intersection.bim",
        hapmap_fam="{dir}/ancestry/{pheno}/{hapmap}.intersection.fam",
    output:
        "{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_genome_v3_maf{maf}.pruned.missnp"
    wildcard_constraints:
        hapmap=config['hapmap']
    shell:
        """
        set +e
        plink --bed {input.hapmap_bed} \
            --bim {input.hapmap_bim} \
            --fam {input.hapmap_fam} \
            --bmerge {input.ukb_bed} {input.ukb_bim} {input.ukb_fam} \
            --merge-mode 6 \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/compare.{wildcards.hapmap}.ukb_imp_genome_v3_maf{wildcards.maf}.pruned
        set -e
        """

rule excludeMissnpUKB:
    input:
        ukb_bed="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.bed",
        ukb_bim="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.bim",
        ukb_fam="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.fam",
        missnp=expand("{{dir}}/ancestry/{{pheno}}/compare.{hapmap}.ukb_imp_genome_v3_maf{{maf}}.pruned.missnp",
            hapmap=config['hapmap'])
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.nomissnp.bed",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.nomissnp.bim",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.nomissnp.fam"
    shell:
        "plink --bed {input.ukb_bed} \
            --bim {input.ukb_bim} \
            --fam {input.ukb_fam} \
            --exclude {input.missnp} \
            --make-bed \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.intersection.nomissnp"

rule excludeMissnpHapmap:
    input:
        hapmap_bed="{dir}/ancestry/{pheno}/{hapmap}.intersection.bed",
        hapmap_bim="{dir}/ancestry/{pheno}/{hapmap}.intersection.bim",
        hapmap_fam="{dir}/ancestry/{pheno}/{hapmap}.intersection.fam",
        missnp=expand("{{dir}}/ancestry/{{pheno}}/compare.{{hapmap}}.ukb_imp_genome_v3_maf{maf}.pruned.missnp",
            maf=config['maf'])
    output:
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.nomissnp.bed",
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.nomissnp.bim",
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.nomissnp.fam",
    wildcard_constraints:
        hapmap=config['hapmap']
    shell:
        "plink --bed {input.hapmap_bed} \
            --bim {input.hapmap_bim} \
            --fam {input.hapmap_fam} \
            --exclude {input.missnp} \
            --make-bed \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/{wildcards.hapmap}.intersection.nomissnp"

rule findMismatches:
    input:
        ukb_bim="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.nomissnp.bim",
        hapmap_bim="{dir}/ancestry/{pheno}/{hapmap}.intersection.nomissnp.bim",
    output:
        flippedAlleles="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_genome_v3_maf{maf}.pruned.flippedAlleles",
        same="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_genome_v3_maf{maf}.pruned.same",
        keep="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_genome_v3_maf{maf}.pruned.keep"
    wildcard_constraints:
        hapmap=config['hapmap']
    shell:
        """
        awk 'FNR==NR {{a[$1$2$4$5$6]; next}} $1$2$4$5$6 in a' {input.ukb_bim} {input.hapmap_bim} > {output.same}
        awk 'FNR==NR {{a[$2$5$6]; next}} $2$6$5 in a {{print $2\t$6\t$5\t$5\t$6}}' {input.ukb_bim} {input.hapmap_bim} > {output.flippedAlleles}
        cut -f 2 {output.same} {output.flippedAlleles} > {output.keep}
        """

rule filterUKB:
    input:
        ukb_bed="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.nomissnp.bed",
        ukb_bim="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.nomissnp.bim",
        ukb_fam="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.nomissnp.fam",
        keep=expand("{{dir}}/ancestry/{{pheno}}/compare.{hapmap}.ukb_imp_genome_v3_maf{{maf}}.pruned.keep",
            hapmap=config['hapmap'])
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.bed",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.bim",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.fam"
    shell:
        "plink --bed {input.ukb_bed} \
            --bim {input.ukb_bim} \
            --fam {input.ukb_fam} \
            --extract {input.keep} \
            --make-bed \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.intersection.matched"

rule filterAndFlipHapmap:
    input:
        hapmap_bed="{dir}/ancestry/{pheno}/{hapmap}.intersection.nomissnp.bed",
        hapmap_bim="{dir}/ancestry/{pheno}/{hapmap}.intersection.nomissnp.bim",
        hapmap_fam="{dir}/ancestry/{pheno}/{hapmap}.intersection.nomissnp.fam",
        keep=expand("{{dir}}/ancestry/{{pheno}}/compare.{{hapmap}}.ukb_imp_genome_v3_maf{maf}.pruned.keep",
            maf=config['maf']),
        flipped=expand("{{dir}}/ancestry/{{pheno}}/compare.{{hapmap}}.ukb_imp_genome_v3_maf{maf}.pruned.flippedAlleles",
            maf=config['maf'])
    output:
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.matched.bed",
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.matched.bim",
        "{dir}/ancestry/{pheno}/{hapmap}.intersection.matched.fam",
    wildcard_constraints:
        hapmap=config['hapmap']
    shell:
        "plink --bed {input.hapmap_bed} \
            --bim {input.hapmap_bim} \
            --fam {input.hapmap_fam} \
            --extract {input.keep} \
            --update-alleles {input.flipped} \
            --make-bed \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/{wildcards.hapmap}.intersection.matched"

rule mergeHapMap:
    input:
        ukb_bed="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.bed",
        ukb_bim="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.bim",
        ukb_fam="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.fam",
        hapmap_bed="{dir}/ancestry/{pheno}/{hapmap}.intersection.matched.bed",
        hapmap_bim="{dir}/ancestry/{pheno}/{hapmap}.intersection.matched.bim",
        hapmap_fam="{dir}/ancestry/{pheno}/{hapmap}.intersection.matched.fam",
    output:
        "{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.bed",
        "{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.bim",
        "{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.fam"
    shell:
        "plink --bed {input.hapmap_bed} \
            --bim {input.hapmap_bim} \
            --fam {input.hapmap_fam} \
            --bmerge {input.ukb_bed} {input.ukb_bim} {input.ukb_fam} \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/{wildcards.hapmap}_ukb_imp_genome_v3_maf{wildcards.maf}.pruned.merged"

rule pcaMerge:
    input:
        bed="{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.bed",
        bim="{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.bim",
        fam="{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.fam"
    output:
        "{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.pca"
    shell:
        "flashpca --bed {input.bed} \
                --fam {input.fam} \
                --bim {input.bim} \
                --ndim 10 \
                --outpc {wildcards.dir}/ancestry/{wildcards.pheno}/{wildcards.hapmap}_ukb_imp_genome_v3_maf{wildcards.maf}.pruned.merged.pca \
                --outvec  {wildcards.dir}/ancestry/{wildcards.pheno}/{wildcards.hapmap}_ukb_imp_genome_v3_maf{wildcards.maf}.pruned.merged.eigenvectors \
                --outload  {wildcards.dir}/ancestry/{wildcards.pheno}/{wildcards.hapmap}_ukb_imp_genome_v3_maf{wildcards.maf}.pruned.merged.SNPloadings \
                --outval  {wildcards.dir}/ancestry/{wildcards.pheno}/{wildcards.hapmap}_ukb_imp_genome_v3_maf{wildcards.maf}.pruned.merged.eigenvalues \
                --outpve  {wildcards.dir}/ancestry/{wildcards.pheno}/{wildcards.hapmap}_ukb_imp_genome_v3_maf{wildcards.maf}.pruned.merged.variance_explained"

rule plotPCA:
    input:
        samples=expand("{{dir}}/ancestry/{{pheno}}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.fam",
            maf=config['maf']),
        pcafile=expand("{{dir}}/ancestry/{{pheno}}/{{hapmap}}_ukb_imp_genome_v3_maf{maf}.pruned.merged.pca",
            maf=config['maf'])
    output:
        "{dir}/ancestry/{pheno}/{hapmap}_{pheno}_pca.png",
        "{dir}/ancestry/{pheno}/European_samples_filtered_by_{hapmap}.txt"
    shell:
        "Rscript 'ancestry/selectPCA.R' \
            --directory={wildcards.dir}/ancestry/{wildcards.pheno} \
            --pcadata={input.pcafile} \
            --samples={input.samples} \
            --name={wildcards.pheno} \
            --reference={wildcards.hapmap}"

rule filterEuropeans:
    input:
        keep=expand("{{dir}}/ancestry/{{pheno}}/European_samples_filtered_by_{hapmap}.txt",
            hapmap=config['hapmap']),
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.European.fam",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.European.bim",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.European.bed"
    shell:
        "plink --bed {input.bed} \
            --fam {input.fam} \
            --bim {input.bim} \
            --keep {input.keep} \
            --make-bed \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European"

rule pcaEuropean:
    input:
        bed="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.European.bed",
        bim="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.European.bim",
        fam="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.European.fam",
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.European.pca"
    shell:
        "flashpca --bed {input.bed} \
                --fam {input.fam} \
                --bim {input.bim} \
                --ndim 50 \
                --outpc {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.pca \
                --outvec  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.eigenvectors \
                --outload  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.SNPloadings \
                --outval  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.eigenvalues \
                --outpve  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.variance_explained"

