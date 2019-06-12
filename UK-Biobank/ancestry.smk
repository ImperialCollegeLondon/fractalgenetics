# how to run fro ./:
# snakemake -s ancestry.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources} -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete

configfile: "config/config_ancestry.yaml"

rule all:
    input:
        expand("{ukb}/ancestry/ukb_imp_genome_v3_maf{maf}.pruned.kinship.rel",
            ukb=config["ukbdir"],
            maf=config["maf"]),
        expand("{ukb}/ancestry/ukb_imp_genome_v3_maf{maf}.pruned.pca",
            ukb=config["ukbdir"],
            maf=config["maf"]),
        expand("{ukb}/ancestry/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.pca",
            ukb=config["ukbdir"],
            maf=config["maf"],
            hapmap=config["hapmap"]),
        expand("{ukb}/ancestry/HapMap_UKBB-FD_pca.png",
            ukb=config["ukbdir"]),
        expand("{ukb}/ancestry/ukb_imp_genome_v3_maf{maf}.pruned.European.pca",
            ukb=config["ukbdir"],
            maf=config["maf"])

# maf in
rule ukbSNPs:
    input:
        bed=expand("{{dir}}/maf{maf}/{{pheno}}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
            maf=config['maf-in'])
        bim=expand("{{dir}}/maf{maf}/{{pheno}}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
            maf=config['maf-in'])
        fam=expand("{{dir}}/maf{maf}/{{pheno}}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
            maf=config['maf-in'])
    output:
        "{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.snplist"
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
                --ndim 10 \
                --outpc {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.pca \
                --outvec  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.eigenvectors \
                --outload  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.SNPloadings \
                --outval  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.eigenvalues \
                --outpve  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.variance_explained"


rule extractUKBfromHapmap:
    input:
        bed=expand("{hapdir}/{hapmap}.bed",
            hapdir=config["hapdir"],
            hapmap=config["hapmap"]),
        bim=expand("{hapdir}/{hapmap}.bim",
            hapdir=config["hapdir"],
            hapmap=config["hapmap"]),
        fam=expand("{hapdir}/{hapmap}.fam",
            hapdir=config["hapdir"],
            hapmap=config["hapmap"]),
        snplist="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.snplist",
    output:
        "{dir}/ancestry/{pheno}/HapMap{hapmap}.intersection.bed",
        "{dir}/ancestry/{pheno}/HapMap{hapmap}.intersection.snplist"
    shell:
        "plink2 --bed {input.bed} \
            --bim {input.bim} \
            --fam {input.fam} \
            --make-bed \
            --extract {input.snplist} \
            --write-snplist \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/HapMap{wildcards.hapmap}.intersection"

rule extractHapmapFromUKB:
    input:
        bed="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        bim="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        fam="{dir}/maf{maf}/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
        snplist=expand("{{dir}}/ancestry/{{pheno}}/{hapmap}.intersection.snplist",
            hapmap=config["hapmap"])
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.bed"
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
        hapmap_bed="{{dir}}/ancestry/{pheno}/{hapmap}.intersection.bed",
        hapmap_bim="{dir}/ancestry/{pheno}/{hapmap}.intersection.bim",
        hapmap_fam="{dir}/ancestry/{pheno}/{hapmap}.intersection.fam",
    output:
        "{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_v3_maf{maf}.pruned.missnp"
    shell:
        """
        set +e
        plink --bed {input.hapmap_bed} \
            --bim {input.hapmap_bim} \
            --fam {input.hapmap_fam} \
            --bmerge {input.ukb_bed} {input.ukb_bim} {input.ukb_fam} \
            --merge-mode 6 \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/compare.{wildcards.hapmap}.ukb_imp_v3_maf{wildcards.maf}.pruned
        set -e
        """

rule excludeMissnpUKB:
    input:
        ukb_bed="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.bed",
        ukb_bim="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.bim",
        ukb_fam="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.fam",
        missnp="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_v3_maf{maf}.pruned.missnp",
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.nomissnp.bed"
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
        missnp="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_v3_maf{maf}.pruned.missnp",
    output:
        "{dir}/ancestry/{pheno}/HapMap{hapmap}.intersection.nomissnp.bed",
    shell:
        "plink --bed {input.hapmap_bed} \
            --bim {input.hapmap_bim} \
            --fam {input.hapmap_fam} \
            --exclude {input.missnp} \
            --make-bed \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/HapMap{wildcards.hapmap}.intersection.nomissnp"

rule findMismatches:
    input:
        ukb_bim="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.nomissnp.bim",
        hapmap_bim="{dir}/ancestry/{pheno}/{hapmap}.intersection.nomissnp.bim",
    output:
        flippedAlleles="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_v3_maf{maf}.pruned.flippedAlleles",
        same="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_v3_maf{maf}.pruned.same",
        keep="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_v3_maf{maf}.pruned.keep"
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
        keep="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_v3_maf{maf}.pruned.keep",
    output:
        "{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.bed"
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
        hapmap_fam="{dir}/ancestry/{pheno}/{hapmap}.intersection.nomissnp.fa,",
        keep="{dir}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_v3_maf{maf}.pruned.keep",
        flipped="{ukb}/ancestry/{pheno}/compare.{hapmap}.ukb_imp_v3_maf{maf}.pruned.flippedAlleles",
    output:
        "{dir}/ancestry/{pheno}/HapMap{hapmap}.intersection.matched.bed",
    shell:
        "plink --bed {input.hapmap_bed} \
            --bim {input.hapmap_bim} \
            --fam {input.hapmap_fam} \
            --extract {input.keep} \
            --update-alleles {input.flipped} \
            --make-bed \
            --out {wildcards.dir}/ancestry/{wildcards.pheno}/HapMap{wildcards.hapmap}.intersection.matched"

rule mergeHapMap:
    input:
        ukb_bed="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.bed",
        ukb_bim="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.bim",
        ukb_fam="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.fam",
        hapmap_bed="{dir}/ancestry/{pheno}/{hapmap}.intersection.matched.bed",
        hapmap_bim="{dir}/ancestry/{pheno}/{hapmap}.intersection.matched.bim",
        hapmap_fam="{dir}/ancestry/{pheno}/{hapmap}.intersection.matched.fam",
    output:
        "{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.bed"
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
        fam="{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged.fam",
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
        samples="{dir}/ancestry/{pheno}/ukb_imp_genome_v3_maf{maf}.pruned.intersection.matched.fam",
        pcafile="{dir}/ancestry/{pheno}/{hapmap}_ukb_imp_genome_v3_maf{maf}.pruned.merged_pca",
    output:
        "{dir}/ancestry/{pheno}/HapMap_UKBB-FD_pca.png",
        "{dir}/ancestry/{pheno}/European_samples.csv"
    shell:
        "Rscript 'ancestry/selectPCA.R' --directory={wildcards.dir}/ancestry \
            --pcadata={input.pcafile} \
            --samples={input.samples} \
            --name=UKBB-FD"

rule filterEuropeans:
    input:
        keep="{dir}/ancestry/{pheno}/European_samples.csv",
        bed="{dir}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.bed",
        bim="{dir}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.bim",
        fam="{dir}/maf{maf}/ukb_imp_genome_v3_maf{maf}.pruned.fam",
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
                --ndim 10 \
                --outpc {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.pca \
                --outvec  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.eigenvectors \
                --outload  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.SNPloadings \
                --outval  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.eigenvalues \
                --outpve  {wildcards.dir}/ancestry/{wildcards.pheno}/ukb_imp_genome_v3_maf{wildcards.maf}.pruned.European.variance_explained"

