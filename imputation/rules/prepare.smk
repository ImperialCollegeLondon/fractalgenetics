def createStringPlink(wildcards):
    if wildcards.chr == "X_PAR1":
        plink_str="--chr X --from-bp 60001 --to-bp 2699520"
    elif wildcards.chr == "X_PAR2":
        plink_str="--chr X --from-bp 154931044 --to-bp 999999999"
    elif wildcards.chr == "X":
        plink_str="--chr X --from-bp 2699521 --to-bp 154931043"
    else:
        plink_str="--chr {}".format(wildcards.chr)
    return plink_str

rule splitChromosomes:
    input:
        bim=expand("{qcdir}/combined/{{name}}.bim",
            qcdir=config['genodir']),
        bed=expand("{qcdir}/combined/{{name}}.bed",
            qcdir=config['genodir']),
        fam=expand("{qcdir}/combined/{{name}}.fam",
            qcdir=config['genodir'])
    params:
        stringPlink=lambda wildcards: createStringPlink(wildcards)
    output:
        bim="{dir}/unphased/{name}.chr{chr}.bim",
        bed="{dir}/unphased/{name}.chr{chr}.bed",
        fam="{dir}/unphased/{name}.chr{chr}.fam",
    shell:
        """
        plink --bed {input.bed} \
              --bim {input.bim} \
              --fam {input.fam} \
              {params.stringPlink} \
              --make-bed \
              --out {wildcards.dir}/unphased/{wildcards.name}.chr{wildcards.chr}
        """
