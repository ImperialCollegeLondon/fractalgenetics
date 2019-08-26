# how to run fro ./:
# snakemake -s association_replication.smk --jobs 5000 --latency-wait 30 --cluster-config config/cluster.json --cluster 'bsub -J {cluster.name} -q {cluster.queue} -n {cluster.n} -R {cluster.resources}
# -M {cluster.memory}  -o {cluster.output} -e  {cluster.error}' --keep-going --rerun-incomplete --use-conda

configfile: "config/config_association.yaml"

rule all:
    input:
        expand("{ukb}/gwas/{pheno}/bgenie_{analysis}_lm_st_chr{chr}.gz",
            ukb=config["ukbdir"],
            analysis=['slices'],
            pheno=config['pheno'],
            chr=range(1,23)),
        expand("{ukb}/gwas/{pheno}/bgenie_{analysis}_lm_pseudomt_{plot}.pdf",
            analysis=['slices'],
            plot=['qqplot', 'manhattanplot'],
            pheno=config['pheno'],
            ukb=config["ukbdir"]),
        expand("{ukb}/gwas/{pheno}/Replication_association_significant_in_discovery.txt",
            pheno=config['pheno'],
            ukb=config["ukbdir"]),


rule association:
    input:
        geno=expand("{geno}/ukb_imp_chr{{chr}}_v3.bgen",
            geno=config["genodir"]),
        pheno=expand("{ukb}/phenotypes/{{pheno}}/FD_{{analysis}}_bgenie.txt",
            ukb=config["ukbdir"]),
        covs=expand("{ukb}/phenotypes/{{pheno}}/FD_covariates_bgenie.txt",
            ukb=config["ukbdir"]),
        rsids=expand("{ukb}/maf0.001/{discovery}/ukb_imp_chr{{chr}}_v3_maf0.001.rsid",
                    discovery=config["discovery"],
                    ukb=config["ukbdir"])
    params:
        n=config["n"]
    output:
        "{dir}/gwas/{pheno}/bgenie_{analysis}_lm_st_chr{chr}.gz"
    shell:
        "bgenie --bgen {input.geno} \
            --pheno {input.pheno} \
            --covar {input.covs} \
            --include_rsids  {input.rsids} \
            --pvals --exc_missing_inds \
            --thread {params.n} \
            --out {wildcards.dir}/gwas/{wildcards.pheno}/bgenie_{wildcards.analysis}_lm_st_chr{wildcards.chr}"

rule results:
    input:
        tags=expand("{{dir}}/tags/{discovery}/European_ukb_imp_chr{chr}_v3_maf{maf}_{kb}kb_r{r2}.tags.list",
                    discovery=config["discovery"],
                    maf=config['maf'],
                    r2=config['r2'],
                    kb=config['kbwindow'],
                    chr=range(1,23)),
        gwas=expand("{{dir}}/gwas/{{pheno}}/bgenie_{{analysis}}_lm_st_chr{chr}.gz",
                    chr=range(1,23)),
        pheno="{dir}/phenotypes/{pheno}/FD_{analysis}_EUnorel.csv"
    output:
        "{dir}/gwas/{pheno}/bgenie_{analysis}_lm_pseudomt_qqplot.pdf",
        "{dir}/gwas/{pheno}/bgenie_{analysis}_lm_pseudomt_manhattanplot.pdf",
        "{dir}/gwas/{pheno}/Pseudomultitrait_{analysis}_sig5e08.txt"
    shell:
        "Rscript association/association-results.R \
            --pheno {input.pheno} \
            --name {wildcards.analysis} \
            --directory {wildcards.dir}/gwas/{wildcards.pheno} \
            --showProgress "

rule prep_replication:
    input:
        dsv=expand("{{dir}}/gwas/{discovery}/Pseudomultitrait_slices_sig5e08.txt",
            discovery=config['discovery']),
        dsv_filter=expand("{{dir}}/gwas/{discovery}/Pseudomultitrait_slices_sig5e08_ldFiltered.txt",
            discovery=config['discovery']),
        rep_multi="{dir}/gwas/{pheno}/bgenie_slices_lm_pseudomt_all.csv",
        rep_uni="{dir}/gwas/{pheno}/bgenie_slices_lm_st_genomewide.csv"
    output:
        rep_multi="{dir}/gwas/{pheno}/bgenie_slices_lm_pseudomt_all_replication.csv",
        rep_uni="{dir}/gwas/{pheno}/bgenie_slices_lm_st_genomewide_replication.csv",
        dsv="{dir}/gwas/{pheno}/Pseudomultitrait_discovery_slices_sig5e08_rsID.txt",
        rep_multi_filter="{dir}/gwas/{pheno}/bgenie_slices_lm_pseudomt_all_replication_ldFiltered.csv"
    shell:
        """
        cut -d " " -f 2 {input.dsv} > {output.dsv}
        awk -F, 'FNR==NR {{a[$0]; next}} $2 in a'{output.dsv} {input.rep_multi} > \
            {output.rep_multi}
        head -n 1 {input.rep_uni} > {output.rep_uni}
        awk -F, 'FNR==NR {{a[$0]; next}} $2 in a'{output.dsv} {input.rep_uni} >> \
            {output.rep_uni}
        awk -F, 'FNR==NR {{a[$2]; next}} $2 in a'{input.dsv_filter} \
            {input.rep_multi} > {output.rep_multi_filter}
        """

rule replication:
    input:
        rep_filter="{dir}/gwas/{pheno}/bgenie_slices_lm_pseudomt_all_replication_ldFiltered.csv"
    output:
        rep="{dir}/gwas/{pheno}/Replication_association_significant_in_discovery.txt"
    params:
        nloci=config['nloci']
    shell:
        """
        Rscript association/replication-analysis.R \
            --nloci {params.nloci} \
            --directory {wildcards.dir}/gwas/{wildcards.pheno}
        """
