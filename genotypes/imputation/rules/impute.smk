def createStringImpute(wildcards):
    if wildcards.chr == "X_PAR1":
        impute_str="-chrX -Xpar"
    elif wildcards.chr == "X_PAR2":
        impute_str="-chrX -Xpar"
    elif wildcards.chr == "X":
        impute_str="-chrX"
    else:
        impute_str=""
    return impute_str

def chunkSettings(wildcards):
    with open ("{}/chunkBoundariesChr{}.txt".format(config['referencedir'],
        wildcards.chr)) as file:
        content = file.readlines()
        content = [x.strip() for x in content]
        for line in content:
            parts = line.split()
            if parts[1] == str(wildcards.chunk):
                if len(parts) > 2:
                    return "{} {}".format(parts[2], parts[3])
        return False

rule imputation:
    input:
        hap="{dir}/phased/{name}.chr{chr}.hap.gz",
        sample="{dir}/phased/{name}.chr{chr}.sample"
    output:
        gen="{dir}/imputed/chr{chr}/{name}.chr{chr}.{chunk}.gen",
        warnings="{dir}/imputed/chr{chr}/{name}.chr{chr}.{chunk}.gen_warnings",
        summary="{dir}/imputed/chr{chr}/{name}.chr{chr}.{chunk}.gen_summary",
        samples="{dir}/imputed/chr{chr}/{name}.chr{chr}.{chunk}.gen_samples",
        info="{dir}/imputed/chr{chr}/{name}.chr{chr}.{chunk}.gen_info"
    params:
        buffer=config['buffer_size'],
        k_hap=config['k_hap'],
        states=config['states'],
        Ne=config['Ne'],
        ref=config['referencedir'],
        map=lambda wildcards: "{}/genetic_map_chr{}_combined_b37.txt".format(
            config['referencedir'], wildcards.chr),
        reflegend=lambda wildcards: "{}/chr{}.{:02d}.shapeit.legend.gz".format(
            config['referencedir'], wildcards.chr, int(wildcards.chunk)),
        refhap=lambda wildcards: "{}/chr{}.{:02d}.shapeit.hap.gz".format(
            config['referencedir'], wildcards.chr, int(wildcards.chunk)),
        boundaries=lambda wildcards: chunkSettings(wildcards),
        stringImpute=lambda wildcards: createStringImpute(wildcards),
    shell:
         """
         impute2 -allow_large_regions \
                    -m {params.map} \
                    -h {params.refhap} -l {params.reflegend} -use_prephased_g \
                    -known_haps_g {input.hap} \
                    -sample_g {input.sample} \
                    -k_hap {params.k_hap} \
                    -Ne {params.Ne} -buffer {params.buffer} -verbose \
                    -int {params.boundaries} \
                    -o {output.gen} {params.stringImpute}
         """
