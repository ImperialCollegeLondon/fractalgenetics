#!/usr/local/bin/bash
###############################################
###                                         ###
###     Impute gencall-called genotypes     ###
###     1. step: phasing via shapeit        ###
###     2. Impute via impute2               ###
###                                         ###
###     * ref: combined UK10K/1000Genomes   ###
###     * impuation chromosome-wise         ###
###     * possible QC:                      ###
###         ** impute info criterion        ###
###         ** HWE (per ethnicity)          ###
###         ** MAF/MAC                      ###
###                                         ###
###  by Hannah Meyer                        ###
###############################################

# function to check for occurrence of element in array
ElementIn () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
  return 1
}

echo "Sourcing config file and making pipeline directories"

### defaults
ethnicity_filter="yes"

### source configuration file ###
. $1

### make pipeline directories
mkdir -p $unphaseddir/log
mkdir -p $phasedir/log
mkdir -p $imputedir/log
mkdir -p $resultdir/log
mkdir -p $imputedir/plots


##############################################
### marker exclusion based on clusterplots ###
##############################################

geno_exclude=$geno
files=`ls $imputedir`
if [[ $files == *"-excludeSNPs.txt"* ]]; then
    echo "Marker exclusion based on cluster plot re-evaluation"
    cat $imputedir/*-excludeSNPs.txt > $imputedir/excludeSNPs.list
    echo -e "plink \
        --bfile $rawdir/$geno \
        --exclude $imputedir/excludeSNPs.list \
        --make-bed \
        --out $rawdir/$geno.marker_exclusion
    mv $rawdir/$geno.marker_exclusion.log \
        $rawdir/plink_log/$geno.marker_exclusion_before_imputation.log" \
        > $imputedir/log/marker_exclusion.cmd
   # bsub -g /impute  -J $geno.marker_exclusion_$iteration -o $imputedir/log/marker_exclusion_$iteration.log -e $imputedir/log/marker_exclusion_$iteration.err < $imputedir/log/marker_exclusion.cmd
    geno_exclude=$geno.marker_exclusion
fi


##############################
### phasing and imputation ###
##############################

## phasing and imputation done in a per-chromosome manner
#for chr in 22; do
for chr in `seq 1 22; echo 'X X_PAR1 X_PAR2'`; do

    # set initial id and formating flags
    id="UK10K1000Genomes"
    format_flag=0

    # imputation info score filtering
    if [[ $imputeThr != "" ]]; then
        imputeStr=imputeThr$imputeThr.
    else
        imputeStr=""
    fi

    # set X chromosome boundaries split up in PAR1, PAR2 and rest of chrX
    if [[ "$chr" == "X_PAR1" ]]; then
        pos=24
        plink_str="--chr X --from-bp 60001 --to-bp 2699520"
        chrX_phase_str=""
        chrX_impute_str="-chrX -Xpar"
    elif [[ "$chr" == "X" ]]; then
        pos=23
        plink_str="--chr X --from-bp 2699521 --to-bp 154931043"
        chrX_phase_str="--chrX"
        chrX_impute_str="-chrX"
    elif [[ "$chr" == "X_PAR2" ]]; then
        pos=25
        plink_str="--chr X --from-bp 154931044 --to-bp 999999999"
        chrX_phase_str=""
        chrX_impute_str="-chrX -Xpar"
    else
        pos=$((chr-1))
        plink_str="--chr $chr"
        chrX_phase_str=""
        chrX_impute_str=""
    fi


    ##########################################################################
    ### step 1: pre-phase using SHAPEIT (Delaneau et al (2013) Nat Methods ###
    ##########################################################################

    if [[ ! -e $phasedir/chr$chr.hap.gz ]]; then

        # keep track of which chromosomes are being phased
        phasing[$pos]=chr$chr 
        
        # create phasing command files: split plink files per chromosome; 
        # phase these files via shapeit
        echo -e "#! /usr/local/bin/bash
        \nplink --noweb \
            --bfile $rawdir/$geno_exclude  $plink_str \
            --make-bed \
            --out $unphaseddir/chr$chr
        \nshapeit --thread $thread \
            --window $window_size \
            --states $states \
            --effective-size 11418 -B $unphaseddir/chr$chr \
            --input-map $refdir/genetic_map_chr${chr}_combined_b37.txt \
            --output-log $phasedir/chr$chr.shapeit \
            --output-max $phasedir/chr$chr.hap.gz \
            $phasedir/chr$chr.sample $chrX_phase_str
        " > $phasedir/log/chr$chr.cmd

        # submit phasing jobs per chromosome
        echo "Submitted phasing chromosome $chr"
        if [[ -e $imputedir/excludeSNPs.list ]] && \
            [[ ! -e $rawdir/${geno}.bed ]]; then
            bsub -J $geno.shapeit.chr$chr \
                -w "ended($geno.marker_exclusion_$iteration)" \ 
            -o $phasedir/log/chr$chr.shapeit.log \
                -e $phasedir/log/chr$chr.shapeit.err \
                -n $thread \
                -R "span[ptile=$thread] select[mem>5000] rusage[mem=5000]" \ 
            -M 5000 < $phasedir/log/chr$chr.cmd >/dev/null
        else
            bsub -J $geno.shapeit.chr$chr \
                -o $phasedir/log/chr$chr.shapeit.log \
                -e $phasedir/log/chr$chr.shapeit.err \
                -n $thread \
                -R "span[ptile=$thread] select[mem>5000] rusage[mem=5000]" \
                -M 5000 < $phasedir/log/chr$chr.cmd >/dev/null
        fi
    fi

    #############################################################
    ### step 2: impute via IMPUTE2  on chunks  per chromosome ###
    #############################################################

    # make directory to store imputation per chunk results
    mkdir -p $imputedir/chr$chr/log

    echo "Imputation based on chunk boundaries of imputation set"
    chunk_num=`ls $refdir | grep -e "chr$chr\..*.shapeit.hap.gz" | tail -n 1 | sed  's/chr[0-9X]*\.\([0-9X]*\)\.shapeit\.hap.gz/\1/'`
    if [[ $chr == "X_PAR1" || $chr == "X_PAR2" ]]; then
        chunk_num=1
    fi
    
    # submit chunk-wise imputation jobs
    echo "Submitting $chunk_num imputation jobs for chromosome $chr"
    for chunk in `seq 1 $chunk_num`; do

        # create chunkStr for file naming
        chunkStr=`printf "%02d" $chunk`

        # if (any) chunk has to be imputed, set format_flag to 1, 
        # so formatting jobs are submitted; else jump to next chunk
        if [[ -e $imputedir/chr$chr/chr$chr.$chunkStr.gen.gz ]]; then
            continue
        fi
        format_flag=1

        # use reference chunk boundaries
        mem=30000
        if [[ $chr == "X_PAR1" ]]; then
            refhap=$refdir/chr$chr.shapeit.hap.gz
            reflegend=$refdir/chr$chr.shapeit.legend.gz
            chunk_begin=60001
            chunk_end=2699520
        elif [[ $chr == "X_PAR2" ]]; then
            refhap=$refdir/chr$chr.shapeit.hap.gz
            reflegend=$refdir/chr$chr.shapeit.legend.gz
            chunk_begin=154931044
            chunk_end=999999999
        else
            refhap=$refdir/chr$chr.${chunkStr}.shapeit.hap.gz
            reflegend=$refdir/chr$chr.${chunkStr}.shapeit.legend.gz
            chunk_begin=`cat $refdir/chr$chr.$chunkStr.shapeit.cmd | grep '\-int' |sed 's/.*-int \([0-9]*\) [0-9]*.*/\1/'`
            chunk_end=`cat $refdir/chr$chr.$chunkStr.shapeit.cmd | grep '\-int' |sed 's/.*-int [0-9]* \([0-9]*\).*/\1/'`
        fi
        
        # send imputation jobs to job group /impute, 
        # remove potential old summary and log-files first
        rm -f $imputedir/chr$chr/chr$chr.$chunkStr.gen_summary
        rm -f $imputedir/chr$chr/chr{$chr}_summary_overview
        rm -f $imputedir/chr$chr/chr${chr}_err_overview
        rm -f $imputedir/chr$chr/chr${chr}_warnings_overview

        # construct imputation command based on sample exclusion file
        if [[ ! -e $exclude_samples ]]; then
            echo -e "impute2 -allow_large_regions \
                -m $refdir/genetic_map_chr${chr}_combined_b37.txt \
                -h $refhap -l $reflegend -use_prephased_g \
                -known_haps_g $phasedir/chr$chr.hap.gz \
                -sample_g $phasedir/chr$chr.sample $extra_str \
                -k_hap $k_hap -int $chunk_begin $chunk_end \
                -Ne 20000 -buffer $buffer_size -verbose \
                -o $imputedir/chr$chr/chr$chr.$chunkStr.gen $chrX_impute_str" \
                > $imputedir/chr$chr/log/chr$chr.$chunkStr.cmd
        else
            echo -e "impute2 -allow_large_regions \
                -m $refdir/genetic_map_chr${chr}_combined_b37.txt \
                -h $refhap -l $reflegend -use_prephased_g \
                -known_haps_g $phasedir/chr$chr.hap.gz \
                -sample_g $phasedir/chr$chr.sample $extra_str \
                -exclude_samples_g $exclude_samples -k_hap $k_hap \
                -int $chunk_begin $chunk_end -Ne 20000 \
                -buffer $buffer_size -verbose \
                -o $imputedir/chr$chr/chr$chr.$chunkStr.gen $chrX_impute_str" \
                > $imputedir/chr$chr/log/chr$chr.$chunkStr.cmd
        fi 
        
        # Imputation depends on phasing being finished...
        # check if phasing jobs finished and imputation file doesn't exist yet
        if [[ ! -e $imputedir/chr$chr/chr$chr.$chunkStr.gen.gz ]] \
            && ! ElementIn "chr$chr" "${phasing[@]}"; then
                bsub -g /impute -J $geno.chr$chr.$chunkStr \
                    -R "select[mem>$mem] rusage[mem=$mem]" \
                    -M${mem} -o $imputedir/chr$chr/log/chr$chr.$chunkStr.log \
                    -e $imputedir/chr$chr/log/chr$chr.$chunkStr.err \
                    < $imputedir/chr$chr/log/chr$chr.$chunkStr.cmd  >/dev/null 
        else
                bsub -g /impute -J $geno.chr$chr.$chunkStr \
                    -w "ended($geno.shapeit.chr$chr)" \
                    -R "select[mem>$mem] rusage[mem=$mem]" -M${mem} \
                    -o $imputedir/chr$chr/log/chr$chr.$chunkStr.log \
                    -e $imputedir/chr$chr/log/chr$chr.$chunkStr.err \
                    < $imputedir/chr$chr/log/chr$chr.$chunkStr.cmd  >/dev/null 
        fi

        # examine sample and snp info files for missing SNPs and compress the chunk-output files
        echo -e "#! /usr/local/bin/bash \
        \nif grep -q 'There are no SNPs in the imputation interval' $imputedir/chr$chr/chr$chr.$chunkStr.gen_summary; then \
            \necho \"chr$chr\tchunk $chunkStr\tThere are no SNPS in the imputation interval $chunk_begin $chunk_end of chromosome $chr, corresponding to chunk number $chunkStr\n\" >> $imputedir/chr$chr/chr${chr}_summary_overview;
        \nelif grep -q 'ERROR' $imputedir/chr$chr/chr$chr.$chunkStr.gen_summary; then
            \necho \"chr$chr chunk $chunkStr\" >> $imputedir/chr$chr/chr${chr}_summary_overview;
            \necho \`grep 'ERROR' $imputedir/chr$chr/chr$chr.$chunkStr.gen_summary\` >> $imputedir/chr$chr/chr${chr}_summary_overview;
        \nelse
            \ngzip -f $imputedir/chr$chr/chr$chr.$chunkStr.gen
            \nN_info=\`awk 'NR>1' $imputedir/chr$chr/chr$chr.$chunkStr.gen_info | wc -l | awk '{printf \$1}'\`
            \nN_gen=\`zcat $imputedir/chr$chr/chr$chr.$chunkStr.gen.gz | wc -l | awk '{printf \$1}'\`
        \nfi
        \necho \"$chr\t$chunkStr\t\" >> $imputedir/chr$chr/chr${chr}_err_overview
        \ncat $imputedir/chr$chr/chr$chr.$chunkStr.err >> $imputedir/chr$chr/chr${chr}_err_overview
        \necho \"$chr\t$chunkStr\t\" >> $imputedir/chr$chr/chr${chr}_warnings_overview
        \ncat $imputedir/chr$chr/chr$chr.$chunkStr.gen_warnings >> $imputedir/chr$chr/chr${chr}_warnings_overview
        " > $imputedir/chr$chr/log/chr$chr.${chunkStr}_format.cmd
        bsub -g /impute -J $geno.chr$chr.${chunkStr}_format \
            -w "ended($geno.chr$chr.$chunkStr)" \
            -o $imputedir/chr$chr/log/chr$chr.$chunkStr.log \
            -e $imputedir/chr$chr/log/chr$chr.$chunkStr.err < \
            $imputedir/chr$chr/log/chr$chr.${chunkStr}_format.cmd >/dev/null

        #########################################################
        ### imputationQC: filter SNPs based on info-criterion ###
        #########################################################
        
        if [[ $imputeThr != "" ]]; then
            imputeStr=imputeThr$imputeThr.
            echo -e "#! /usr/local/bin/bash
            \n awk -v imputeThr=$imputeThr '{if(\$5 <= $imputeThr)  print NR-1}' $imputedir/chr$chr/chr${chr}.${chunkStr}.gen_info  > $imputedir/chr$chr/chr${chr}.${chunkStr}.gen_imputeqc_fail
            \n perl -e 'open(IN, \$ARGV[0]); while(<IN>) {chomp \$_; \$line{\$_} =1}; close(IN); open(IN, \"gunzip -c \$ARGV[1] |\"); while(<IN>) {unless (exists(\$line{\$.})) {print \$_}}; close(IN)' $imputedir/chr$chr/chr${chr}.${chunkStr}.gen_imputeqc_fail $imputedir/chr$chr/chr${chr}.${chunkStr}.gen.gz > $imputedir/chr$chr/chr${chr}.${chunkStr}.${imputeStr}gen
            \n gzip -f $imputedir/chr$chr/chr$chr.$chunkStr.${imputeStr}gen"  > $imputedir/chr$chr/log/chr$chr.${chunkStr}_imputeThr.cmd
            bsub -g /impute -J $geno.chr$chr.${chunkStr}_imputeThr \
                -w "ended($geno.chr$chr.${chunkStr}_format)" \
                -o $imputedir/chr$chr/log/chr$chr.$chunkStr.imputeThr.log \
                -e $imputedir/chr$chr/log/chr$chr.$chunkStr.imputeThr.err \
                < $imputedir/chr$chr/log/chr$chr.${chunkStr}_imputeThr.cmd \
                >/dev/null
        fi
    done

    # concatenate all chunk-wide output files per chromosome into one output file per chromosome
    if [[ $format_flag -eq 1 ]]; then
        if [[ $imputeThr != "" ]]; then
            bsub -g /impute -J "$geno.chr${chr}.gz" \
                -o $resultdir/log/chr${chr}.concatenate.log \
                -e $resultdir/log/chr${chr}.concatenate.err \
                -w "ended($geno.chr$chr.*_imputeThr)" \
                "zcat $imputedir/chr$chr/chr${chr}.*.${imputeStr}gen.gz > \
                $resultdir/chr${chr}.${imputeStr}gen.$id" >/dev/null
        else
            bsub -g /impute -J "$geno.chr${chr}.gz" \
                -o $resultdir/log/chr${chr}.concatenate.log \
                -e $resultdir/log/chr${chr}.concatenate.err \
                -w "ended($geno.chr$chr.*_format)" \ 
                "zcat $imputedir/chr$chr/chr${chr}.*.${imputeStr}gen.gz > \
                    $resultdir/chr${chr}.${imputeStr}gen.$id" >/dev/null
        fi
        
        ##############################
        ### filter for MAF and HWE ### 
        ##############################
        
        if [[ $postqc_filter == "yes" ]]; then
            # if ethnicity filter, then do simple MAF/HWE filtering
            if [[ $ethnicity_filter == "yes" ]]; then
                echo -e "#! /usr/local/bin/bash
                if [[ ! -e $imputedir/snptest_samples.txt ]]; then
                    #awk '{if (NR==2){\$7=\"P\"} print \$0}'  $imputedir/chr$chr/chr$chr.10.gen_samples > $imputedir/snptest_samples.txt
                fi
                #snptest -summary_stats_only -hwe -data $resultdir/chr${chr}.${imputeStr}gen.$id.gz $imputedir/snptest_samples.txt -o $resultdir/chr${chr}.gen.stats
                #awk  -F' ' '!/^#/ &&  (\$19 < $MAFThr || \$21 < $HWEThr) {print NR-10,\$2,\$4, \$19,\$21}' $resultdir/chr${chr}.gen.stats > $resultdir/chr${chr}.gen.toFilter
                awk 'FNR==NR{a[\$2\$3]++; next} !a[\$2\$3]'  $resultdir/chr${chr}.gen.toFilter $resultdir/chr${chr}.${imputeStr}gen.$id > $resultdir/chr$chr.${imputeStr}gen.snptest
                gzip -f $resultdir/chr$chr.${imputeStr}gen.snptest"  > $resultdir/log/chr${chr}_postimpute.cmd
                
                bsub -g /impute -J $geno.chr$chr.postimpute -w "ended($geno.chr${chr}.gz)" -o $resultdir/log/chr$chr.postimpute.log -e $resultdir/log/chr$chr.postimpute.err < $resultdir/log/chr${chr}_postimpute.cmd >/dev/null
            else
                # else, split cohort based on ethnicities and 
                # do separate HWE filtering
                echo -e "#! /usr/local/bin/bash
                if [[ ! -e $imputedir/snptest_samples.txt ]]; then
                    awk '{if (NR==2){\$7=\"P\"} print \$0}'  $imputedir/chr$chr/chr$chr.10.gen_samples > $imputedir/snptest_samples.txt
                fi
                for file in caucasian african mexican; do
                    snptest -summary_stats_only -hwe -data $resultdir/chr${chr}.${imputeStr}gen.$id.gz $imputedir/snptest_samples.txt -exclude_samples $rawdir/gencall.\${file}_exclude.txt -o $resultdir/chr${chr}.\${file}.gen.stats
                    awk  -F' ' '!/^#/ &&  \$21 < $HWEThr {print NR-10,\$2,\$4,\$19,\$21}' $resultdir/chr${chr}.\${file}.gen.stats > $resultdir/chr${chr}.\${file}.gen.toFilter
                    cat $resultdir/chr${chr}.\${file}.gen.toFilter >> $resultdir/chr${chr}.gen.toFilter
                done
                snptest -summary_stats_only -hwe -data $resultdir/chr${chr}.${imputeStr}gen.$id.gz $imputedir/snptest_samples.txt  -o $resultdir/chr${chr}.alleth.gen.stats
                awk  -F' ' '!/^#/ &&  \$19 < $MAFThr  {print NR-10,\$2,\$4,\$19,\$21}' $resultdir/chr${chr}.alleth.gen.stats > $resultdir/chr${chr}.alleth.gen.toFilter
                cat $resultdir/chr${chr}.alleth.gen.toFilter >> $resultdir/chr${chr}.gen.toFilter
                awk 'FNR==NR{a[\$2\$3]++; next} !a[\$2\$3]'  $resultdir/chr${chr}.gen.toFilter <(gunzip -c $resultdir/chr$chr.${imputeStr}gen.$id.gz) > $resultdir/chr$chr.${imputeStr}gen.snptest
                gzip -f $resultdir/chr$chr.${imputeStr}gen.snptest"  > $resultdir/log/chr${chr}_postimpute.cmd
            
                bsub -g /impute -J $geno.chr$chr.postimpute -w "ended($geno.chr${chr}.gz)" -o $resultdir/log/chr$chr.postimpute.log -e $resultdir/log/chr$chr.postimpute.err < $resultdir/log/chr${chr}_postimpute.cmd >/dev/null
            fi
            id="snptest"
        fi

        #############################
        ### reformat output files ###
        #############################
        
        echo -e "
        \nzcat $resultdir/chr$chr.${imputeStr}gen.$id.gz | awk -v chr=$chr '{printf chr; \$1=\"\"; print \$0}' | gzip -f   > $resultdir/chr${chr}.${imputeStr}gen.${id}.format
        \nzcat $resultdir/chr$chr.${imputeStr}gen.$id.gz | sed 's/,rs/|rs/g' | awk -v chr=$chr '{ snp=(NF-5)/3; if(\$2 ~/^rs/) s=\$2;else s=\"NA\"; printf \"chr\"chr\":\"\$3\"-\"s\"-\"\$4\"-\"\$5\",\" \$4 \",\" \$5; for(i=1; i<=snp; i++) printf \",\" \$(i*3+3)*2+\$(i*3+4); printf \"\\\n\" }' | gzip -f   > $resultdir/chr${chr}.impute2.UK10K1000Genomes.${imputeStr}$id.gz
        \nzcat $resultdir/chr$chr.${imputeStr}gen.$id.gz | awk -v chr=$chr '{if(\$2 ~/^rs/) s=\$2;else s=\"NA\"; print \"chr\"chr\":\"\$3\"-\"s\"-\"\$4\"-\"\$5\",\" \$3 \",\" chr}' > $resultdir/chr${chr}.impute2.UK10K1000Genomes.${imputeStr}$id.pos " > $resultdir/log/chr$chr.$id.cmd

        if [[ $postqc_filter == "yes" ]]; then
            bsub -J $geno.step2.chr$chr -w "ended($geno.chr${chr}.postimpute)" \
                -o $resultdir/log/chr$chr.$id.LOG \
                -e $resultdir/log/chr$chr.$id.ERR \
                -R "select[mem>300] rusage[mem=300]" \
                -M 300 < $resultdir/log/chr$chr.$id.cmd >/dev/null
        else
            bsub -J $geno.step2.chr$chr -w "ended($geno.chr${chr}.gz)" \
                -o $resultdir/log/chr$chr.$id.LOG \
                -e $resultdir/log/chr$chr.$id.ERR \
                -R "select[mem>300] rusage[mem=300]" \
                -M 300 < $resultdir/log/chr$chr.$id.cmd >/dev/null
        fi
    fi
done

##################################################################
### combine chromsome (autosomes) files into genome file (exX) ###
##################################################################

echo -e "#! /usr/local/bin/bash
if [[ -e $resultdir/genomeExX.impute2.UK10K1000Genomes.${imputeStr}$id.gz ]]; \ 
    then
    rm $resultdir/genomeExX.impute2.UK10K1000Genomes.${imputeStr}$id.gz
    rm $resultdir/genomeExX.impute2.UK10K1000Genomes.${imputeStr}$id.pos
fi
for chr in \`seq 1 22\`; do
    cat $resultdir/chr\${chr}.impute2.UK10K1000Genomes.${imputeStr}$id.gz >> \
        $resultdir/genomeExX.impute2.UK10K1000Genomes.${imputeStr}$id.gz
    cat $resultdir/chr\${chr}.impute2.UK10K1000Genomes.${imputeStr}$id.pos >> \
        $resultdir/genomeExX.impute2.UK10K1000Genomes.${imputeStr}$id.pos
done" > $resultdir/log/combine_impute2.UK10K1000Genomes.${imputeStr}${id}.cmd
bsub -J $geno.combine \
    -w "ended($geno.step2.chr*)" \
    -o $resultdir/log/combine_impute2.UK10K1000Genomes.${imputeStr}${id}.log \
    -e $resultdir/log/combine_impute2.UK10K1000Genomes.${imputeStr}${id}.err \
    -R "select[mem>300] rusage[mem=300]" -M 300 \
    < $resultdir/log/combine_impute2.UK10K1000Genomes.${imputeStr}${id}.cmd \
    >/dev/null

###########################
### create files for QC ###
###########################

bsub -J $geno.imputeQC \
    -w "ended($geno.combine)" \
    -o $resultdir/log/imputeQC.log \
    -e $resultdir/log/imputeQC.err \
    -R "select[mem>300] rusage[mem=300]" -M 300 \
    "bash ~/GWAS/analysis/genotyping/imputeControl.sh \
            $resultdir/imputeQC_config"
