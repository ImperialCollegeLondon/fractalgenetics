###############################################################
###############################################################
###                                                         ###
###     Script for:                                         ###
###         * analysing differential missingness rates      ###
###           between batches                               ###
###                                                         ###
###     by Hannah Meyer                                     ###
###                                                         ###
###############################################################
###############################################################

#############################################
### 1. analysis parameters and input data ###
#############################################

### a) source and display input parameters:
dir=$1
batchnames=$2
alg=$3
lmissTh=$4
commonprobes=$5

echo "dir (Output directory): $dir"
echo "batchnames (genotyping batches to be combined): $batchnames"
echo "alg (genotype calling algorithm): $alg"
echo "lmissTh (threshold for exclusion based on missing genotype rate): $alg"
echo "commonprobes (file with rsID of common probes across batches): $commonprobes"

### b) Copy genotype files from respective batch directory ###
batches=($batchnames)
combined="$dir/combined"

for b in "${!batches[@]}"; do
    batch=${batches[$b]}
    batchdir="$dir/$batch"
    cp $batchdir/$alg.${batch}.bed $combined/$alg.${batch}.raw.bed
    cp $batchdir/$alg.${batch}.bim $combined/$alg.${batch}.raw.bim
    cp $batchdir/$alg.${batch}.fam $combined/$alg.${batch}.raw.fam
    batches[$b]=$alg.${batches[$b]}
done

######################################################
### Differential missingness rates between batches ###
######################################################

### a) Filter SNPs
## i) extract SNPs with same probe ID
for file in "${batches[@]}"; do
    plink --bfile $combined/$file.raw --extract $commonprobes \
        --make-bed --out $combined/$file.commonProbes
    mv $combined/${file}.commonProbes.log $combined/log/${file}.keepcommonProbes.log
done

## ii) Filter for SNPs present in all batches
awk '{a[$2]++} a[$2] == 3 {print $2}' $combined/${batches[0]}.commonProbes.bim \
    $combined/${batches[1]}.commonProbes.bim \
    $combined/${batches[2]}.commonProbes.bim \
    > $combined/SNPs_in_all_batches

for file in "${batches[@]}"; do
    plink --bfile $combined/$file.commonProbes \
        --extract $combined/SNPs_in_all_batches \
        --make-bed --out $combined/$file
    mv $combined/${file}.log $combined/log/${file}.keepcommonSNPs.log
done

### b)  Analyse SNPs for differential missing genotypes (batch effects)
union12=${batches[0]}.${batches[1]}.union
union13=${batches[0]}.${batches[2]}.union
union23=${batches[1]}.${batches[2]}.union

# i) create pair-wise merge files
plink --bmerge $combined/${batches[1]} --bfile $combined/${batches[0]} \
    --make-bed --out $combined/$union12
mv $combined/$union12.log $combined/log/$union12.log
plink --bmerge $combined/${batches[2]} --bfile $combined/${batches[0]} \
    --make-bed --out $combined/$union13
mv $combined/$union13.log $combined/log/$union13.log
plink --bmerge $combined/${batches[2]} --bfile $combined/${batches[1]} \
    --make-bed --out $combined/$union23
mv $combined/$union23.log $combined/log/$union23.log

# ii) create case files
cut -d " " -f 1,2 $combined/${batches[0]}.fam > $combined/cases.${batches[0]}
cut -d " " -f 1,2 $combined/${batches[1]}.fam > $combined/cases.${batches[1]}

# iii) create case/control setting (one batch case, other control)
# use plink's test-missing to find genos with differential missingness rate
for file in $union12 $union13 $union23; do
    if [[ $file == *"${batches[0]}"* ]]; then
        cases=cases.${batches[0]}
    else
        cases=cases.${batches[1]}
    fi

    plink --bfile $combined/$file --make-pheno $combined/$cases '*' \
        --geno $lmissTh --make-bed --out $combined/$file.pheno
    mv $combined/$file.pheno.log $combined/log/$file.createpheno.log

    # find differential missingness rate with p < 1e-5
    plink --bfile $combined/$file.pheno --test-missing midp --pfilter 0.00001 \
        --out $combined/$file
    mv $combined/$file.log $combined/log/$file.differentialmissingsness.log

done

# iii) create union of differentially missing SNPs and remove from datasets
awk '{print $2}'  $combined/$union12.missing \
    $combined/$union13.missing \
    $combined/$union23.missing |sort |uniq | tail -n +2 \
    > $combined/differentialmissing.SNP

# iv) Remove SNPs with differential genotype missingness rate
for file in "${batches[@]}" ; do
    plink --bfile $combined/$file --exclude $combined/differentialmissing.SNP \
        --make-bed --out $combined/$file
    mv $combined/${file}.log $combined/log/${file}.removemissnp.log
done

### c) Merge cleaned datasets into one common set
# i) write list of files to merge
echo -e "$combined/${batches[0]}\n$combined/${batches[1]}\n$combined/${batches[2]}" \
    > $combined/merge_list
# ii) merge files
plink --merge-list $combined/merge_list  --make-bed \
    --out $combined/$alg.combined
mv $combined/$alg.combined.log $combined/log/$alg.combined.mergelist.log
