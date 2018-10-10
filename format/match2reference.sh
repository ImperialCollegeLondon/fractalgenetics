####################################
### 1. input parameters and data ###
####################################

alg=$1
qcdir=$2
UK10K1KGdir=$3
ensemblNotInUK10K1KG=$4

mkdir -p $qcdir/plink_log

##################################################
### 2. match input files to reference datasets ###
##################################################

# reference data sets: imputation (UK10K1KG) and ENSEMBL
# (if not present in UK10K1000Genomes)
# (ENSEMBL generated via format/format_20160415_GRCh37p13_human_vcf.sh)

### a) remove variants that are neither ensembl nor UK10K1KG annotated
# (important for later fusion of genotypes from different batches)
awk 'FNR==NR {a[$2]++; next} a[$2] ' $qcdir/$alg.raw.bim  \
    $UK10K1KGdir/genome.legend.format > $qcdir/${alg}_genome.legend.format
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]++; next} a[$3]  {print $1,$3,0,$2,$4,$5}' \
    $qcdir/$alg.raw.bim  $ensemblNotInUK10K1KG \
    > $qcdir/${alg}_GRCh37.p13_not_in_UK10K1KG.txt

awk 'FNR==NR {a[$2]++; next} !a[$2]' $UK10K1KGdir/genome.legend.format \
    $qcdir/$alg.raw.bim  > $qcdir/${alg}_not_in_genome.legend.format
awk 'FNR==NR {a[$3]++; next} !a[$2] {print $2} ' $ensemblNotInUK10K1KG \
    $qcdir/${alg}_not_in_genome.legend.format  \
    > $qcdir/${alg}_not_in_GRCh37p13_or_UK10K1KG

plink --bfile $qcdir/$alg.raw \
      --exclude $qcdir/${alg}_not_in_GRCh37p13_or_UK10K1KG \
      --make-bed --out $qcdir/$alg.onlyannotated
mv $qcdir/${alg}.onlyannotated.log $qcdir/plink_log

### b) check if SNPs are concordant with UK10K1KG information
# perfect match between SNPs in batch and UK10K1KG
awk 'FNR==NR {a[$1$2$4$5$6]; next} $1$2$4$5$6 in a {print $2}' \
    $qcdir/$alg.onlyannotated.bim $qcdir/${alg}_genome.legend.format > \
    $qcdir/${alg}.match_UK10K1KG

# chromosome mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/$alg.onlyannotated.bim $qcdir/${alg}_genome.legend.format > \
    $qcdir/${alg}.toUpdateChr_UK10K1KG

plink --bfile $qcdir/$alg.onlyannotated \
    --update-chr $qcdir/$alg.toUpdateChr_UK10K1KG 1 2\
      --make-bed \
      --out $qcdir/$alg.updateChr_UK10K1KG
mv $qcdir/$alg.updateChr_UK10K1KG.log \
    $qcdir/plink_log/$alg.updateChr_UK10K1KG.log

# position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $qcdir/$alg.updateChr_UK10K1KG.bim $qcdir/${alg}_genome.legend.format > \
    $qcdir/$alg.toUpdatePos_UK10K1KG

plink --bfile $qcdir/$alg.updateChr_UK10K1KG \
      --update-map $qcdir/$alg.toUpdatePos_UK10K1KG 1 2\
      --make-bed \
      --out $qcdir/$alg.updatePos_UK10K1KG
mv $qcdir/$alg.updatePos_UK10K1KG.log \
    $qcdir/plink_log/$alg.updatePos_UK10K1KG.log

# allele order mismatch
#awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
#    ($1$2$4 in a && a[$1$2$4] == $6$5)  {print $2,$5,$6,$6,$5}' \
#    $qcdir/$alg.updatePos_UK10K1KG.bim $qcdir/${alg}_genome.legend.format > \
#    $qcdir/${alg}.toUpdateAlleles_UK10K1KG

# possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$alg.updatePos_UK10K1KG.bim $qcdir/${alg}_genome.legend.format > \
    $qcdir/${alg}.toFlip_UK10K1KG

# Update and flip alleles
#plink --bfile $qcdir/${alg}.updatePos_UK10K1KG \
#      --update-alleles $qcdir/${alg}.toUpdateAlleles_UK10K1KG \
plink --bfile $qcdir/${alg}.updatePos_UK10K1KG \
      --flip $qcdir/${alg}.toFlip_UK10K1KG \
      --make-bed \
      --out $qcdir/$alg.UK10K1KG_flipped
mv $qcdir/${alg}.UK10K1KG_flipped.log \
    $qcdir/plink_log/${alg}.UK10K1KG_flipped.log

# allele order mismatch of flipped alleles
#awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
#    ($1$2$4 in a && a[$1$2$4] == $6$5)  {print $2,$5,$6,$6,$5}' \
#    $qcdir/$alg.UK10K1KG_flipped.bim $qcdir/${alg}_genome.legend.format > \
#    $qcdir/${alg}.toUpdate_flipped_UK10K1KG

#plink --bfile $qcdir/${alg}.UK10K1KG_flipped \
#      --update-alleles $qcdir/${alg}.toUpdate_flipped_UK10K1KG \
#      --make-bed \
#      --out $qcdir/$alg.UK10K1KG_flipped.updated
#mv $qcdir/${alg}.UK10K1KG_flipped.updated.log \
#    $qcdir/plink_log/${alg}.UK10K1KG_flipped.updated.log

# non-matching alleles
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$alg.UK10K1KG_flipped.bim \
    #$qcdir/$alg.UK10K1KG_flipped.updated.bim \
    $qcdir/${alg}_genome.legend.format > \
    $qcdir/${alg}.toRemove_UK10K1KG

#plink --bfile $qcdir/${alg}.UK10K1KG_flipped.updated \
plink --bfile $qcdir/${alg}.UK10K1KG_flipped \
      --exclude $qcdir/${alg}.toRemove_UK10K1KG \
      --make-bed \
      --out $qcdir/$alg.UK10K1KG
mv $qcdir/${alg}.UK10K1KG.log $qcdir/plink_log/${alg}.UK10K1KG.log


### c) check if SNPs are concordant with ENSEMBL reference info
# check if there are SNPs which not in UK10K1KG reference set

# chromosome mismatch (ignore sex chromosomes)
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$2]=$1; next}} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/$alg.UK10K1KG.bim $qcdir/${alg}_GRCh37.p13_not_in_UK10K1KG.txt |\
    sed -n '/^[XY]/!p' | sort |uniq > $qcdir/${alg}.UK10K1KG.toUpdateChr_ensembl

plink --bfile $qcdir/${alg}.UK10K1KG \
      --update-chr $qcdir/${alg}.UK10K1KG.toUpdateChr_ensembl 1 2 \
      --make-bed \
      --out $qcdir/$alg.UK10K1KG.updateChr_ensembl
mv $qcdir/$alg.UK10K1KG.updateChr_ensembl.log \
    $qcdir/plink_log/$alg.UK10K1KG.updateChr_ensembl.log

# position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$2]=$4; next}} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $qcdir/$alg.UK10K1KG.updateChr_ensembl.bim \
    $qcdir/${alg}_genome.legend.format > \
    $qcdir/$alg.UK10K1KG.toUpdatePos_ensembl

plink --bfile $qcdir/$alg.UK10K1KG.updateChr_ensembl \
      --update-map $qcdir/$alg.UK10K1KG.toUpdatePos_ensembl 1 2\
      --make-bed \
      --out $qcdir/$alg.UK10K1KG.updatePos_ensembl
mv $qcdir/$alg.UK10K1KG.updatePos_ensembl.log \
    $qcdir/plink_log/$alg.UK10K1KG.updatePos_ensembl.log

# allele order mismatch
#awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
#    ($1$2$4 in a && a[$1$2$4] == $6$5)  {print $2,$5,$6,$6,$5}' \
#    $qcdir/$alg.UK10K1KG.updatePos_ensembl.bim \
#    $qcdir/${alg}_GRCh37.p13_not_in_UK10K1KG.txt > \
#    $qcdir/${alg}.UK10K1KG.toUpdateAlleles_ensembl

# possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$alg.UK10K1KG.updatePos_ensembl.bim \
    $qcdir/${alg}_GRCh37.p13_not_in_UK10K1KG.txt > \
    $qcdir/${alg}.UK10K1KG.toFlip_ensembl

#plink --bfile $qcdir/${alg}.UK10K1KG.updatePos_ensembl \
#      --update-alleles $qcdir/${alg}.UK10K1KG.toUpdateAlleles_ensembl \
plink --bfile $qcdir/${alg}.UK10K1KG.updatePos_ensembl \
      --flip $qcdir/${alg}.UK10K1KG.toFlip_ensembl \
      --make-bed \
      --out $qcdir/$alg.UK10K1KG.ensembl_flipped
mv $qcdir/${alg}.UK10K1KG.ensembl_flipped.log \
    $qcdir/plink_log/${alg}.UK10K1KG.ensembl_flipped.log

# allele order mismatch of flipped alleles
#awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
#    ($1$2$4 in a && a[$1$2$4] == $6$5)  {print $2,$5,$6,$6,$5}' \
#    $qcdir/$alg.UK10K1KG.ensembl_flipped.bim \
#    $qcdir/${alg}_GRCh37.p13_not_in_UK10K1KG.txt > \
#    $qcdir/${alg}.UK10K1KG.toUpdate_flipped_ensembl

#plink --bfile $qcdir/${alg}.UK10K1KG.ensembl_flipped \
#      --update-alleles $qcdir/${alg}.UK10K1KG.toUpdate_flipped_ensembl \
#      --make-bed \
#      --out $qcdir/$alg.UK10K1KG.ensembl_updated_flipped
#mv $qcdir/${alg}.UK10K1KG.ensembl_updated_flipped.log \
#    $qcdir/plink_log/${alg}.UK10K1KG.ensembl_updated_flipped.log

# snps that do not match between dataset and reference
#$qcdir/$alg.UK10K1KG.ensembl_updated_flipped.bim \
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$alg.UK10K1KG.ensembl_flipped.bim \
    $qcdir/${alg}_GRCh37.p13_not_in_UK10K1KG.txt > \
    $qcdir/${alg}.toRemove_ensembl

# dataset matched to UK10K1KG and ensembl reference
plink --bfile $qcdir/${alg}.UK10K1KG.ensembl_flipped \
      --exclude $qcdir/${alg}.toRemove_ensembl \
      --make-bed \
      --out $qcdir/$alg
mv $qcdir/${alg}.log $qcdir/plink_log/${alg}.log


