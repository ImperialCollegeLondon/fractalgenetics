###############################################################################
###                                                                         ###
### Sample ancestry check based on genotypes                                ###
###                                                                         ###
###     * numbers in QC indicate steps described in the protocol            ###
###      "Anderson et al. (2010) Data quality control in genetic            ###
###       case-control association studies." Nature protocols               ###
###       5(9):1564-73                                                      ###
###                                                                         ###
###                                                                         ###
###     by Hannah Meyer                                                     ###
###                                                                         ###
###############################################################################

####################################
### 1. input parameters and data ###
####################################

### a) parameters
alg=$1
qcdir=$2
highld=$3
reference=$4
hapmap=$5
hapmapdir=$6

mkdir -p $qcdir/plink_log

##########################
### 2. Ancestry check  ###
##########################

# 15. merge study genotypes to HapMap3 data
# Filter hapmap3 data for the same SNPset as in study
plink --bfile  $hapmapdir/$hapmap --extract $qcdir/$alg.prune.in \
      --make-bed --out $qcdir/$hapmap.pruned
mv  $qcdir/$hapmap.pruned.log $qcdir/plink_log/$hapmap.pruned

# check if SNPs which are in the HapMapIII (position matched to UK10K1KG)
# reference sets are concordant with HapMapIII info

# chromosome mismatch (ignore sex chromosomes)
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/$alg.pruned.bim $qcdir/$hapmap.pruned.bim | \
    sed -n '/^[XY]/!p' > $qcdir/${hapmap}.toUpdateChr

plink --bfile $qcdir/${hapmap}.pruned \
      --update-chr $qcdir/${hapmap}.toUpdateChr 1 2 \
      --make-bed \
      --out $qcdir/$hapmap.updateChr
mv $qcdir/$hapmap.updateChr.log $qcdir/plink_log/$hapmap.updateChr.log

# position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $qcdir/$alg.pruned.bim $qcdir/$hapmap.pruned.bim > \
    $qcdir/${hapmap}.toUpdatePos

# possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$alg.pruned.bim $qcdir/$hapmap.pruned.bim > \
    $qcdir/$hapmap.toFlip

# allele order mismatch
#awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
#    ($1$2$4 in a && a[$1$2$4] == $6$5)  {print $2,$5,$6,$6,$5}' \
#    $qcdir/$alg.pruned.bim $qcdir/$hapmap.pruned.bim > \
#    $qcdir/$hapmap.toUpdateAlleles

#--update-alleles $qcdir/${hapmap}.toUpdateAlleles \
plink --bfile $qcdir/$hapmap.updateChr \
      --update-map $qcdir/${hapmap}.toUpdatePos 1 2 \
      --flip $qcdir/${hapmap}.toFlip \
      --make-bed \
      --out $qcdir/${hapmap}.flipped
mv $qcdir/$hapmap.flipped.log \
    $qcdir/plink_log/$hapmap.flipped.log

# snps that do not match between dataset and reference
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$alg.pruned.bim \
    $qcdir/$hapmap.flipped.bim > \
    $qcdir/$hapmap.mismatch

plink --bfile $qcdir/$hapmap.flipped \
      --exclude $qcdir/$hapmap.mismatch \
      --make-bed \
      --out $qcdir/$hapmap.clean
mv $qcdir/$hapmap.clean.log \
    $qcdir/plink_log/$hapmap.clean.log

# merge study genotypes to genotypes of HapMap populations
plink --bfile $qcdir/$alg.pruned  \
      --bmerge $qcdir/$hapmap.clean.bed $qcdir/$hapmap.clean.bim \
         $qcdir/$hapmap.clean.fam  \
      --make-bed --out $qcdir/$alg.merge.$hapmap
mv $qcdir/$alg.merge.$hapmap.log $qcdir/plink_log

# 17. Conduct PCA on the merged data
plink --bfile $qcdir/$alg.merge.$hapmap --pca \
      --out $qcdir/$alg.$reference
mv $qcdir/$alg.$reference.log $qcdir/plink_log



