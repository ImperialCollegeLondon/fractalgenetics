###############################################################################
###                                                                         ###
### Genotype Quality Control                                                ###
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
hapmap=$4
hapmapdi=$5

mkdir -p $qcdir/plink_log

##############
### 2. QC  ###
##############

### a) Identification of Individuals with discordant sex information
# 4.  calculate mean homozygosity rate across X for each individual
plink --bfile $qcdir/$alg --check-sex --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.sexcheck
grep 'PROBLEM' $qcdir/$alg.sexcheck > $qcdir/$alg.failsex

### b) Identification of individuals with elevated missing data rates or
### outlying heterozygosity rate
# 7.  creates files .imiss (4th column: number of missing SNPs N_MISS,
# 6th column: missing SNPs %/individum F_miss) and .lmiss
plink --bfile  $qcdir/$alg --missing --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.missing

# 8.  creates the file .het (3rd column: number of homozygous gt O(Hom),
# 5th column: number of nonmissing genotypes per individual N(NM)
plink --bfile  $qcdir/$alg --het --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.hr

### c) Identification of duplicated or related individuals
# 11.  minimize computational complexity by reducing the number of SNPs in
# IBS matrix through removal of SNPs within LD (r2 threshold 0.05) and MAF
# threshold results in files: .prune.in --> SNPs kept in the analyses;
plink --bfile $qcdir/$alg --exclude range $highld --indep-pairwise 50 5 0.2 \
      --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.prune

# 12.  create pairwise IBS for all pairs of individuals
plink --bfile  $qcdir/$alg --extract $qcdir/$alg.prune.in --make-bed \
      --out $qcdir/$alg.pruned
mv $qcdir/$alg.pruned.log $qcdir/plink_log/$alg.pruned
plink --bfile $qcdir/$alg.pruned --maf 0.1 --genome --out $qcdir/$alg
mv $qcdir/$alg.log $qcdir/plink_log/$alg.pruned.IBS

### d) Identification of individuals of divergent ancestry
# 15. merge study genotypes to HapMap3 data
# Filter hapmap3 data for the same SNPset as in study
plink --bfile  $hapmapdir/$hapmap --extract $qcdir/$alg.prune.in \
      --make-bed --out $qcdir/$hapmap.pruned
mv  $qcdir/$hapmap.pruned.log $qcdir/plink_log/$hapmap.pruned

# check if SNPs which are in the HapMapIII (position matched to UK10K1KG)
# reference sets are concordant with HapMapIII info
# perfect match between SNPs in batch and HapMap
awk 'FNR==NR {{a[$1$2$4$5$6]; next}} $1$2$4$5$6 in a {print $2}' \
    $qcdir/$hapmap.pruned.bim $qcdir/${alg}.pruned.bim > \
    $qcdir/${alg}.match_hapmap

# allele order mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] == $6$5)  {print $2,$5,$6,$6,$5}' \
    $qcdir/$hapmap.pruned.bim $qcdir/${alg}.pruned.bim > \
    $qcdir/${alg}.toUpdate_hapmap

# possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$hapmap.pruned.bim $qcdir/${alg}.pruned.bim > \
    $qcdir/${alg}.toFlip_hapmap

# all snps considered in first match to reference
cat $qcdir/${alg}.match_hapmap $qcdir/${alg}.toFlip_hapmap \
    $qcdir/${alg}.toUpdate_hapmap > $qcdir/${alg}.keep_hapmap

plink --bfile $qcdir/${alg}.pruned \
      --update-alleles $qcdir/${alg}.toUpdate_hapmap \
      --extract $qcdir/${alg}.keep_hapmap \
      --flip $qcdir/${alg}.toFlip_hapmap \
      --make-bed \
      --out $qcdir/$alg.hapmap_flipped
mv $qcdir/${alg}.hapmap_flipped.log \
    $qcdir/plink_log/${alg}.hapmap_flipped.log

# allele order mismatch of flipped alleles
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] == $6$5)  {print $2,$5,$6,$6,$5}' \
    $qcdir/$alg.hapmap_flipped.bim $qcdir/$hapmap.pruned.bim > \
    $qcdir/${alg}.toUpdate_flipped_hapmap

plink --bfile $qcdir/${alg}.hapmap_flipped \
      --update-alleles $qcdir/${alg}.toUpdate_flipped_hapmap \
      --make-bed \
      --out $qcdir/$alg.hapmap_updated_flipped
mv $qcdir/${alg}.hapmap_updated_flipped.log \
    $qcdir/plink_log/${alg}.hapmap_updated_flipped.log

# snps that do not match between dataset and reference
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] != $5$6)  {print $2}' \
    $qcdir/$alg.hapmap_updated_flipped.bim \
    $qcdir/${hapmap}.pruned.bim > \
    $qcdir/${alg}.mismatch_hapmap

plink --bfile $qcdir/${alg}.hapmap_updated_flipped \
      --exclude $qcdir/${alg}.mismatch_hapmap \
      --make-bed \
      --out $qcdir/$alg.hapmap.clean
mv $qcdir/${alg}.hapmap.clean.log \
    $qcdir/plink_log/${alg}.hapmap.clean.log

# merge HapMap-mapped genotypes to genotypes of HapMap populations
plink --bfile $qcdir/$hapmap.pruned  \
      --bmerge $qcdir/$alg.hapmap.clean.bed $qcdir/$alg.hapmap.clean.bim \
         $qcdir/$alg.hapmap.clean.fam  \
      --make-bed --out $qcdir/$alg.merge.$hapmap
mv $qcdir/$alg.merge.$hapmap.log $qcdir/plink_log

# 17. Conduct PCA on the merged data
plink --bfile $qcdir/$alg.merge.$hapmap --pca \
      --out $qcdir/$alg.HapMapIII.pruned.pca
mv $qcdir/$alg.HapMapIII.pruned.pca.log $qcdir/plink_log

