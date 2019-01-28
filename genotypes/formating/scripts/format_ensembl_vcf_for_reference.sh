###################################################################
###################################################################
###                                                             ###
###     Script to generate ENSEMBL-annotated reference genome   ###
###     for genotyping QC                                       ###
###                                                             ###
###     * downloads ensembl ref version of choice               ###
###     * compares to UK10K1KG reference panel                  ###
###     * generates ref panel of ENSEMBL but not                ###
###       UK10K1KG annotated SNPS                               ###
###                                                             ###
###     by Hannah Meyer                                         ###
###                                                             ###
###################################################################
###################################################################

### 1. source config file with input parameters
. $1

echo -e "* datadir (directory to download annotation and save formated file to): $datadir\n
         * ensemblfile (filename under which ensembl file will be saved): $ensemblfile\n
         * ensembl_ftp_url (url for ensembl download): $ensembl_ftp_url\n
         * UK10K1KGdir (path to UK10K1KG reference dir containing genome.legend file): $UK10K1KGdir"

### 2. Ensembl human variation vcf file
mkdir -p $datadir/log

# a) Download
wget $ensembl_ftp_url -O ${ensemblfile}_human.vcf.gz -o \
    $datadir/log/get_$ensemblfile.log
gunzip ${ensemblfile}_human.vcf.gz

# b) Format VCF file to be in .bim format
awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2, $3, $4, $5}' \
    $datadir/${ensemblfile}_human.vcf  > \
    $datadir/${ensemblfile}_human_SNP_tmp.txt

# c) Remove IDs that only have one designated allele (only 12 variants for GrCh37p13)
awk '$5 != ""'  $datadir/${ensemblfile}_SNP_tmp.txt > \
    $datadir/${ensemblfile}_SNP.txt

### 3. Compare ensembl annotation to UK10K1KG reference annotation
# Takes forever, only rerun if totally necessary
plink1.07 --noweb --id-match $datadir/${ensemblfile}_SNP.txt SNP,3 \
    $UK10K1KGdir/genome.legend.format SNP,2 + complete \
    --out $datadir/UK10K_EnsemblGRCh37p13

# a) check IDs where chromosomes don't match
awk '$1 != $6 {print $3}' $datadir/UK10K_EnsemblGRCh37p13.matched > \
    $datadir/UK10K_EnsemblGRCh37p13.chromosome.mismatch

# b) check IDs where positions don't match
awk '$2 != $9 {print $3}' $datadir/UK10K_EnsemblGRCh37p13.matched > \
    $datadir/UK10K_EnsemblGRCh37p13.position.mismatch

# c) check IDs where alleles don't match
awk '($4$5 != $10$11) && ($4$5 != $10$11) {print $3}' \
    $datadir/UK10K_EnsemblGRCh37p13.matched > \
    $datadir/UK10K_EnsemblGRCh37p13.allele.mismatch

# d) combine mismatch IDs
cat $datadir/UK10K_EnsemblGRCh37p13.chromosome.mismatch \
    $datadir/UK10K_EnsemblGRCh37p13.position.mismatch \
    $datadir/UK10K_EnsemblGRCh37p13.chromosome.mismatch | sort | \
    uniq > $datadir/mismatch_rsIDs

# e) remove mismatch IDs from ensembl data
awk 'FNR==NR {a[$1]++; next} !a[$3]' \
    $datadir/mismatch_rsIDs $datadir/${ensemblfile}_SNP.txt > \
    $datadir/${ensemblfile}_UK10K1000Genomes_filtered.txt

### 4. generate file with SNPs annotated in ensembl but not in UK10K1KG ref set
awk 'FNR==NR {a[$3]++; next} !a[$3]' $datadir/UK10K_EnsemblGRCh37p13.matched \
    $datadir/${ensemblfile}_UK10K1000Genomes_filtered.txt > \
    $datadir/${ensemblfile}_not_in_UK10K1000Genomes_all.txt

### 5. Remove any rsIDs that have double sex chromosome annotations
cut -f 1-3  $datadir/${ensemblfile}_not_in_UK10K1000Genomes_all.txt | \
    sort -k3 | uniq -f 2 -D > $datadir/${ensemblfile}_doubleSex.txt

awk 'FNR==NR {a[$3]++; next} !($3 in a)' $datadir/${ensemblfile}_doubleSex.txt \
    $datadir/${ensemblfile}_not_in_UK10K1000Genomes_all.txt > \
    $datadir/${ensemblfile}_not_in_UK10K1000Genomes.txt




