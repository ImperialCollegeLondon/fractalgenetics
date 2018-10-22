#!/usr/local/bin/bash

### a) parameters
alg=$1
qcdir=$2
rawdir=$3
center=$4
rawdata=$5
sample=$6

# change working directory to data dir
mkdir -p $qcdir/plink_log

### b)  format raw data and transfer to QC directory ###
if [[ ! -s $qcdir/${alg}.raw.fam ]]; then
    ### Reformat input files if parental information is encoded incorrectly:
    if [[ $center == *"singapore"* ]]; then
        echo "Reformating Singapore files..."
        perl scripts/parse_singapore_data.pl \
            --inped $rawdir/$rawdata \
            --samples $rawdir/${rawdata}_sample_sex.txt \
            --outped $rawdir/$alg.raw \
            --incorrect $sample
    elif [[ $center == *"sanger"* ]]; then
        echo "Reformating Sanger files..."

        # Remove sample prefix "urn:wtsi:[\d\w]_A01_"
        sed -e 's/urn:wtsi:[0-9_ABCDEFGH]*_//g' $rawdir/${rawdata}.fam |
            cut -d " " -f 1,2 | \
            paste <(cut -d " " -f 1,2 $rawdir/${rawdata}.fam) - > \
            $rawdir/${rawdata}.shortIDs
        plink --bfile $rawdir/$rawdata \
              --update-ids $rawdir/${rawdata}.shortIDs \
              --make-bed \
              --out $rawdir/$rawdata.shortIDs
        mv $rawdir/$rawdata.shortIDs.log $rawdir/plink_log

        # Replace sanger IDs (internal genotype IDs) with BRU Ids
        awk 'FNR==NR {a[$1]=$2; next} $1 in a \
            {id1=$1;id2=$2;$1=a[id1];$2=a[id1]; print id1,id2,$1,$2}' \
            $sample $rawdir/${rawdata}.shortIDs.fam > \
            $rawdir/${rawdata}.updateIDs
        plink --bfile $rawdir/$rawdata.shortIDs \
              --update-ids $rawdir/${rawdata}.updateIDs \
              --make-bed \
              --out $rawdir/$rawdata.updateIDs
        mv $rawdir/$rawdata.shortIDs.log $rawdir/plink_log

        # Update parental information
        awk '{print $1,$2,0,0}' $rawdir/${rawdata}.updateIDs.fam > \
            $rawdir/${rawdata}.updateParents
        plink --bfile $rawdir/$rawdata.updateIDs \
              --update-parents $rawdir/${rawdata}.updateParents \
              --make-bed \
              --out $rawdir/$alg.raw
        mv $rawdir/$alg.raw.log $rawdir/plink_log
    else
        echo "Center $center no known, break"
        break
    fi
    ### Copy reformated data from rawdir to qcdir
    cp $rawdir/$alg.raw.bim $qcdir/$alg.raw.bim
    cp $rawdir/$alg.raw.fam $qcdir/$alg.raw.fam
    cp $rawdir/$alg.raw.bed $qcdir/$alg.raw.bed
fi

### c) Check for duplicate IDs in dataset ###
cut -d " " -f 1 $qcdir/$alg.raw.fam | sort | uniq -d > $qcdir/$alg.duplicate.IDs
if [[ -s $qcdir/$alg.duplicate.IDs ]]; then
    echo "Remove duplicate ID's..."
    plink --bfile $qcdir/$alg.raw --recode --out $qcdir/$alg.raw
    mv $qcdir/${alg}.raw.log $qcdir/plink_log/${alg}.raw.recode.log
    awk  'FNR==NR {a[$1]; next} !($1 in a)' $qcdir/$alg.duplicate.IDs \
        $qcdir/$alg.raw.ped  > $qcdir/$alg.raw.duplicate.remove.ped
    plink --ped $qcdir/$alg.raw.duplicate.remove.ped \
        --map $qcdir/$alg.raw.map --make-bed --out $qcdir/${alg}.raw
    mv $qcdir/${alg}.raw.log $qcdir/plink_log/${alg}.raw.duplicateremove.log
fi

### d)  Remove 'questionable' IDs if existent (files generated via IDmatching.R)
if [[ -s $qcdir/$alg.preQCfail.IDs ]]; then
    echo "Remove preQCfail ID's..."
    plink --bfile $qcdir/$alg.raw \
          --remove $qcdir/$alg.preQCfail.IDs \
          --make-bed --out $qcdir/$alg.raw
    mv $qcdir/${alg}.raw.log $qcdir/plink_log/${alg}.raw.questionableIDs.log
fi

