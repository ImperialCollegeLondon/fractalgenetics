#!/nfs/research2/software/prefix/usr/bin/perl

#######################################################################
### 	script based on gwas_parse.pl,								###
###  	kindly provided by Roddy Walsh to format raw data for plink ###
#######################################################################

################
### modules  ###
################

use strict;
use warnings;
use diagnostics;

use Getopt::Long;


#################
### variables ###
#################
my $in_correct="";
my ($out_ped, $samples, $in_ped, $plink_core);
my (%sample_sex, %sample_swap);

############################
### example command line ###
############################

# perl ~/GWAS/analysis/genotyping/parse_singapore_data.pl --inped ~/GWAS/data/genotype/MRI_genotype/omnix_NHCS_20160217/Heart_StudyStuartCook_Omniexpress_08052015  --samples ~/GWAS/data/genotype/MRI_genotype/omnix_NHCS_20160217/Heart_StudyStuartCook_Omniexpress_08052015_sample_sex.txt --outped ~/GWAS/data/genotype/MRI_genotype/QC/gencall.singapore12/gencall.singapore12.raw

# perl ~/GWAS/analysis/genotyping/parse_singapore_data.pl --inped ~/GWAS/data/genotype/MRI_genotype/omnix_NHCS_20160217/Edmund_NHCS_Heart_Failure  --samples ~/GWAS/data/genotype/MRI_genotype/omnix_NHCS_20160217/Edmund_NHCS_Heart_Failure_sample_sex.txt --outped ~/GWAS/data/genotype/MRI_genotype/QC/gencall.singapore3/gencall.singapore3.raw --incorrect ~/GWAS/data/genotype/MRI_genotype/omnix_NHCS_20160217/corrected_sample_id.txt

###############
### program ###
###############

# get command line parameters
GetOptions ("outped=s" => \$out_ped, "samples=s" => \$samples, "inped=s" => \$in_ped, "incorrect=s" => \$in_correct);

# Output file location
open (RES, ">$out_ped".".ped") or die "can't write to output file $out_ped.ped, reason: $!";
open (MISS, ">$out_ped".".missing") or die "can't write to output file $out_ped.missing, reason: $!";

# Sample ID and sex file (format: "10AB01234 1")
open (FILE, "$samples") or die "can't open sample file $samples, reason $!";
while (<FILE>)	{
  if (/^(\S+)\s([12FM]{1})/) {
    my $id = $1;
    my $sex = $2;
    if ($sex =~/(F|M)/) {
        $sex=2 if $sex eq "F";
        $sex=1 if $sex eq "M";
    }
    $sample_sex{$id}=$sex;
    #print $sex, "\t", $id, "\n";
  }
}
close FILE;

if ($in_correct ne "") {
    open (CORR, "$in_correct") or die "can't open ped file $in_correct, reason $!";
    while (<CORR>) {
        if (/^\d{1,3}\t969_(\d{2}[A-Z]{2}\d+)_\S+\t\d{1,3}\t(\d{2}[A-Z]{2}\d+)/) {
            $sample_swap{$1}=$2;
            #print "$1\t$2\n";
        }
    }
    close CORR;
}
# Parse through PED file and add sex data for each sample, as well as cleaning up the sample ID
open (FILE, "$in_ped".".ped") or die "can't open ped file $in_ped.ped, reason $!";
while (<FILE>) {
  if (/^(\S+)\t969_(\d{2}[A-Z]{2}\d+)_\S+\t0\t0\t0\t0\t(.+)/) {
    my $sample_id=$2;
    my $snps=$3;
    if(exists($sample_swap{$sample_id})) {
        #print $sample_id, "\t", $sample_swap{$sample_id}, "\n";
        $sample_id = $sample_swap{$sample_id};
    }
    if (exists($sample_sex{$sample_id})) {
        print RES "$sample_id\t$sample_id\t0\t0\t$sample_sex{$sample_id}\t0\t$snps\n";
    } else {
        print MISS $sample_id, "\tsample look-up $samples \n";
    }
  } else {
    print $_, " Sample ID doesn't match format\n";
  }
}
close FILE;
close MISS;

# rename .map file according to output .ped file
system("cp $in_ped.map $out_ped.map");
system("plink --noweb --file $out_ped --make-bed --out $out_ped");

__END__
  
  
  
