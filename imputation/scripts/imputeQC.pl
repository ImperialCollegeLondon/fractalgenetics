#! /usr/bin/perl

#################
#### packages ###
#################

use strict;
use warnings;
use diagnostics;

use Getopt::Long;

#################
### variables ###
#################

my ($indir, $refdir, $outdir, $name, $nosnp, $verbose, $concordance_flag,
    $analysis_flag, $SNP_ref, $SNP_ref_sample, $SNP_sample, $SNP_total,
    $concordance, $concordance_overall, $chrStr);
my (@data);

############################
### example command line ###
############################

# imputeQC.pl --indir ~/data/genotype/imputation/combined --outdir ~/data/genotype/imputation/combined/counts --refdir ~/data/hmeyer/UK10K1000Genomes --nosnp 1.48,9.22,21.1 --verbose --name gencall.combined.clean.related
###############
### program ###
###############

GetOptions ("indir=s" => \$indir, "refdir=s" => \$refdir, "outdir=s" => \$outdir,
            "nosnp=s" => \$nosnp, "name=s" => \$name, "verbose" => \$verbose);


open (OUT_COUNTSALL, ">$outdir/SNPsPerChr.txt") or die "Can't write to filehandle $outdir/SNPsPerChr.txt, reason: $!";
print OUT_COUNTSALL "Chr\tgenotyped\timputed\timputed QC\n";
foreach  my $chr (1..22, "X") {
    #foreach  my $chr (1..22, "X_PAR1", "X_PAR2", "X") {
    $chrStr = (length($chr) == 1 && $chr ne "X") ? "0$chr" : $chr;
    open (OUT_NUMBERS, ">$outdir/chr$chrStr".".chunkStats.txt") or die "Can't write to filehandle $outdir/chr$chrStr"."_numbers.txt, reason: $!";
    open (OUT_DATA, ">$outdir/chr$chrStr".".chunkConcordance.txt") or die "Can't write to filehandle $outdir/chr$chrStr"."_data.txt, reason: $!";
    print  OUT_NUMBERS "Chr\tChunk\tSNP_total\tSNP_ref\tSNP_sample\tSNP_ref_sample\tconcordance_overall\n";
    print  OUT_DATA "Chr\tChunk\tStart\tEnd\tGenotypes\tConcordance\n";

    # Chunks per chromosome
    open (REFFILE, "$refdir/chunkBoundariesChr$chr.txt") or die "Can't open filehandle $refdir/chunkBoundariesChr$chr.txt, reason $!\n";
    my @chunk_num = ();
    while (defined(my $line =<REFFILE>)) {
        chomp $line;
        my @array = split (/\t/, $line);
        if (scalar @array > 2) {
            if (defined $nosnp) {
                my %missing = split /[,\.]/, $nosnp;
                unless ($missing{$array[0]} && $missing{$array[0]} == $array[1]) {
                        push @chunk_num, $array[1];
                }
            } else {
                push @chunk_num, $array[1];
            }
        }
    }
    # 1. Imputation quality per chunk
    print "Imputation quality for chr"."$chr\n" if $verbose;
    foreach my $chunk (@chunk_num) {
        my $file = "$indir/imputed/chr$chr/$name.chr$chr.$chunk.gen_summary";
        open (INFILE, "$file") or die "Can't open filehandle $file, reason $!\n";
        ($analysis_flag, $SNP_ref, $SNP_sample, $SNP_ref_sample, $SNP_total, $concordance_flag, $concordance, $concordance_overall) = (0,0,0,0,0,0,0,0);
        @data =();
        while (defined(my $line =<INFILE>)) {
            chomp $line;
            if ($analysis_flag == 1) {
                $SNP_ref = $line =~ s/\s*--(\d*)\s*type\s*0\s*SNPs/$1/r if $line =~ m/\s*--(\d*)\s*type\s*0\s*SNPs/;
                $SNP_ref_sample = $line  =~ s/\s*--(\d*)\s*type\s*2\s*SNPs/$1/r if $line =~ m/\s*--(\d*)\s*type\s*2\s*SNPs/;
                $SNP_sample = $line =~ s/\s*--(\d*)\s*type\s*3\s*SNPs/$1/r if $line =~ m/\s*--(\d*)\s*type\s*3\s*SNPs/;
                $SNP_total = $line =~ s/\s*--(\d*)\s*total\s*SNPs/$1/r if $line =~ m/\s*--(\d*)\s*total\s*SNPs/;
                if ($SNP_total ne "0") {
                   $analysis_flag = 0;
                }
            }
            if ($line =~m/\s*\[0.0-0.1\]/) {
                $concordance_flag = 1;
            }
            if ($concordance_flag == 1) {
                if ($line =~ m/.*0\.0\]\s*[\d\.]*\s*([\d\.]*).*/) {
                     $concordance_overall = $line =~ s/.*0\.0\]\s*[\d\.]*\s*([\d\.]*).*/$1/r;
                }
                $concordance = $line =~ s/\s*\[([\d\.]*)-([\d\.]*)\]\s*(\d*)\s*([\d\.]*).*/$1\t$2\t$3\t$4\n/r;
                push @data, "$chr\t$chunk\t$concordance";
            }
            if ($line =~m/\s*\[0.9-1.0\]/) {
                $concordance_flag = 0;
            }
            if ($line =~ m/.*Analysis region.*/) {
                $analysis_flag = 1;
            }
        }
        print  OUT_NUMBERS "$chr\t$chunk\t$SNP_total\t$SNP_ref\t$SNP_sample\t$SNP_ref_sample\t$concordance_overall\n";
        print OUT_DATA @data;
        close (INFILE);
    }
    # 2. SNP counts per chromosome: genotyped, imputation, QC
    print "Marker counts for chr"."$chr\n" if $verbose;
    open (OUT_COUNTSCHR, ">$outdir/chr$chr".".SNPsPerChunk.txt") or die "Can't write to filehandle $outdir/chr$chr".".SNPsPerChunk.txt, reason: $!";
    print OUT_COUNTSCHR "Chunk\timputed\n";
    my $perChrImpute = 0;
    foreach my $chunk  (@chunk_num) {
        my $perChunkImpute = `wc -l $indir/imputed/chr$chr/*chr$chr.$chunk.gen | cut -d " "  -f 1`;
        chomp $perChunkImpute;
        $perChrImpute = $perChrImpute + $perChunkImpute;
        print OUT_COUNTSCHR "$chunk\t$perChunkImpute\n";
    }
    my $perChrGenotyped = `wc -l $indir/unphased/*chr$chr.bim | cut -d " "  -f 1`;
    my $perChrImputeQC = `zcat $indir/genotypes/*chr$chr.qc.gen.gz | wc -l | cut -d " " -f 1`;
    chomp $perChrGenotyped;
    chomp $perChrImputeQC;
    print OUT_COUNTSALL "$chr\t$perChrGenotyped\t$perChrImpute\t$perChrImputeQC\n";

    close (REFFILE);
    close (OUT_DATA);
    close (OUT_NUMBERS);
    close (OUT_COUNTSCHR);
}
close (OUT_COUNTSALL);
