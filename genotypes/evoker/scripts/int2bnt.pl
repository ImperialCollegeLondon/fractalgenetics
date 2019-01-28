#!/software/bin/perl

## Description: This script generates binary intensity files in the format required by Evoker
## Note: We assume that the input intensity file has the same snps and samples as the bim and fam files for the mathing bed file
##
## Usage: >./int2bnt.pl --intensities collection.chr.int --filetype="illuminus" --out output_filename {--samples samples.fam}
## Input: Intensity file in one of the accepted formats named in the form collection.chromosome.int
## Output: Binary Intensity file collection.chromosome.bnt
## Arguments: 
## -i --input input intensity file 
## -f --filetype [chaimo | affy | illuminus]
## -o --output output file name and path
## -s --samples .fam file for sample to include in output
## default format:
##	A matrix of intensities with SNPs as rows and individuals as pairs of whitespaceÐseparated columns. 
##Êchaimo input format:  
##Ê	Tab-delimited plain text, one line per SNP, consisting of AffyID, RSID, position, alleleA, alleleB and one pair of intensities per sample for each of the two alleles
## affy birdsuite format: 
##	Birdsuite allele_summary file, which has the intensities of each allele of each SNP in matrix format. (two lines per SNP, one for each allele)
## illuminus and beagle format: 
##	Tab-delimited plain text, one line per SNP, consisting of ID, position, alleles and one pair of intensities per sample for each of the two alleles
## 
##ÊTODO: deal with NA values in chiamo and affy files

use strict;
use Getopt::Long;

my $inputfile = '';
my $samples   = '';
my $filetype  = '';
my $outfile   = '';

GetOptions(	'input=s'    => \$inputfile,
			'samples=s'  => \$samples,
			'filetype=s' => \$filetype,
			'output=s'   => \$outfile 
		  );
		  
unless ($inputfile && $outfile && $filetype) {
	die "Missing required arguments: -i <intensity file> -o <output file> -f <filetype>";	
}

my $out_fh;
if ($outfile =~ /\.bnt$/i) {
	open ($out_fh, ">$outfile") or die "Can't open output '$outfile': $!";
} else {
	open ($out_fh, ">$outfile.bnt") or die "Can't open output '$outfile.bnt': $!";	
}
	
## magic number to ensure the binary is a real evoker file not just garbage
print $out_fh pack('B*',"0001101000110001");
	
if ($inputfile =~ /\.gz$/){
	open (IN, "zcat $inputfile |") or die "Can't open '>zcat $inputfile': $!";
}else{
	open (IN, $inputfile) or die "Error: Can't open '$inputfile': $!";
}

if ($samples) {
	if ($filetype =~ /illuminus/i || $filetype =~ /beagle/i) {
  	
  		my $aSamples = parse_sample_file($samples);
		my $header   = <IN>;  		
  		  		
  		my $aSampPos = sample_position_array($header, $aSamples);
  	
		while (my $line = <IN>){
  			chomp($line);
  			my @fields = split(/\s+/, $line);
  		
  			for my $pos (@$aSampPos) {
				my $int = $fields[$pos];
   	 			if ($int eq 'NaN') {
   	 				print $out_fh pack('f*', -1);
   	 			} else {
   	 				print $out_fh pack('f*', $int);	
   	 			}
			}
		}
	} else {
		die "Sample filtering is not supported for file type '$filetype'\n";
	}	
}

if ($filetype =~ /chiamo/i) {	 
	my $header = <IN>;
	while (my $line = <IN>){
  		chomp($line);
  		my @fields = split(/\s+/, $line);
  		for (my $i = 5; $i < scalar(@fields); $i++){
   	 		print $out_fh pack('f*', $fields[$i]);
  		}  		
	}
} elsif ($filetype =~ /affy/i) {
	my %header_info;
	my @allele_a;
	my @allele_b;
	my $header_line = <IN>;
	while($header_line =~ /^#%(.+)\n$/) {
		my ($param, $val) = split(/=/, $1);
		$header_info{$param} = $val; 
		$header_line = <IN>;	
	}
	## the line which broke out of the while loop contains the column headings
	my $col_headings = $header_line;
	## generate a family file for testing
	
	
	my $SAMPLE_NUM   = $header_info{'affymetrix-algorithm-param-apt-opt-cel-count'};

	while (my $line = <IN>) {
		chomp($line);
		my @allele = split(/\t/, $line);
		my $snp_id = $allele[0];
		if ($snp_id =~ /-A$/i) {
			@allele_a = @allele;
		} elsif ($snp_id =~ /-B$/i) {
			@allele_b = @allele;
		}
		## when the A and B allele data for a SNP is parsed print the intensity values out, this is tested by both allele arrays being the correct size  (+1 is the SNP id)
		if (scalar(@allele_a) == $SAMPLE_NUM + 1 && scalar(@allele_b) == $SAMPLE_NUM + 1) {
			my $snp_id_a = shift(@allele_a);
			my $snp_id_b = shift(@allele_b);
			## remove the trailing allele
			$snp_id_a =~ s/-A$//;
			$snp_id_b =~ s/-B$//;
			if ($snp_id_a eq $snp_id_b) {
				for (my $i=0; $i<@allele_a; $i++) {
					print $out_fh pack('f*', ($allele_a[$i],$allele_b[$i]));
				}
				@allele_a = ();
				@allele_b = ();	
			} else {
				die "Affy Error: IDs do not match '$snp_id_a' '$snp_id_b'";
			}
		}
	}
} elsif ($filetype =~ /illuminus/i || $filetype =~ /beagle/i) {
  	my $header = <IN>;
	while (my $line = <IN>){
  		chomp($line);
  		my @fields = split(/\s+/, $line);
  		for (my $i = 3; $i < scalar(@fields); $i++){
   	 		my $int = $fields[$i];
   	 		if ($int eq 'NaN') {
   	 			print $out_fh pack('f*', -1);
   	 		} else {
   	 			print $out_fh pack('f*', $int);	
   	 		}
  		}  		
	}						
} elsif ($filetype =~ /default/) {
	<IN>;
	while (<IN>){
  		my @fields = split; 		
  		for (my $i = 1; $i <= $#fields; $i++){
   	 		print $out_fh pack('f*',$fields[$i]);
  		}	
	}
}

close IN;
close OUT;
  	
sub parse_sample_file {
	my $file = shift;
	open(FAM, $file) or die "Can't open '$file': $!\n";
	my @samples;
	while (my $line = <FAM>) {
		chomp $line;
		my @vals = split(/\s+/, $line);
		push(@samples, $vals[1]);
	}
	close FAM;
	return \@samples;
}  	

##Êcreate a sample_position array to match each of the samples to a position on the vcf line
sub sample_position_array {
	my $header   = shift;
	my $aSamples = shift;
	
	my @sample_position;
	my @hsamples = split( /\s+/, $header );
	## create a hash of all the samples with their position on the line
	my %sample_pos;
	for (my $i = 0; $i <= $#hsamples; $i++) {
		$sample_pos{$hsamples[$i]} = $i;
	}
			
	for my $sample (@$aSamples) {
		## get both A and B
		for my $channel ('A','B') {
			my $test = $sample.$channel;
			if (exists $sample_pos{$test}) {
				push(@sample_position, $sample_pos{$test});
				print "$test $sample_pos{$test}\n";	
			} else {
				warn "Sample '$test' not in data\n";
			}	
		}
	}
	return \@sample_position;
}
