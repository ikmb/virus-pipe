#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $assembly = undef;
my $assembly_stats = undef;
my $coverage = undef;

my $max_missing = 5;
my $min_target_cov = 95;

my $help;

GetOptions(
    "help" => \$help,
    "assembly=s" => \$assembly,
    "assembly_stats=s" => \$assembly_stats,
    "coverage=s" => \$coverage);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

########################
## PARSE MOSDEPTH REPORT
#######################

open (my $IN, '<', $coverage) or die "FATAL: Can't open file: $coverage for reading.\n$!\n";

my $target_cov = "undetermined";

chomp(my @lines = <$IN>);

foreach my $line (@lines) {

        chomp($line);
	
	# Get the fraction of the genome covered at 20X
	if ($line =~ /^total	20	.*/) {
		my ($t,$c,$target_coverage) = split(/\t/,$line);
		$target_cov = $target_coverage*100;
	}

}
close($IN);

#########################
## PARSE ASSEMBLY STATS
#########################

open (my $IN, '<', $assembly_stats) or die "FATAL: Can't open file: $assembly_stats for reading.\n$!\n";

my $assembly_length = "undetermined";
my $assembly_gaps = "undetermined";
my $genome_fraction = "undetermined";

while (<$IN>) {

        chomp;
        my $line = $_;

        if ($line =~ /.*Nb of nucleotides \(counting.*/) {
		($assembly_length) = $line =~ /(\d+)/;
	} elsif ( $line =~ /.Nb of Ns./) {
		($assembly_gaps) = $line =~ /(\d+)/;
		$genome_fraction = ($assembly_gaps/$assembly_length) ;
	}
}
close($IN);

my $rounded = 100*(sprintf "%.2f", $genome_fraction);

if ($rounded <= $max_missing && $target_cov >= $min_target_cov) {

	system("cp $assembly 00PASS/");

} else {
	system("cp $assembly 00FAIL/");
}
	
