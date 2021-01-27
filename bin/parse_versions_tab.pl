#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd qw(getcwd);

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $help;

GetOptions(
    "help" => \$help,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

my $directory = getcwd;

opendir (DIR, $directory) or die $!;

my $version = undef;
my $tool = undef;

while (my $file = readdir(DIR)) {

	next unless ($file =~ /^v_.*\.txt$/);
	open (my $IN, '<', $file) or die "FATAL: Can't open file: $file for reading.\n$!\n";
	
	chomp(my @lines = <$IN>);	
	
	if ($file =~ /v_gatk.*/) {
		my $line = @lines[0];
		$tool = "GATK4" ;
		$version = (split " ", $line)[-1];
		
	} elsif ($file =~ /^v_nextflow\.txt$/ )  {
		my $line = @lines[0];
		$tool = "Nextflow";
		$version = $line;
		next;
	} elsif ($file =~ /^v_ikmb_virus_pipe\.txt$/) {
		my $line = @lines[0];
                $tool = "Sars-CoV2 Pipeline";
		$version = $line;
	} elsif ($file =~ /v_freebayes.*/) {
                my $line = @lines[0];
                my @elements = split(" ",$line);
                $version = @elements[-1];
                $tool = "Freebayes";
	} elsif ($file =~ /^v_picard\.txt/) {
		my $line = @lines[0];
                $tool = "Picard";
                $version = (split " ", $line)[-1];
	} elsif ($file =~ /^v_bwa\.txt/) {
		my $line = @lines[2];
                $tool = "BWA";
                $version = (split " ", $line)[-1];
	} elsif ($file =~ /^v_bowtie2.txt/) {
		my $line = @lines[0];
		my @elements = (split " ",$line);
		$version = @elements[-1];
		$tool = "Bowtie2";
	} else {
		my $line = @lines[0];
		my @elements = (split " ",$line);
		$tool = @elements[0];
		$tool =~ s/\,//;
		$version = @elements[-1];
	}

	printf $tool . "\t" . $version . "\n";	
	
	close($IN);
	
}

close(DIR);


