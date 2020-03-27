#!/usr/bin/env perl
# convert Kraken reports to MultiQC yaml format

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--report filename]
		The name of the file to read. 
  Ouput:    
    [--gff filename]
        The name of the output file. By default the output is the
        standard output
};

my $report = undef;
my $help;

GetOptions(
    "help" => \$help,
    "report=s" => \$report);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

open (my $IN, '<', $report) or die "FATAL: Can't open file: $report for reading.\n$!\n";

while (<$IN>) {
	chomp;
	my $line = $_;

	my @elements = split("\t", $line);

}
