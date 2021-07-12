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

my $outfile = undef;
my $infile = undef;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

my $n_count = 0;
my $total_count = 0;

while (my $line = <$IN>) {

	chomp;

	next if ($line =~ /^>.*/);

	my $total_length = length($line);
	my $noNs=$line;

  	$noNs =~ s/[Nn]//g; # remove Ns and ns	
	my $Ns = ($total_length - length($noNs) );

	$total_count += $total_length;
	$n_count += $Ns;
}

printf "Total bases: ${total_count}\n";
printf "N  bases: ${n_count}\n";
