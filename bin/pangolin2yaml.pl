#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd;
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

my %data;
my $dir = getcwd;

my %lookup = ( "AY" => "B.1.617.2",
		"Q" => "B.1.1.7",
		"AZ" => "B.1.1.318",
		"C" => "B.1.1.1.",
		"D" => "B.1.1.25",
		"G" => "B.1.258.2",
		"K" => "B.1.1.277",
		"M" => "B.1.1.294",
		"N" => "B.1.1.33",
		"P" => "B.1.1.28",
		"R" => "B.1.1.316",
		"S" => "B.1.1.217",
		"U" => "B.1.177.60",
		"V" => "B.1.177.54",
		"W" => "B.1.177.53",
		"Y" => "B.1.177.52",
		"Z" => "B.1.177.50",
		"AA" => "B.1.177.15",
		"AB" => "B.1.160.16",
		"AC" => "B.1.1.405",
		"AD" => "B.1.1.315",
		"AE" => "B.1.1.306",
		"AF" => "B.1.1.305",
		"AH" => "B.1.1.241",
);

my $header = qq(
id: 'pangolin_reports'
section_name: 'Pangolin lineage assignment'
plot_type: 'html'
description: '- Assign lineages to libraries using Pangolin'
data: |\n  <dl class="dl-horizontal">
);

printf $header . "\n";

my $entry = "<dt>Library</dt><dd><samp>Lineage</samp></dd>" ;

printf "    $entry\n";

foreach my $file (glob("$dir/*.csv")) {

	my $fh = IO::File->new();
	$fh->open( $file );

	my $f = (split "/" , $file)[-1];
	my $lib = (split /\./ , $f)[0];
	#my $entry = "<dt>Library</dt><dd><samp>$lib</samp></dd>" ;

        #printf "    $entry\n";

	chomp(my @lines = <$fh>);

	my $header = shift @lines ;
	
	foreach my $line (@lines) {

		chomp($line);
		# taxon,lineage,probability,pangoLEARN_version,status,note
		# NODE_1_length_29902_cov_249.978980,B,1.0,2021-01-16,passed_qc,
		# taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note

                my ($seq,$lineage,$conflict,$ambig,$scorpio_call,$scorpio_support,$scorpio_conflict,$vers,$p_vers,$p_learn_vers,$p_vers,$status,$note) = split(",", $line);

		next unless ($status eq "passed_qc");
		chomp($lineage);

		my $trunk = (split /\./, $lineage)[0] ;

		if (exists $lookup{$trunk}) {
			$lineage = $lookup{$trunk};
		}

		my $entry = "<dt>$lib</dt><dd><samp>$lineage</samp></dd>" ;
                printf "    $entry\n";
		
	}

	close($fh);
}

printf "  </dl>\n";

