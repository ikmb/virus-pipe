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

		my $entry = "<dt>$lib</dt><dd><samp>$lineage</samp></dd>" ;
                printf "    $entry\n";
		
	}

	close($fh);
}

printf "  </dl>\n";

