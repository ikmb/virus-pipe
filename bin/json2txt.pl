#!/usr/bin/env perl

use strict;
use Getopt::Long;
use JSON;
use Data::Dumper;

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

my $json_text = do {
   open(my $json_fh, "<:encoding(UTF-8)", $infile)
      or die("Can't open \"$infile\": $!\n");
   local $/;
   <$json_fh>
};

my $json = JSON->new;
my $data = $json->decode($json_text);

my $variants = $data->{'Variants'};

printf "Gene;Position;Effect;RefBase;AltBase;Annotation;HGVS_P\n";

foreach my $v (@$variants) {
	printf $v->{'gene_name'} . ";" . $v->{'genomic_position'} . ";". $v->{'effect'} . ";" . $v->{'ref_base'} . ";" . $v->{'alt_base'} . ";" . $v->{'annotation'} . ";" . $v->{'hgvs_p'}  . "\n";
}
