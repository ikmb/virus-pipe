#!/usr/bin/env perl

use strict;
use Cwd;
use Getopt::Long;
use Excel::Writer::XLSX;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;

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

my $help;

GetOptions(
    "help" => \$help,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

#if ($outfile) {
#    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
#}

die "Must specify an outfile (--outfile)" unless (defined $outfile);

my $dir = getcwd;

# Initiate the XLS workbook
my $workbook = Excel::Writer::XLSX->new($outfile);

# Add a new sheet
my $worksheet = $workbook->add_worksheet();

my $row = 0;
my @bucket;

my @h = ( "Referenz", "Panel","K-Nummer", "Pangolin-Typisierung", "TechnischeValidierung" );
push(@bucket, "Referenz;Panel;K-Nummer;Pangolin-Typisierung;TechnischeValidierung");

&write_xlsx($worksheet, $row, @h);
++$row;

foreach my $file (glob("$dir/*.csv")) {

        my $fh = IO::File->new();
        $fh->open( $file );

	printf STDERR "Reading $file\n";

	my $f = (split "/" , $file)[-1];
        my $lib = (split /\./ , $f)[0];

        chomp(my @lines = <$fh>);

        my $header = shift @lines ;

        foreach my $line (@lines) {

                chomp($line);
                # taxon,lineage,probability,pangoLEARN_version,status,note
		# NEW: taxon,lineage,conflict,pangoLEARN_version,pango_version,status,note
		# NEWER: taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note

                my ($seq,$lineage,$conflict,$ambig,$scorpio_call,$scorpio_support,$scorpio_conflict,$vers,$p_vers,$p_learn_vers,$p_vers,$status,$note) = split(",", $line);
		
                #my ($seq,$lineage,$conflict,$p_vers,$vers,$status,$note) = split(",", $line);

                next unless ($status eq "passed_qc");

		# shorten call into main lineage only
		chomp($lineage);
                my $trunk = (split /\./, $lineage)[0] ;
                if (exists $lookup{$trunk}) {
                        $lineage = $lookup{$trunk};
                }

                my @ele = ( "NC_045512.2", "QIASeq-SARS-CoV-2_Illumina", $seq, $lineage,  "OK" );
                &write_xlsx($worksheet, $row, @ele);
		++$row;
		my $entry = "NC_045512.2;QIASeq-SARS-CoV-2_Illumina;" . $seq . ";" . $lineage . ";" . "OK" ;
		push(@bucket, $entry);
        }

        close($fh);

}

$workbook->close();

foreach my $b (@bucket) {
	printf $b . "\n";
}
sub write_xlsx{
    my ($worksheet, $tem_row, @ele) = @_;
    for(my $i = 0; $i < @ele; ++$i){
        $worksheet->write( $tem_row, $i, $ele[$i]);
    }
}
