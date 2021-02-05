#!/usr/bin/env perl

use strict;
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

die "Must specify an outfile (--outfile)" unless (defined $outfile);

# Initiate the XLS workbook
my $workbook = Excel::Writer::XLSX->new($outfile);

# Add a new sheet
my $worksheet = $workbook->add_worksheet();

my $row = 0;

my $fh = IO::File->new();
$fh->open( $infile );

foreach my $file (glob("$dir/*.csv")) {

        my $fh = IO::File->new();
        $fh->open( $file );

	my $f = (split "/" , $file)[-1];
        my $lib = (split /\./ , $f)[0];

        chomp(my @lines = <$fh>);

        my $header = shift @lines ;

        foreach my $line (@lines) {

                chomp($line);
                # taxon,lineage,probability,pangoLEARN_version,status,note
                # NODE_1_length_29902_cov_249.978980,B,1.0,2021-01-16,passed_qc,
                my ($seq,$lineage,$prob,$vers,$status,$note) = split(",", $line);

                next unless ($status eq "passed_qc");

                my @ele = ( $lib, $lineage, $prob );
                &write_xlsx($worksheet, $row, @ele);
		++$row;
        }

        close($fh);

}

sub write_xlsx{
    my ($worksheet, $tem_row, @ele) = @_;
    for(my $i = 0; $i < @ele; ++$i){
        $worksheet->write( $tem_row, $i, $ele[$i]);
    }
}
