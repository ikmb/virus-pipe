#!/usr/bin/env perl

use strict;
use Cwd;
use Getopt::Long;
use Excel::Writer::XLSX;
use JSON;

my $alias = undef;

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
    "alias=s" => \$alias,
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

open my $fh, '<', $alias or die "Can't open file $!";

my $file_content = do { local $/; <$fh> };

my $lookup = decode_json($file_content);

my $dir = getcwd;

# Initiate the XLS workbook
my $workbook = Excel::Writer::XLSX->new($outfile);

# Add a new sheet
my $worksheet = $workbook->add_worksheet();

my $row = 0;
my @bucket;

my @h = ( "Referenz", "Panel","K-Nummer", "Pangolin-Typisierung", "Pangolin-Sublinie", "TechnischeValidierung","IMS_ID" );
push(@bucket, "Referenz;Panel;K-Nummer;Pangolin-Typisierung;Pangolin-Sublinie;TechnischeValidierung;IMS_ID");

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

		my ($seq,$lineage,$conflict,$ambig,$scorpio_call,$scorpio_support,$scorpio_conflict,$scorpio_notes,$vers,$p_vers,$s_vers,$c_vers,$designated,$qc_status,$qc_notes,$note) = split(",", $line);
		my $lineage_full = $lineage;

                next unless ($qc_status eq "pass");

		# shorten call into main lineage only
		chomp($lineage);
                my $trunk = (split /\./, $lineage)[0] ;
                if (exists $lookup->{$trunk}) {
			my $match = $lookup->{$trunk};
			if (length $match > 0) {
				# The reference can sometimes be a list of options. We just join it. 
				if (ref($match) eq 'ARRAY') {
					$lineage = (join ",", $match);
				} else {
		                        $lineage = $match;
				}
			}
                }

                my @ele = ( "NC_045512.2", "QIASeq-SARS-CoV-2_Illumina_v2", $seq, $lineage, $lineage_full, "OK" );
                &write_xlsx($worksheet, $row, @ele);
		++$row;
		my $entry = "NC_045512.2;QIASeq-SARS-CoV-2_Illumina_v2;" . $seq . ";" . $lineage . ";" . $lineage_full . ";" . "OK" . ";" ;
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

