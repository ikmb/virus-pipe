#!/usr/bin/env perl

use strict;
use Getopt::Long;
use PDF::Create;
use PDF::Create::Page;

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
my $pangolin = undef;
my $kraken = undef;
my $bam_stats = undef;
my $assembly_stats = undef;
my $vcf = undef;
my $infile = undef;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "kraken=s" => \$kraken,
    "bam_stats=s" => \$bam_stats,
    "vcf=s" => \$vcf,
    "assembly_stats=s" => \$assembly_stats,
    "pangolin=s" => \$pangolin,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

my $patient = (split /\./, $kraken)[0];

########################
## PARSE PANGOLIN REPORT
########################

open (my $IN, '<', $pangolin) or die "FATAL: Can't open file: $pangolin for reading.\n$!\n";

# parse pangolin report and get lineage assignment
my $global_lineage = "undetermined";

chomp(my @lines = <$IN>);

my $header = shift @lines ;

foreach my $line (@lines) {

	chomp($line);
        # taxon,lineage,probability,pangoLEARN_version,status,note
        # NODE_1_length_29902_cov_249.978980,B,1.0,2021-01-16,passed_qc,
        my ($seq,$lineage,$prob,$vers,$status,$note) = split(",", $line);

        next unless ($status eq "passed_qc");

	$global_lineage = $lineage;	

}

close($IN);

######################
## PARSE KRAKEN REPORT
######################

open (my $IN, '<', $kraken) or die "FATAL: Can't open file: $kraken for reading.\n$!\n";

my $status = "negativ";

while (<$IN>) {

	chomp;
	my $line = $_;

	my ($perc,$num_frag_tax,$num_frag_ass,$rank,$tax_id,$name) = split(/\t+/,$line);

        next unless ($rank =~ /S.*/);

        $name =~ s/^\s+|\s+$//g ;
        $rank =~ s/^\s+|\s+$//g ;

        if ($name =~ "Severe acute respiratory syndrome coronavirus 2") {

 		$status = "positiv";

        }
}

close($IN);

########################
## PARSE SAMTOOLS REPORT
########################

open (my $IN, '<', $bam_stats)  or die "FATAL: Can't open file: $bam_stats for reading.\n$!\n";

# SN      reads mapped and paired:        38084507        # paired-end technology bit set + both mates mapped
my $mapped_reads = "undetermined";

while (<$IN>) {

        chomp;
        my $line = $_;

	if ($line =~ /.*reads mapped and paired.*/) {
		my @elements = split(/\t/, $line);
		$mapped_reads = @elements[2];
	}
}

close($IN);

#########################
## PARSE ASSEMBLY STATS
#########################

open (my $IN, '<', $assembly_stats) or die "FATAL: Can't open file: $assembly_stats for reading.\n$!\n";

my $assembly_length = "undetermined";
my $reference_length = "undetermined";
my $genome_fraction = "undetermined";


while (<$IN>) {

        chomp;
        my $line = $_;

        if ($line =~ /.*Genome fraction.*/) {

		my ($key,$gf) = split(/\t/, $line);
		$genome_fraction = $gf;
	}

}
close($IN);

#########################
## PARSE VCF FILE
#########################

open (my $IN, '<', $vcf) or die "FATAL: Can't open file: $vcf for reading.\n$!\n";

# MN908947.3      21      .       CAGGTAACAA      GACGGCCAGT      5925.49 

my @records;

while (<$IN>) {

        chomp;
        my $line = $_;

	next if ($line =~ /^#.*/);

	my @elements = split(/\t/, $line);

	my $record = @elements[1] . ":" . @elements[3] . ">" . @elements[4] ;

	push(@records,$record);
}

close($IN);


##################
## BUILD REPORT
##################

my $report_name = $outfile;

my $pdf = PDF::Create->new(
    'filename'     => $report_name,
    'Author'       => 'IKMB Diagnostik',
    'Title'        => 'SARS-Cov2 Sequenzierung',
    'CreationDate' => [ localtime ]
);
 
# Add a A4 sized page
my $root = $pdf->new_page('MediaBox' => $pdf->get_page_size('A4'));
 
# Add a page which inherits its attributes from $root
my $page1 = $root->new_page;
 
# Prepare a font
my $font = $pdf->font('BaseFont' => 'Helvetica');

my $logo = $pdf->image("ikmb_bfx_logo.jpg");

$page1->image(
    'image'  => $logo,
    'xscale' => 0.8,
    'yscale' => 0.8,
    'xpos'   => 50,
    'ypos'   => 750
);
 
# Write some text
$page1->line(50, 740,   550, 740);
$page1->stringl($font, 14, 50, 720, 'SARS-Cov2 Sequenzierung');
$page1->stringl($font, 10, 50, 700, 'Analyse zur Identifikation und Typisierung von SARS-Cov2 Viren mittels Genomsequenzierung');
$page1->line(50, 680,   550, 680);

$page1->stringl($font, 12, 50, 660, "Patient/Library:");
$page1->stringl($font, 12, 250, 660, "$patient");

$page1->stringl($font, 12, 50, 640, "COVID19 Status:");
$page1->stringl($font, 12, 250, 640, "$status");


$page1->stringl($font, 12, 50, 620, "Typ/Lineage:"); 
$page1->stringl($font, 12, 250, 620, "${global_lineage} (Pangolin)");

$page1->stringl($font, 10, 50, 600, 'Pangolin typisiert SARS-Cov2 Genomdaten basierend auf vorab definierten genomischen Varianten');
$page1->stringl($font, 10, 50, 585, 'https://github.com/cov-lineages/pangolin');

$page1->stringl($font, 12, 50, 545, 'Metriken');
$page1->stringl($font, 10, 50, 525, 'Reads sequenziert und gemapped:');
$page1->stringl($font, 10, 250, 525, "$mapped_reads");

$page1->stringl($font, 10, 50, 505, 'Genomabdeckung (Assembly):');
$page1->stringl($font, 10, 250, 505, $genome_fraction);

$page1->stringl($font, 12, 50, 480, 'Beobachtete Varianten');

my $pos = 480;
foreach my $r (@records) {
	printf $r . "\n";
	$pos -= 15 ;
	$page1->stringl($font, 10, 50, $pos, $r);

}

$pdf->close;

