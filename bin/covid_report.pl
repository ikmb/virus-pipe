#!/usr/bin/env perl

use strict;
use Getopt::Long;
use PDF::API2;
use PDF::Table;

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
## PARSE SnpEff VCF FILE
#########################

open (my $IN, '<', $vcf) or die "FATAL: Can't open file: $vcf for reading.\n$!\n";

# MN908947.3      21      .       CAGGTAACAA      GACGGCCAGT      5925.49 

my @records;
my @fields;

my @table = (
	[ "Position","Effekt","Feature_ID","HGVS_c" ]
);

while (<$IN>) {

        chomp;
        my $line = $_;
	if ($line =~ /.<ID=ANN.*/) {
		my $field_string = (split "=",$line)[-1];
		$field_string =~  s/Functional annotations\: //g ;
		$field_string =~ s/[',"]//g ;
		my @temp_fields = split(/\|/, $field_string);
		foreach my $f (@temp_fields) {
			chomp($f);
			$f =~ tr/ //ds;
			push(@fields,$f);
		}
	}

	next if ($line =~ /^#.*/);

	my @elements = split(/\t/, $line);
	my $info_field =  @elements[7] ;
	my @info = split(";", $info_field);	

	my %data;
	foreach my $i (@info) {

		my ($key,$values) = split("=",$i);

		$data{$key} = $values ;
	}

	my @annotations = (split /\|/ ,  $data{"ANN"});

	my %annot;

	foreach my $f (@fields) {
		my $a = shift @annotations;
		$annot{$f} = $a ;
	}

	my $effect = $annot{"Annotation_Impact"} ;
	my $gene_id = $annot{"Feature_ID"} ;
	my $hgvs_p = $annot{"HGVS.c"} ;

	my @te = [ @elements[1] . ":" . @elements[3] . ">" . @elements[4] , 
		$effect ,
		$gene_id , 
		$hgvs_p
	] ;

	push(@table, @te);

}

close($IN);


##################
## BUILD REPORT
##################

my $report_name = $outfile;

my $pdf = PDF::API2->new();
 
$pdf->info(
	'Author' => "IKMB Bioinformatics platform",
	'CreationDate' => localtime,
	'Title' => 'Sars2-CoV2 Report',
);

# Add a page 
my $page = $pdf->page();
$page->mediabox('A4');

 # Prepare a font
my $font = $pdf->corefont('Helvetica');
my $b_font = $pdf->corefont('Helvetica-Bold');

my $gfx = $page->gfx();

my $image = $pdf->image_jpeg("ikmb_bfx_logo.jpg");

$gfx->image($image,50,730,$image->width,$image->height);

my $text = $page->text();
$text->font($font,10);

# Write some text

my $step = 700;

$text->translate(50,$step);
$text->font($b_font,14);
$text->text("Sars-CoV2 Sequenzierung");

$step -= 20;
$text->font($font,10);
$text->translate(50,$step);
$text->text("Analyse zur Identifikation und Typisierung von SARS-Cov2 Viren mittels Genomsequenzierung");

$step -= 30;
$text->font($b_font,10);
$text->translate(50,$step);
$text->text("Patient/Library:");

$text->font($font,10);
$text->translate(250,$step);
$text->text($patient);

$step -= 20;
$text->font($b_font,10);
$text->translate(50,$step);
$text->text("COVID19 Status:");

$text->font($font,10);
$text->translate(250,$step);
$text->text($status);

$step -= 20;
$text->font($b_font,10);
$text->translate(50,$step);
$text->text("Typ/Lineage:");

$text->font($font,10);
$text->translate(250,$step);
$text->text("${global_lineage} (Pangolin)");

$step -= 20;
$text->font($font,8);
$text->translate(50,$step);
$text->text("Pangolin typisiert SARS-Cov2 Genomdaten basierend auf vorab definierten genomischen Varianten");

$step -= 10;
$text->translate(50,$step);
$text->text("https://github.com/cov-lineages/pangolin");

$step -= 40;
$text->font($b_font,12);
$text->translate(50,$step);
$text->text("Metriken");

$text->font($b_font,10);
$step -= 20;
$text->translate(50,$step);
$text->text("Reads sequenziert und gemapped:");

$text->font($font,10);
$text->translate(250,$step);
$text->text($mapped_reads);

$step -= 20;
$text->font($b_font,10);
$text->translate(50,$step);
$text->text("Assembly komplett:");

$text->font($font,10);
$text->translate(250,$step);
$text->text($genome_fraction . "%");

$step -= 80;
$text->font($b_font,12);
$text->translate(50,$step);
$text->text("Beobachtete Varianten");

$text->font($font,8);

# Variant table

$step -= 20;
$text->translate(50,$step);
$text->text("Sars-CoV2 reference: NC_045512.2.");

my $pdftable = new PDF::Table;

my $table_ref = \@table;

my $left_edge_of_table = 50;

$step -= 20;

$pdftable->table(
     $pdf,
     $page,
     $table_ref,
     'x' => $left_edge_of_table,
     'w' => 500,
     'y' => $step,
     'h' => 300,
     'next_y'          => 750,
     'next_h'          => 500,
     'padding'         => 2,
     'padding_right'   => 2,
     'font_size'       => 8,
     'bg_color_odd'    => "gray",
     'bg_color_even'   => "lightgray", 
     'max_word_length' => 50, # 50 between forced splits
);

## Footer ##

$text->font($font,8);

$text->translate(50,70);
$text->text("Dieser Report wurde erstellt mit der IKMB Virus-pipe Pipeline - version 1.0. Weitere Details unter: https://github.com/ikmb/virus-pipe");

$text->translate(50,50);

my $date = localtime;

$text->text("Report erstellt: $date");

$pdf->saveas($report_name);

