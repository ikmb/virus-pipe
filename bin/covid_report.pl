#!/usr/bin/env perl

use strict;
use Getopt::Long;
use PDF::API2;
use PDF::Table;
use JSON;
use  POSIX;

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
my $depth = undef;
my $bam_stats = undef;
my $software = undef;
my $assembly_stats = undef;
my $patient = undef;
my $vcf = undef;
my $plot = undef;
my $infile = undef;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "kraken=s" => \$kraken,
    "patient=s" => \$patient,
    "depth=s" => \$depth,
    "software=s" => \$software,
    "bam_stats=s" => \$bam_stats,
    "vcf=s" => \$vcf,
    "assembly_stats=s" => \$assembly_stats,
    "pangolin=s" => \$pangolin,
    "plot=s" => \$plot,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

my $library = (split /\./, $kraken)[0];

my $date = localtime;

my %data; # hold information for JSON reporting

$data{"Sample"}= {"library" => $library, "patient" => $patient} ;
$data{"Reference"} = "NC_045512.2.";
$data{"timestamp"} = $date;

##########################
## PARSE SOFTWARE VERSIONS
##########################

open (my $IN, '<', $software) or die "FATAL: Can't open file: $software for reading.\n$!\n";

my %tools;
my $pipeline_version;

while (<$IN>) {

        chomp;
        my $line = $_;
	my ($tool,$version) = split(/\t/,$line);
	if ($tool =~ /.*Pipeline.*/) {
		$pipeline_version = $version;
	} else {
		$tools{$tool} = $version;
	}

}

close($IN);

$data{"Software"} = { %tools };
$data{"Version"} = $pipeline_version;

########################
## PARSE MOSDEPTH REPORT
#######################

open (my $IN, '<', $depth) or die "FATAL: Can't open file: $depth for reading.\n$!\n";

my $target_cov = "0";

chomp(my @lines = <$IN>);

foreach my $line (@lines) {

        chomp($line);
	
	# Get the fraction of the genome covered at 20X
	if ($line =~ /^total	20	.*/) {
		my ($t,$c,$target_coverage) = split(/\t/,$line);
		if ($target_coverage eq "undetermined") {
			$target_coverage = 0 ;
		} else {
			$target_cov = $target_coverage*100;
		}
	}

}
close($IN);

########################
## PARSE PANGOLIN REPORT
########################

open (my $IN, '<', $pangolin) or die "FATAL: Can't open file: $pangolin for reading.\n$!\n";

# parse pangolin report and get lineage assignment
my $global_lineage = "undetermined";
my $voc_call = undef;

chomp(my @lines = <$IN>);

my $header = shift @lines ;

foreach my $line (@lines) {

	chomp($line);

	my ($seq,$lineage,$conflict,$ambig,$scorpio_call,$scorpio_support,$scorpio_conflict,$scorpio_notes,$vers,$p_vers,$s_vers,$c_vers,$designated,$qc_status,$qc_notes,$note) = split(",", $line);

        next unless ($qc_status eq "pass");

	$global_lineage = $lineage;
	if (length $scorpio_call > 0) {
		$voc_call = $scorpio_call;
	}
}

close($IN);

$data{"Pangolin"}= {"lineage" => $global_lineage, "voc" => $voc_call } ;

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

$data{"Sars-CoV2"}= {"Status" => $status} ;

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

$data{"reads"}= {"mapped" => $mapped_reads , "coverage_20X" => $target_cov} ;

#########################
## PARSE ASSEMBLY STATS
#########################

open (my $IN, '<', $assembly_stats) or die "FATAL: Can't open file: $assembly_stats for reading.\n$!\n";

my $assembly_length = "undetermined";
my $assembly_gaps = "undetermined";
my $genome_fraction = "undetermined";

while (<$IN>) {

        chomp;
        my $line = $_;

        if ($line =~ /^Total.*/) {
		($assembly_length) = $line =~ /(\d+)/;
	} elsif ( $line =~ /^N.*/) {
		($assembly_gaps) = $line =~ /(\d+)/;
		$genome_fraction = ($assembly_gaps/$assembly_length) ;
	}
}
close($IN);

my $rounded = 100*(sprintf "%.2f", $genome_fraction);

$data{"Assembly"}= {"Anteil_Ns" => $rounded, "Laenge" => $assembly_length, "Gaps" => $assembly_gaps} ;

#########################
## PARSE SnpEff VCF FILE
#########################

open (my $IN, '<', $vcf) or die "FATAL: Can't open file: $vcf for reading.\n$!\n";

# Array of hashes to old effect predictions for each position
my @variant_data;

my @records;
my @fields;

my @table = (
	[ "Position","Annotation","Effekt","Gen Name","Transkript","HGVS_p" ]
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
	my $gene_name = $annot{'Gene_Name'};
	my $gene_id = $annot{"Feature_ID"} ;
	my $hgvs_p = $annot{"HGVS.p"} ;
	my $hgvs_c = $annot{"HGVS.c"} ;
	my $a = $annot{'Annotation'};

	my @te = [ @elements[1] . ":" . @elements[3] . ">" . @elements[4] , 
		$a,
		$effect ,
		$gene_name,
		$gene_id ,
		$hgvs_p
	] ;

	push(@table, @te);

	my %entry = (
		"genomic_position" => @elements[1],
		"ref_base" => @elements[3],
		"alt_base" => @elements[4],
		"annotation" => $a,
		"effect" => $effect,
		"gene_name" => $gene_name,
		"transcript_id" => $gene_id,
		"hgvs_c" => $hgvs_c,
		"hgvs_p" => $hgvs_p
	);
	
	push(@variant_data,\%entry);		

}

$data{"Variants"} = \@variant_data;

close($IN);


##################
## Build JSON
##################

my $json = encode_json \%data;

print $json . "\n";


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
$step -= 20;

$text->translate(50,$step);
$text->font($b_font,10);
$text->text("Pipeline version:");

$text->font($font,10);
$text->translate(250,$step);
$text->text($pipeline_version);

$step -= 30;
$text->font($b_font,10);
$text->translate(50,$step);
$text->text("Alternative ID");

$text->font($font,10);
$text->translate(250,$step);
$text->text($patient);

$step -= 20;
$text->font($b_font,10);
$text->translate(50,$step);
$text->text("Library ID");

$text->font($font,10);
$text->translate(250,$step);
$text->text($library);

$step -= 20;
$text->font($b_font,10);
$text->translate(50,$step);
$text->text("Sars-CoV2 Status:");

$text->font($font,10);
$text->translate(250,$step);
$text->text($status);

$step -= 20;
$text->font($b_font,10);
$text->translate(50,$step);
$text->text("Typ/Lineage (Pangolin):");

$text->font($font,10);
$text->translate(250,$step);
$text->text("${global_lineage}");

if (defined $voc_call) {
	$step -= 15;
	$text->font($b_font,10);
	$text->translate(50,$step);
	$text->text("Kritische Variante");
	
	$text->font($font,10);
	$text->translate(250,$step);
	$text->text("Ja (${voc_call})");
}

$step -= 30;
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
$text->text("Anteil Ns am Assembly:");

$text->font($font,10);
$text->translate(250,$step);

if ($rounded > 5) {
        $text->fillcolor('red');
}

$text->text( $rounded . "%");
$text->fillcolor('black');

$step -= 20;
$text->font($b_font,10);
$text->translate(50,$step);
$text->text("Abdeckung mit mind. 20X:");

$text->font($font,10);
$text->translate(250,$step);

if ($target_cov < 95) {
        $text->fillcolor('red');
}

$text->text($target_cov . "%");
$text->fillcolor('black');

$step -= 20;

my $gfx_plot = $page->gfx();

my $image_plot = $pdf->image_jpeg($plot);

$gfx_plot->image($image_plot,330,$step,$image_plot->width/10,$image_plot->height/10);

$step -= 10;
$text->font($font,8);;
$text->translate(330,$step);
$text->text("Abdeckung der Genomsequenz (max. 200X)");

$step -= 20;
$text->font($b_font,12);
$text->translate(50,$step);
$text->text("Beobachtete Varianten");

$text->font($font,8);

# Variant table

$step -= 20;
$text->translate(50,$step);
$text->text("Sars-CoV2 Referenz: NC_045512.2");

my $pdftable = new PDF::Table;

my $table_ref = \@table;

my $left_edge_of_table = 50;

$step -= 10;

$pdftable->table(
     $pdf,
     $page,
     $table_ref,
     'x' => $left_edge_of_table,
     'w' => 500,
     'start_y' => $step,
     'start_h' => 300,
     'next_y'          => 750,
     'next_h'          => 500,
     'padding'         => 1,
     'padding_right'   => 1,
     'font_size'       => 7,
     'background_color_odd'    => "white",
     'background_color_even'   => "lightgray", 
     'max_word_length' => 30, # 50 between forced splits
);

## Footer ##

$text->font($font,8);

$text->translate(50,75);
$text->text("IKMB pipeline version $pipeline_version - see: https://github.com/ikmb/virus-pipe | Software tools:");

$text->translate(50,65);

# Write the software versions

my $audit = "";

foreach my $k (keys %tools) {

	my $version = $tools{$k};	
	$audit .= "$k:$version "

}

$text->text($audit);

$text->translate(50,50);


$text->text("Report erstellt: $date");

$pdf->saveas($report_name);

