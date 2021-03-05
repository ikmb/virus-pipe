#!/usr/bin/env perl 

use strict;
use warnings;
use REST::Client;
use Time::HiRes;
use JSON;
use JSON::Parse 'parse_json';
use Cwd;
use Data::Dumper;

my $server = 'http://172.21.99.123/IBDBase_Interface/api';
my $client = REST::Client->new();

my $global_headers = { 'Content-Type' => 'application/json' };

my $dir = getcwd;

printf "SENDING_LAB;DATE_DRAW;SEQ_TYPE;SEQ_REASON;SAMPLE_TYPE;PUBLICATION_STATUS;OWN_FASTA_ID\n";

foreach my $file (glob("$dir/*.fasta")) {

	my $file_name = (split "/",$file)[-1];

	my $knumber = (split ".fasta", $file_name)[0];	

	my $url = $server . "/get_diagnostic_job_date";
	my %data = ( "knumber" => $knumber );
 	my $json = encode_json \%data;
 
	my $ret = $client->request("get",$url,$json,$global_headers);
	my $r = parse_json($ret->responseContent() );
	my $date_full = $r->{'order_date'};
	my $date = (split "T", $date_full)[0];
	$date =~ s/-//g ;
	
	printf "10337" . ";" . $date . ";" . "ILLUMINA" . ";" . "X" . ";" . "X" . ";" . "P" . ";" . $knumber . "\n" ;
}


