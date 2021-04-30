#!/usr/bin/env perl
#
#This script creates a map file for astral, assuming a text list
#of gene names in the form of SpCode_OGnumber_DNnumber


use warnings;
use strict;

my $filename = $ARGV[0];

open(FH, '<', $filename) or die $!;

my @otu;
while(<FH>){
	chomp($_);
	@otu = split(/\s+/, $_);
	if($otu[0] =~ /\d/) {
		#OTU is digit so first line containing ntax and nchar
	}else{
		print $otu[0]."\t";
		print "morphology\t";
		print $otu[0]."morph\t";
		print $otu[1]."\n";
	}
}

close(FH);


