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
	@otu = split(/OG/, $_);
	chop($otu[0]);	#remove trailing _
	print $_."\t".$otu[0]."\n";
}

close(FH);


