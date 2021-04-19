#!/usr/bin/env perl

use strict;
use warnings;

my $filename = $ARGV[0];
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";

system "rm -rf clocktrees";
system "mkdir clocktrees";

system "rm -rf clockgenes";
system "mkdir clockgenes";

while (my $row = <$fh>) {
        chomp $row;
	my @field = split(' ', $row);
	unless($field[0] eq 'name'){
		my @og = split('\.', $field[0]);
#		print $og[0]."\n";
		system "cp sortadate/$field[0] clocktrees/";
		system "cp sortadate/$og[0].raln.fa clockgenes/";
		system "cp sortadate/$og[0].pppft.treefile clocktrees/";
	}
}


