#!/usr/bin/env perl
#
#does search and replace of sample names replacing with species names in a file

use warnings;
use strict;
use File::Slurp;

my $tabfile = $ARGV[0];
my $changefile = $ARGV[1];

open(FH, '<', $tabfile) or die $!;

my $changetxt = read_file($changefile);

my @line;
while(<FH>){
	chomp($_);
	@line = split(/\t/, $_);
	my $species = $line[0];
	my $country = $line[1];
	my $sample = $line[2];
	if($changetxt =~ m/$sample/){
#		print $species."\t".$country."\t".$sample."\n";
		$species =~ s/ /\_/g ;
		$country =~ s/ /\_/g ;
		my $fullotu = $species."__".$country;
#		$changetxt =~ s/$sample/$fullotu/g ;
		$changetxt =~ s/$sample/$species/g ;
	}else{
#		print "\t\t**No $sample\n";
	}
}

close(FH);
print "$changetxt\n";

