#!/usr/bin/env perl
#
#This script takes sequences from the mysql database and writes as fasta to screen.
#Also produces a gene map file, compatible with astral-pro
#usage ./ogsql2fasta <OG>
#


#select a.sampleid, o.assemblyid, o.og from aaseqs a, orthogroup o WHERE a.assemblyid=o.assemblyid AND o.og LIKE 'OG';

use strict;
use Bio::SeqIO;
use DBI;
my $dbh = DBI->connect ("DBI:mysql:ellis2021","root","9DogsinaTree");

my $ogsearch = $ARGV[0];
my $assemblyid="";
my $sampleid="";
my $geneid="";
my $fullseq;
my $query;
my $query_handle;
my $genename;
my $og;

#Create a hash to look up whether geneid already exists because many duplicated geneids from trinity
my %geneidLUT;

	# PREPARE THE QUERY

	$query = "select a.sampleid, o.og, a.geneid, a.aa from aaseqs a, orthogroup o WHERE a.assemblyid=o.assemblyid AND o.og LIKE '$ogsearch'; ";

	$query_handle = $dbh->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$sampleid, \$og, \$geneid, \$fullseq);

# LOOP THROUGH RESULTS FROM THE QUERY
	while($query_handle->fetch()) {
		$fullseq =~ s/\*//g;
		chomp($fullseq);

		#Assemble gene name as Sample_Ortholog_Gene
		$genename = join("_", $sampleid, $og, $geneid);
#print $genename."\n";
		if (exists $geneidLUT{$genename}) {
			print "#".$genename." is duplicated ";
			if (length($geneidLUT{$genename}) == length($fullseq) ) {
#			if (length($fullseq) > 10 ) {
				print "and equal length \n";
			}else{
				print "and unequal length ";
				if(length($geneidLUT{$genename}) > length($fullseq) ) {
					print "keeping longer\n";
				}else{
					$geneidLUT{$genename} = $fullseq;
					print " \n";
				}
			}
		}else{
			$geneidLUT{$genename} = $fullseq;
		}
	}


for(keys %geneidLUT){
	print(">$_\n$geneidLUT{$_}\n");
}
