#!/usr/bin/env perl
#
#This script takes all orthogroups from the mysql database and writes as text file
#


#select og from orthogroup;

use strict;
use DBI;
my $dbh = DBI->connect ("DBI:mysql:ellis2021","root","9DogsinaTree");

my $og;
my $query;
my $query_handle;

	# PREPARE THE QUERY

	$query = "select DISTINCT og from orthogroup;";

	$query_handle = $dbh->prepare($query);
	$query_handle->execute();
	$query_handle->bind_columns(undef, \$og);



# LOOP THROUGH RESULTS FROM THE QUERY
	while($query_handle->fetch()) {
		chomp($og);
		print $og."\n";
	}
