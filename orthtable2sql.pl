#!/usr/bin/env perl

use strict;
use warnings;
use DBI;

my $filename = $ARGV[0];
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";


my $dbh = DBI->connect ("DBI:mysql:ellis2021","root","9DogsinaTree");
while (my $row = <$fh>) {
        chomp $row;
        if ($row eq ""){
                exit;
        }
        my @field = split /\t/, $row;
	if ($field[0] eq 'Orthogroup') { 	#still on first line with header 
	}else{
		foreach ( my $i = 0; $i < @field; $i++) {
			if($i > 0) {
				$field[$i] =~ s/\t//g;
				if($field[$i] =~ m/\,/) {
					chomp($field[$i]);
				        my @paralog = split /\,/, $field[$i];
					foreach ( my $j = 0; $j < @paralog; $j++) {
						$paralog[$j] =~ s/\t//;
						$paralog[$j] =~ s/ //;


					        $dbh->do("insert into orthogroup (og,assemblyid)
                					values
               					 ('$field[0]','$paralog[$j]' )");

#						print $field[0]."\t".$paralog[$j]."\n";
					}
				}elsif($field[$i] =~ m/\_/) {	#all non-empty names will have _ so only print non-empty

					$dbh->do("insert into orthogroup (og,assemblyid)
                				values
               			 	('$field[0]','$field[$i]' )");
#					print $field[0]."\t".$field[$i]."\n";
				}
			}
		}
	}
}


