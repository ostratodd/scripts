#!/usr/bin/env perl
use Bio::SeqIO;		#requires bioperl library


    my $filename = $ARGV[0];
    open(my $fh, '<:encoding(UTF-8)', $filename)
        or die "Could not open file '$filename' $!";

    my $in  = Bio::SeqIO->new(-file => $filename ,
                           -format => 'fasta');
    my $count;
    my $sampleid;
    while ( my $seq = $in->next_seq() ) {
	print $seq->id."\t";
	print "concatmatrix"."\t";
	print "matrix".$seq->id."\t";
	my $sequence = $seq->seq;
#next line keeps only substring for easier debugging comment to remove and print full sequence
#	$sequence = substr($sequence, 4,100);
	print $sequence."\n";
    }

