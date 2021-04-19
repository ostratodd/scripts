#!/usr/bin/env perl
use Bio::SeqIO;		#requires bioperl library
#get the length of the longest sequence in a fasta file
    my $filename = $ARGV[0];
    open(my $fh, '<:encoding(UTF-8)', $filename)
        or die "Could not open file '$filename' $!";

    my $in  = Bio::SeqIO->new(-file => $filename ,
                           -format => 'fasta');
    while ( my $seq = $in->next_seq() ) {
	my $curlen = length($seq->seq);
	if($curlen > $maxlen) {
		$maxlen = $curlen;
		$curlen = 0;
	}
    }
    print $maxlen."\n";
