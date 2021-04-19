#!/usr/bin/env perl
use Bio::SeqIO;		#requires bioperl library

    my $filename = $ARGV[0];
    #assuming first part of file name is ortholog name
    my @og = split('\.', $filename);
    #strip anything before / for directory names
    $og[0] =~ s/.+\///g ;

    open(my $fh, '<:encoding(UTF-8)', $filename)
        or die "Could not open file '$filename' $!";

    my $in  = Bio::SeqIO->new(-file => $filename ,
                           -format => 'fasta');
    my $count;
    my $sampleid;
    while ( my $seq = $in->next_seq() ) {
        my $description = $seq->desc;
	my $fullid = $seq->id;
	print $fullid."\t".$og[0]."\t".$fullid."_".$og[0]."\t".$seq->seq."\n";
    }

