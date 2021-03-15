#!/usr/bin/env perl
use Bio::SeqIO;		#requires bioperl library
use DBI;		#requires DBI to interface with mysql
#Assumes TRINITY assembly, followed by transdecoder; with sample/transcriptome name appended to beginning of each fasta header
#
#Relies on mysql table configured as follows
#
# CREATE TABLE aaseqs
# (
# id int unsigned not null auto_increment primary key,
# sampleid varchar(100),
# geneid varchar(50),
# assemblyid varchar(300),
# cluster varchar(5),
# gene varchar(4),
# isoform varchar(4),
# gstatus varchar(50),
# aa text
# );

    my $filename = $ARGV[0];
    open(my $fh, '<:encoding(UTF-8)', $filename)
        or die "Could not open file '$filename' $!";

    my $in  = Bio::SeqIO->new(-file => $filename ,
                           -format => 'fasta');
    my $count;
    my $sampleid;
    while ( my $seq = $in->next_seq() ) {
        my $description = $seq->desc;
	$description =~ s/\:\:/\:/g ;	#remove double ::
        my @descs = split(/:/, $description);


        my $gstatusString = $descs[1];
        my $gstatus = $gstatusString;
        $gstatus =~ s/ len//;   #remove ' len' at end of string don't need length as sql can calculate

        #Parse ID, cluster, gene, and isoform numbers
	my $fullid = $seq->id;
	my @sid = split(/\|/, $fullid);
	$sampleid = $sid[0];

        my @ids = split(/\|/, $fullid);
        my $geneid=$ids[1];

	#Not transcriptome so these do not exist
        my $cluster="";
        my $gene = "";
        my $isoform = "";

        my $aaseq = $seq->seq."\n";
        my $assemblyid = $fullid;
	$assemblyid =~ s/\:\:/\_\_/g ;
	$assemblyid =~ s/ //g ;
	$assemblyid =~ s/\|/\_/g ;
	$assemblyid =~ s/\./\_/g ;

	#into sql
        my $dbh = DBI->connect ("DBI:mysql:ellis2021","root","9DogsinaTree");
        $dbh->do("insert into aaseqs (sampleid, geneid, assemblyid, cluster, gene, isoform, gstatus, aa)
                values
                 ('$sampleid', '$geneid', '$assemblyid', '$cluster', '$gene', '$isoform', '$gstatus', '$aaseq')");
#print "SID: $sampleid GID:$geneid ASID: $assemblyid CLU: $cluster GENE:$gene ISO:$isoform GST:$gstatus \nAA:$aaseq\n";
        $count++;
    }
    print "\n\tADDED $count sequences from $sampleid to mysql in ellis2021 database, aaseqs table\n";

