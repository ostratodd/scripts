#! /usr/bin/perl -w

use strict;
use warnings;
#phylocatenator.pl
#Written by Todd H. Oakley UCSB
#Version 1.0 May, 2012
#
#phyloconcatenator.pl [infile] [also requires arguments 1-9 detailed below]
#
#Version 1.0.1 Sept 2012 -- fixed a bug that led sometimes to species with no data being retained.
#This bug would not affect the data that was retained, so if raxml ran all was fine. However, 
#raxml will not run with species that have completely missing data. Added extra section to remove 
#those species.

#For debugging command line pass, uncomment next
#for(my $i=0; $i < @ARGV; $i++){
#	print "Arg $i:  $ARGV[$i] \n\n";
#}
#exit;

die "Check arguments" unless @ARGV == 10;

#This was created for use in Galaxy so all 10 arguments must be passed in order
#Obtain Arguments
my $infile1=shift(@ARGV);		#0 input file
my $mingf_persp=shift(@ARGV);		#1 minimum number of genefamiles to keep species
my $minsp_pergf=shift(@ARGV);		#2 minimum number of species to keep genefamily
my $min_gf_len =shift(@ARGV);		#3 minimum gene family length
my $speciesfile=shift(@ARGV);		#4 optional species file 'None' if false
my $modelsfile=shift(@ARGV);		#5 models file
my $outfile=shift(@ARGV);		#6 data outfile
my $partfile=shift(@ARGV);		#7 partition file
my $htmlfile=shift(@ARGV);		#8 html file
my $htmlfile2=shift(@ARGV);		#9 html file for accession numbers

#Open 2 output files
open (OUT, ">$outfile") or die "Cannot create $outfile \n";
open (PART, ">$partfile") or die "Cannot create $partfile\n";
open (TABLE, ">$htmlfile") or die "Cannot create $htmlfile\n";
open (TABLE2, ">$htmlfile2") or die "Cannot create $htmlfile2\n";

my %HashAccession;
my %HoAdatatable;
my %HoGF;
my %lengthof;
my %ModelFor;
my $line = "";

my @specieslist;
my @genelist;
my @nsp;
my $nsl=1;
my @uniquemodels;
my @nullmodels;
my @retainedspecies;

# read file with species
unless($speciesfile eq 'None'){
	open SPFILE, $speciesfile or die "ERROR: Cannot open file $speciesfile\n\n";
	$nsl=0;
while (<SPFILE>) {
        chomp;
	my $currentinput = "$_";
        if($currentinput =~m /\w/){     #must have a word otherwise empty
                push(@nsp, $currentinput);
        }else{
                die "ERROR: Fewer than 4 species meet your specified criteria.\n";
        }
#	if(@nsp < 4){
#		die "ERROR: Species file must have more than 4 species.\n";
#	}
}
#print "\n\n";
} #end of unless

#Determine models used for each partition, make hash
unless($modelsfile eq 'None'){
	open MODFILE, $modelsfile or die "ERROR: Cannot open file $modelsfile\n\n";
	while (<MODFILE>) {
	        chomp;
		my $currentinput = "$_";
	        if($currentinput =~m /\t/){     #must have a tab otherwise wrong file format
			if($currentinput =~m /\t\t/){
				print OUT "ERROR: file contains 2 tabs in a row.  Check phytab format.\n";
				die;
			}else{
				my @genemodel = split(/\t/, $currentinput);
				my $genefamily=$genemodel[0];
				my $curmodel = $genemodel[1];
				if (exists $ModelFor{$genefamily}) {
					print OUT "ERROR: Model specification for $genefamily is duplicated\n";
					die;
				}else{
					$ModelFor{$genefamily}=$curmodel;
					push(@uniquemodels,$curmodel);
				}
			}
		}else{
			die "ERROR: Model LUT must be genefamily\tmodel and contain no blank lines\n";
		}
	}
}
@uniquemodels = uniq(@uniquemodels); 	#remove redundant models - uniq is subroutine at end of script

#Now check that models are valid raxml models NOT DONE YET
checkraxmlmodel();

#INPUT ALL DATA
open(INFILE1, $infile1);
foreach $line(<INFILE1>) {
	chop($line);
	my $getline = $line;
	my @column = split(/\t/, $getline);
	my $species = $column[0];
	my $genefamily = $column[1];
	my $genename = $column[2];
	my $sequence = $column[3];

	if (exists $HoAdatatable{$species}{$genefamily}) {
		print OUT "ERROR: $species $genefamily is duplicated\n";
		die;
	}else{
		$HashAccession{$species}{$genefamily} = $genename;
		$HoAdatatable{$species}{$genefamily} = $sequence;
		$HoGF{$genefamily}{$species}=$sequence;
	}
}

#If no models file selected, set every gene family to GTR model
if($modelsfile eq 'None'){
	foreach my $gfkey (sort keys %HoGF){
		$ModelFor{$gfkey} = 'GTR';
	}
	push(@uniquemodels, 'GTR');
}

#First, keep all species with enough total partitions present
foreach my $specieskey (sort keys %HoAdatatable)
{
	#Count species with minimum gfs
	my $ngf_persp=0;
	foreach my $genefamilykey (sort keys %{$HoAdatatable{$specieskey}})
	{
		$ngf_persp++;
	}
	unless($ngf_persp < $mingf_persp){	#too few genes for this species, delete the sp
		if($nsl){	#No species list supplied, push all species into list
			push(@specieslist,$specieskey);
		}else{
			#See if current specieskey is in inputted species list
			my $nsp;
			foreach $nsp(@nsp){
        			if (index($nsp,$specieskey) ge 0){
					push(@specieslist,$specieskey);
				}
        		}
		}
	} 
}
print OUT "\n";
unless (@specieslist){
	print OUT "ERROR: No species with more than $mingf_persp genes\n";
	die;
}
if(@specieslist < 4){
	print OUT "ERROR: Less than 4 species with more than $mingf_persp genes\n";
	die;
}


my $oldgenelen = 0;
my $currentseqlen = 0;
foreach my $gfkey (sort keys %HoGF)
{
	$oldgenelen=0;
	#Count gfs  with minimum species
	my $nsp_pergf=0;
	for(my $j=0; $j<@specieslist; $j++){
		if(exists $HoGF{$gfkey}{$specieslist[$j]}){
			$nsp_pergf++ ; ###if exists $HoGF{$gfkey}{$specieslist[$j]};


	##get length of gene and check it is consistent
			if(exists $lengthof{$gfkey}){
				$currentseqlen = length($HoGF{$gfkey}{$specieslist[$j]});
			#	$lengthof{$gfkey} = $currentseqlen;
				if($currentseqlen == $oldgenelen){
					#$oldgenelen = $lengthof{$gfkey};
					#okay
				}else{
					print OUT "ERROR: $specieslist[$j] $gfkey sequences ". 
						"different lengths than previous. Sequences must be aligned. If ".
						"sequences are aligned, check that the line ". 
						"does not have an extra data column before the ". 
						"sequence.\n\n";
					die "ERROR: $specieslist[$j] $gfkey sequences ".
						"different lengths than previous. Sequences must be aligned. If ".
						"sequences are aligned, check that the line ". 
						"does not have an extra data column before the ". 
						"sequence.\n\n";
				}
			}else{
				$currentseqlen = length($HoGF{$gfkey}{$specieslist[$j]});
				if($currentseqlen == 0){
					die "ERROR: Zero length sequence in file.\n"
				}
				$lengthof{$gfkey} = $currentseqlen;
				$oldgenelen = $lengthof{$gfkey};
			}
		}
	}
	if($nsp_pergf < $minsp_pergf){
		#too few species for this gene family
	}else{
		if(exists $lengthof{$gfkey}){
			if($lengthof{$gfkey}>$min_gf_len){
				push(@genelist,$gfkey);
			}
		}
	}
}
if (@genelist==0){
	print OUT "ERROR: No gene families/partitions meet the specified criteria\n";
	die "ERROR: No gene families/partitions meet the specified criteria\n";
}

#Now must delete species that lack any *retained* partitions, ie some species may be all missing
foreach(@specieslist)
{
	my $curspecies=$_;
	#Count species with minimum gfs
	my $ngf_persp=0;
	my $nonmissing=0;
	foreach (@genelist)
	{
		if (exists $HoAdatatable{$curspecies}{$_}){
			my $cursequence = $HoAdatatable{$curspecies}{$_};
			$cursequence =~ s/\?//g;
			$cursequence =~ s/\-//g;
			#Will not remove N's assumes those are not missing data. Revisit?
			my $curlen = length($cursequence);
			$nonmissing = $nonmissing + $curlen;
		}
	}
	if($nonmissing > 0)	#must remove species as it contains none of the genes retained
	{
		print "$curspecies has $nonmissing Non-Missing characters\n";
		push(@retainedspecies, $curspecies);
	}
}
@specieslist=@retainedspecies;

#Remove blankspecies from
#print phylip file
#calculate n characters
my $nchar=0;	#total characters
my $ncharPs =0;	#start count for characters in current partition
my $ncharPe =0;	#end count for characters in current partition

#First count total characters for first line
for(my $k=0; $k<@genelist; $k++){
	if (exists $lengthof{$genelist[$k]}){			#check that length was calculated
		$nchar = $nchar + $lengthof{$genelist[$k]};
	}else{
		die "ERROR: $genelist[$k] LENGTH MISSING\n";
	}
}
print OUT @specieslist." ".$nchar."\n";
htmlheader();

#Need to determine gene list order, which will change due to partitioning
#then write header line of gene names, hopefully in correct order
#print TABLE "<td bgcolor=white></td>"; #Blank line in species column
print TABLE "<td style'width:2pc'><font size='-3'>Partition:</td>";
print TABLE2 "<td style'width:2pc'><font size='-3'>Partition:</td>";
for(my $part=0; $part < @uniquemodels; $part++){
	for(my $k=0; $k<@genelist; $k++){
		#First check if current gf matches current partition
		if(exists $ModelFor{$genelist[$k]}){
			if($ModelFor{$genelist[$k]} eq $uniquemodels[$part]){
				if($ModelFor{$genelist[$k]} eq $uniquemodels[$part]){
					print TABLE "<td style'width:2pc'><font size='-3'>$genelist[$k]</td>";
					print TABLE2 "<td style'width:2pc'><font size='-3'>$genelist[$k]</td>";
				}
			}
		}
	}
}
print TABLE "</tr>";
print TABLE2 "</tr>";

#print TABLE "<td bgcolor=white></td>"; #Blank line in species column
print TABLE "<td style'width:2pc'><font size='-3'>Model:</td>";
print TABLE2 "<td style'width:2pc'><font size='-3'>Model:</td>";
for(my $part=0; $part < @uniquemodels; $part++){
	for(my $k=0; $k<@genelist; $k++){
		#First check if current gf matches current partition
		if(exists $ModelFor{$genelist[$k]}){
			if($ModelFor{$genelist[$k]} eq $uniquemodels[$part]){
				if($ModelFor{$genelist[$k]} eq $uniquemodels[$part]){
					print TABLE "<td style'width:2pc'><font size='-3'>$uniquemodels[$part]</td>";
					print TABLE2 "<td style'width:2pc'><font size='-3'>$uniquemodels[$part]</td>";
				}
			}
		}
	}
}
#End of htmlheader printing

for(my $j=0; $j<@specieslist; $j++){
	print OUT "$specieslist[$j]\t";
	print TABLE "
		<tr>
		<td>$specieslist[$j]</td>";
	print TABLE2 "
		<tr>
		<td>$specieslist[$j]</td>";
	for(my $part=0; $part < @uniquemodels; $part++){
		for(my $k=0; $k<@genelist; $k++){
			#First check if current gf matches current partition
			if(exists $ModelFor{$genelist[$k]}){
				if($ModelFor{$genelist[$k]} eq $uniquemodels[$part]){
					if (exists $HoAdatatable{$specieslist[$j]}{$genelist[$k]}){
						print OUT $HoAdatatable{$specieslist[$j]}{$genelist[$k]};
						print TABLE "<td bgcolor=black></td>";
						my $acc = $HashAccession{$specieslist[$j]}{$genelist[$k]};
						print TABLE2 "<td bgcolor=lightgray>".$acc."</td>";
						$ncharPe = $ncharPe + $lengthof{$genelist[$k]};
					}else{
						if (exists $lengthof{$genelist[$k]}){
							for(my $gap=0; $gap<$lengthof{$genelist[$k]}; $gap++){
								print OUT "?";
							}
							$ncharPe = $ncharPe + $lengthof{$genelist[$k]};
							print TABLE "<td bgcolor=white></td>";
							print TABLE2 "<td bgcolor=white></td>";
							#print TABLE "<td style'width:2pc'><font size='-2'>$genelist[$k]</td>";
						}else{
							die "ERROR: BUG!! $genelist[$k] LENGTH MISSING\n";
						}
					}
				}
			}else{
				die "ERROR: $genelist[$k] is not assigned a model.  Check model LUT input.\n";
			}

		}
		if($j==0){		#print partitions first time through gene lists
			if(($ncharPe==0) || ($ncharPs > $ncharPe)){
				push(@nullmodels, $uniquemodels[$part]);
				#print "NOTE: no partitions under the model $uniquemodels[$part] made final dataset\n";
			}else{
				print PART "$uniquemodels[$part], $uniquemodels[$part] = $ncharPs - $ncharPe \n";
			}
			$ncharPs=$ncharPe + 1;
		}
	}
	print TABLE "</tr>";
	print TABLE2 "</tr>";
	print OUT "\n";
}
if($ncharPe/@genelist != $nchar){
	# Have to account for multiple times through gene lists -- ncharPe is summed multiple times but printed once
	#die "ERROR:  BUG!! Last partition number doesn't match total n characters\n";
}

#print Statistics to screen can be redirected for Log File
print "\nSPECIES:\n";
unless($speciesfile eq 'None'){
	print "Used species file to select species. \n";
}
print "Number of species with $mingf_persp or more genefamilies/partitions: ".@specieslist."\n";
print "Species list: @specieslist\n\n\n";
print "\nPARTITIONS/GENE FAMILIES:\n";
print "Number genefamilies/partitions longer than $min_gf_len characters and present in at least $minsp_pergf genefamilies/partitions: ".@genelist."\n";
for(my $i=0; $i < @genelist; $i++){
	print "$genelist[$i] Length:$lengthof{$genelist[$i]}\n";
}

	#Printing Model stats
print "\nMODELS:\n";
if($modelsfile eq 'None'){
	print "All partitions set to GTR (no model LUT supplied)\n\n";
}else{
	print "A LUT of models was supplied.\n";
	print "The following models were present in the model LUT file:\n";
	print join(" ", @uniquemodels), "\n";
	print "\n";
	for(my $i=0; $i < @nullmodels; $i++){
		print "NOTE: no partitions under the model $nullmodels[$i] made final dataset\n";
	}
}
close PART;
close OUT;
close SPFILE;
close TABLE;
close TABLE2;

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

sub checkraxmlmodel {
	my @raxmlmodels = ("DNA","BIN","MULTI","DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG", "GTR", "MTART", "MTZOA", "FLU","PMB", "HIVB","HIVW","JTTDCMUT");
	return(1);	#Not yet implemented
}
sub htmlheader {
	print TABLE '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" 
		"http://www.w3.org/TR/html4/loose.dtd">
		<html>
		<head>
		<title>Dataset Presences and Absences</title>
		<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
		</head>
		<hr><table border="1">
		<tr>
			<th align="center">Species<br>';
	print TABLE2 '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" 
		"http://www.w3.org/TR/html4/loose.dtd">
		<html>
		<head>
		<title>Accession Numbers</title>
		<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
		</head>
		<hr><table border="1">
		<tr>
			<th align="center">Species<br>';
	for(my $i=1; $i < @genelist+1; $i++){
#Genes are printing in wrong order
#		print TABLE "<td style'width:2pc'><font size='-2'>$genelist[$i]</td>";
		print TABLE "<td style'width:2pc'><font size='-2'>$i</td>";
		print TABLE2 "<td style'width:2pc'><font size='-2'>$i</td>";
	}
	print TABLE "</th>
		     </tr>";
	print TABLE2 "</th>
		     </tr>";
}
sub htmltail {
	print TABLE '
		</body>
		</html>';
	print TABLE2 '
		</body>
		</html>';
}

