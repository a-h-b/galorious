#!/usr/bin/env perl

use strict;
use Getopt::Long;

my ($gf_file,$out_file,$smash_files);
GetOptions('s|smash=s' => \$smash_files,
	   'a|annotationgff=s' => \$gf_file,
	   'o|out=s' => \$out_file,);

if ($gf_file eq $out_file) {
	die "Input and ouput file should not be the same!\n";
}


my %data=();
my @all = ();
my $line = ();

print STDOUT "# gff: $gf_file\n";
open(FH, "$gf_file") or die "Failed to open $gf_file\n";
@all = <FH>;
while (@all) {
        $line = shift @all;
        chomp($line);
        next if $line=~/^#/;
        my ($contig, $tool, $type, $start,$end,$point,$sense,$zero,$description) = split(/\t/, $line);
        my @descriptors = split(/;/, $description);
        my $product = $descriptors[0];
        $product =~ s/ID=//;
        $data{$product}{'tabs'} = join("\t",($contig, $tool, $type, $start,$end,$point,$sense,$zero));
	@{$data{$product}{'att'}} = @descriptors;
}
close(FH);


my @smashfiles = split(",",$smash_files);

my $ID_cnt = 1;
my %boring_cats = map { $_ => 1 } ("ID","Name","inference","locus_tag","partial","phase","product","source","transl_table","translation");

my $outfile=$out_file;
open RESULT,">$outfile";

foreach my $infile (@smashfiles){
	print STDOUT "Reading ".$infile."\n";
	my $region;
	if($infile =~ m/[0-9]+\.region[0-9]{3}/){
		$region = $1;
	}
	my $corstart = 0;
        my @neighbour;
        my $offset = 0;
	open (INFILE, "$infile") or die $!;
	while(my $line=<INFILE>){
		chomp($line);
		next if $line =~/#/;
		my @tabs = split("\t",$line);
		my $attr = pop @tabs;
		my @attrs = split(";",$attr);
		my $id = "";
		my @tmpattr = ();
		if($tabs[2] eq "protocluster"){
			my @proto_attr = ();
			foreach my $anno (@attrs){
				if($anno =~ /^ID=/){
					my @cat = split("=",$anno);
					$anno = $cat[0]."=".$region =~ s/[.]/_/r;
				}elsif($anno =~ /^core_location/){
					my @corloc=split("=",$anno);
					$corstart = $corloc[1] =~ s/\[([0-9]+):([0-9]+)]/$1/r;
				}elsif($anno =~ /^neighbourhood/){
					@neighbour = split("=",$anno);
				}
				push @proto_attr,$anno; 
			}
			$tabs[8] = join(";",@proto_attr);
			$offset = $corstart-$neighbour[1];
		}elsif($tabs[2] ne "CDS"){
			my @proto_attr = ();
                       	foreach my $anno (@attrs){
				if($anno =~ /^ID=/){
                               	        $anno = "ID=".$region."_".sprintf("%04d",$ID_cnt) =~ s/[.]/_/r;
					$ID_cnt += 1;
                               	}
				push @proto_attr,$anno;
			}
			$tabs[8] = join(";",@proto_attr);
		}
		if($tabs[2] ne "CDS"){
			$tabs[1] = "antiSMASH";
			$tabs[3] = $tabs[3]+$offset;
			$tabs[4] = $tabs[4]+$offset;
			print RESULT join("\t",@tabs)."\n";			
		}else{
			my @proto_attr = ();
			my @ID;
			foreach my $anno (@attrs){
				my @cat=split("=",$anno);
				if($cat[0] eq "ID"){
					@ID=split("=",$anno);
				}
				unless(exists $boring_cats{$cat[0]}){
					push @proto_attr,$anno;
					print STDOUT $anno."\n";
				}
			}
			if(@proto_attr){
				push @{$data{$ID[1]}{'att'}},@proto_attr;		
			}
		}
	}
	close(INFILE);
}

	
foreach my $id (keys(%data)){
	my $tabs = $data{$id}{'tabs'};
	my @unique = do { my %seen; grep { !$seen{$_}++ } @{$data{$id}{'att'}} };
	print RESULT $tabs."\t".join(";", @unique)."\n";
}

close(RESULT)



