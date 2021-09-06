#!/usr/bin/perl
use strict;
use warnings;

#define variables
my($path1,$vcfsfile,$id,$variants,$freq,$regid)=@ARGV;

my $Line;
my $countf=0;
my %vars=();

##open output file to print Catalogue of VARIANTs
my $output_file = $path1."/u_".$id."_V9_May2021_Variant_List_".$regid."_".$freq."_ALLvars_GRCh38.txt";

open(OUT, ">$output_file") || die "Cannot open file \"$output_file\" to write to!\n";
print OUT "PKID\tVARIANT\tGENOTYPE\tGENE\tGNOMAD\tPOPN\tCLINVAR\n";

##open variants file and read list to check into %vars
open INPUT4, $variants or die "Cannot open $variants\n";
	my $head=<INPUT4>;
        varloop: foreach my $v (<INPUT4>){
                chomp $v;
		$v=~s/\"//g;
		my @linesplit4=split(/\t/,$v);
	#if(!exists $pidgenes{$linesplit4[4]}){next varloop;}
	if(!exists $vars{$linesplit4[0]}{$linesplit4[1]}{$linesplit4[2]}{$linesplit4[3]}){
		$vars{$linesplit4[0]}{$linesplit4[1]}{$linesplit4[2]}{$linesplit4[3]}{'hets'}=0;
		$vars{$linesplit4[0]}{$linesplit4[1]}{$linesplit4[2]}{$linesplit4[3]}{'homs'}=0;
		$vars{$linesplit4[0]}{$linesplit4[1]}{$linesplit4[2]}{$linesplit4[3]}{'gene'}=$linesplit4[4];
		}
	}
close INPUT4;

##open each vcf file and add variants into %vars
open INPUT2, $vcfsfile or die "Cannot open $vcfsfile\n";
        lineloop2: foreach my $file (<INPUT2>){
                chomp $file;
		if($file!~/\/\S+\//){$file=$path1."/norm_".$regid."/".$file;}
	## test if file is gz compressed
	if($file =~/.gz$/){
		open(INPUT, "gunzip -c $file |") or die "Cannot open compressed file $file\n";
		}
	else{open INPUT, $file or die "Cannot open $file\n";}
	## loop through each vcf file line
	my $ptid="\.";
	lineloop: foreach $Line (<INPUT>){
                        chomp $Line;
				if($Line=~/^\#CHROM/){
                                my @linesplit3=split(/\t/,$Line);
				$ptid=$linesplit3[9];
				next lineloop;
                                }
				if($Line=~/^\#/){
                                next lineloop;
                                }
				
		##Get PID gene symbol for variant
		my $Gene="\.";

		##split line on tab
		my @linesplit2=split(/\t/,$Line);
		
		### Filter on FILTER == PASS
		if($linesplit2[6]!~/PASS/){next lineloop;}
		
		##Only include variants specified in imported file
		if(!exists $vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$linesplit2[4]}){
                                next lineloop;
				}#include specific variants only
		$Gene=$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$linesplit2[4]}{'gene'};
		##separate out comma separated variant alleles into @Alleles
		my @Alleles=();
		$Alleles[0]=$linesplit2[3];
		if($linesplit2[4]=~/\,/){
			my @alts=split(/\,/,$linesplit2[4]);
			push(@Alleles,@alts);
			}elsif($linesplit2[4]!~/\,/){
				$Alleles[1]=$linesplit2[4];
				}
				
		##loop through each genotype
                loop: for(my $i=9;$i<scalar(@linesplit2);$i++){
			if($linesplit2[$i]=~/^([0123456789])\/([123456789])\:/){
				my $first=$1;
				my $second=$2;
				if($first eq $second){
						my $var = $linesplit2[0]."_".$linesplit2[1]."_".$linesplit2[3]."_".$Alleles[$second];
						my $genoHOM = $Alleles[$second]."\/".$Alleles[$second];
						print OUT "$ptid\t$var\t$genoHOM\t$Gene\t\.\t\.\t\.\n";
					next loop;
					}
				if($first ne $second){
							if($first!~/0/){$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$first]}{'hets'}++;}
							my $var	= $linesplit2[0]."_".$linesplit2[1]."_".$linesplit2[3]."_".$Alleles[$second];
        	                                        my $genoHET = $linesplit2[3]."\/".$Alleles[$second];
	                                                print OUT "$ptid\t$var\t$genoHET\t$Gene\t\.\t\.\t\.\n";
					next loop;
					} #heterozygote
			}# if genotype non-REF
				#chrX or Y males with single allele genotype call
				if($linesplit2[$i]=~/^([123456789])\:/){
					my $second=$1;
						my $var = $linesplit2[0]."_".$linesplit2[1]."_".$linesplit2[3]."_".$Alleles[$second];
                                                my $genoHEM = $Alleles[$second];
                                                print OUT "$ptid\t$var\t$genoHEM\t$Gene\t\.\t\.\t\.\n";				
					next loop;
					}# if genotype non-REF chrX or Y males
		}#for loop - genotypes
        }#lineloop
        close INPUT;
  }#foreach vcf file
close INPUT2;
close OUT;

exit;
