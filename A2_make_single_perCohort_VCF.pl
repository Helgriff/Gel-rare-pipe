#!/usr/bin/perl
use strict;
use warnings;

#define variables
my $Line;
my %vars=();

my($path1,$vcfsfile,$id)=@ARGV;
my $countf=0;

##open each vcf file and add variants into %vars
open INPUT2, $vcfsfile or die "Cannot open $vcfsfile\n";
        lineloop2: foreach my $file (<INPUT2>){
                chomp $file;
	## test if file is gz compressed
	if($file =~/.gz$/){
		open(INPUT, "gunzip -c $file |") or die "Cannot open compressed file $file\n";
		}
	else{open INPUT, $file or die "Cannot open $file\n";}
	## loop through each vcf file line
	lineloop: foreach $Line (<INPUT>){
                        chomp $Line;
			my $id="\.";
				if($Line=~/^\#CHROM/){
                                my @linesplit3=split(/\t/,$Line);
				$id=$linesplit3[9];
				next lineloop;
                                }
				if($Line=~/^\#/){
                                next lineloop;
                                }
                my @linesplit2=split(/\t/,$Line);
		##separate out comma separated variant alleles into @Alleles
		my @Alleles=();
		$Alleles[0]=$linesplit2[3];
		if($linesplit2[4]=~/\,/){
			my @alts=split(/\,/,$linesplit2[4]);
			push(@Alleles,@alts);
			}elsif($linesplit2[4]!~/\,/){
				$Alleles[1]=$linesplit2[4];
				}
		for(my $a=0;$a<scalar(@Alleles);$a++){		#for each REF and ALT allele
                        if(!exists $vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$a]}){
                                $vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$a]}{'total'}=0;
				$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$a]}{'alleles'}=0;
                                $vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$a]}{'hets'}=0;
                                $vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$a]}{'homs'}=0;
				}#initialise each variant
			}#for each alt allele
                
		##loop through each genotype
                loop: for(my $i=9;$i<scalar(@linesplit2);$i++){
			   for(my $a=0;$a<scalar(@Alleles);$a++){
				$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$a]}{'total'}+=2;
				} #for each alt allele
			if($linesplit2[$i]=~/^0\/0\:/ or $linesplit2[$i]=~/^\.\/\.\:/ or $linesplit2[$i]=~/^0\:/ or $linesplit2[$i]=~/^\.\:/){
				$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[0]}{'alleles'}+=2;
				$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[0]}{'homs'}++;
				next loop;
				} #increment total called genotypes (N.B chrX and Y male non-calls appear to be coded ./.)
			if($linesplit2[$i]=~/^([0123456789])\/([123456789])\:/){
				my $first=$1;
				my $second=$2;
				if($first eq $second){
					$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$second]}{'alleles'}+=2;
					$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$second]}{'homs'}++;
					next loop;
					} #homozygote for variant ... ref genos 0/0 'looped' out line 60
			if($first ne $second){
					$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$first]}{'alleles'}++;
					$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$first]}{'hets'}++;
					$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$second]}{'alleles'}++;
					$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$second]}{'hets'}++;
					next loop;
					} #heterozygote
			}# if genotype non-REF
				#chrX or Y males with single allele genotype call
				if($linesplit2[$i]=~/^([123456789])\:/){
					my $second=$1;
					$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$second]}{'alleles'}+=2;
					$vars{$linesplit2[0]}{$linesplit2[1]}{$linesplit2[3]}{$Alleles[$second]}{'homs'}++;
					next loop;
					}# if genotype non-REF chrX or Y males
		}#for loop - genotypes
        }#lineloop
        print "$file\n";
	$countf++;
        close INPUT;
  }#foreach vcf file
 close INPUT2;

##open output file and print header
my $output_file = $path1."/InHouse_MAFs_".$id."_x".$countf.".vcf";
open(OUT, ">$output_file") || die "Cannot open file \"$output_file\" to write to!\n";
print OUT "\#\#fileformat\=VCFv4\.1\n\#\#FILTER\=\<ID\=PASS\,Description\=\"All filters passed\"\>\n";
print OUT "\#\#contig\=\<ID\=chr1\,length\=248956422\>\n";
print OUT "\#\#contig\=\<ID\=chr2\,length\=242193529\>\n";
print OUT "\#\#contig\=\<ID\=chr3\,length\=198295559\>\n";
print OUT "\#\#contig\=\<ID\=chr4\,length\=190214555\>\n";
print OUT "\#\#contig\=\<ID\=chr5\,length\=181538259\>\n";
print OUT "\#\#contig\=\<ID\=chr6\,length\=170805979\>\n";
print OUT "\#\#contig\=\<ID\=chr7\,length\=159345973\>\n";
print OUT "\#\#contig\=\<ID\=chr8\,length\=145138636\>\n";
print OUT "\#\#contig\=\<ID\=chr9\,length\=138394717\>\n";
print OUT "\#\#contig\=\<ID\=chr10\,length\=133797422\>\n";
print OUT "\#\#contig\=\<ID\=chr11\,length\=135086622\>\n";
print OUT "\#\#contig\=\<ID\=chr12\,length\=133275309\>\n";
print OUT "\#\#contig\=\<ID\=chr13\,length\=114364328\>\n";
print OUT "\#\#contig\=\<ID\=chr14\,length\=107043718\>\n";
print OUT "\#\#contig\=\<ID\=chr15\,length\=101991189\>\n";
print OUT "\#\#contig\=\<ID\=chr16\,length\=90338345\>\n";
print OUT "\#\#contig\=\<ID\=chr17\,length\=83257441\>\n";
print OUT "\#\#contig\=\<ID\=chr18\,length\=80373285\>\n";
print OUT "\#\#contig\=\<ID\=chr19\,length\=58617616\>\n";
print OUT "\#\#contig\=\<ID\=chr20\,length\=64444167\>\n";
print OUT "\#\#contig\=\<ID\=chr21\,length\=46709983\>\n";
print OUT "\#\#contig\=\<ID\=chr22\,length\=50818468\>\n";
print OUT "\#\#contig\=\<ID\=chrX\,length\=156040895\>\n";
print OUT "\#\#contig\=\<ID\=chrY\,length\=57227415\>\n";
print OUT "\#\#contig\=\<ID\=chrM\,length\=16569\>\n";
print OUT "\#\#INFO\=\<ID\=TOTC\,Number\=A\,Type\=Float\,Description\=\"Total called genotypes that are heterozygous or alternative allele homozygous \"\>\n";
print OUT "\#\#INFO\=\<ID\=HETS\,Number\=A\,Type\=Float\,Description\=\" Number of heterozygous genotype calls\"\>\n";
print OUT "\#\#INFO\=\<ID\=HOMS\,Number\=A\,Type\=Float\,Description\=\"  Number of alternative allele homozygous genotype calls\"\>\n";
print OUT "\#\#INFO\=\<ID\=MAF\,Number\=A\,Type\=Float\,Description\=\" Frequency of alternative allele in all alternative allele gentype calls\"\>\n";
print OUT "\#\#FORMAT\=\<ID\=GT\,Number\=1\,Type=String\,Description\=\"Genotype\"\>\n";
print OUT "\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_ID1\n";

foreach my $c (sort keys %vars){
        foreach my $p (sort {$a<=>$b} keys %{$vars{$c}}){
                foreach my $r (sort keys %{$vars{$c}{$p}}){
		      vloop: foreach my $v (sort keys %{$vars{$c}{$p}{$r}}){
				if($r =~/^$v$/){next vloop;}
				my $maf="\.";
				my $totalgenotypes=0;
				unless($vars{$c}{$p}{$r}{$v}{'total'}==0){
					$maf=$vars{$c}{$p}{$r}{$v}{'alleles'}/$vars{$c}{$p}{$r}{$v}{'total'};
					$totalgenotypes=$vars{$c}{$p}{$r}{$v}{'total'}/2;
					}#unless total=0
				print OUT "$c\t$p\t\.\t$r\t$v\t250\tPASS\tTOTC\=$totalgenotypes\;HETS\=$vars{$c}{$p}{$r}{$v}{'hets'}\;HOMS\=$vars{$c}{$p}{$r}{$v}{'homs'}\;MAF\=$maf\tGT\t0\/1\n";
				}
			}
		}
	}
close OUT;

exit;
