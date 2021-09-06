#!/usr/bin/perl
use strict;
use warnings;

##Script to Create text outputs of rare variants per GENE combining participant genotypes and annotation.
##Second output includes counts of rare variants per family (proband and parents when available), per gene, allowing filtering for potential compound heterozygous genotypes.

my $path1="/file/path/to/analysis/directory";
my $path2="/file/path/to/analysis/directory/Reports";

## output files ########################################################
my $outfile2=$path2."/List_of_UTR_Exonic-splice-and-Intronic_variants_with_PKids_etc_NCFBr-probands_MAF-001_and_NotHOMinControls.txt";
my $outfile3=$path2."/List_of_UTR_Exonic-splice-and-Intronic_HetORHom_PROBAND-FAMILY_VARIANT-COUNTS_MAF-001_and_NotHOMinControls.txt";
## input files #########################################################
my $file="formatted_Labkey_Participant_Info_BRONCHplus_V9_unique_PartIds_and_GMCs_preNEW.txt"; ##family proband info file (including 'SOLVED status' column)
##
my $file2="u_ALL_NCFB_proband_relative_PCDx92_V9_May2021_Variant_List_MAF-noFilt_ALLvars_GRCh37_and_GRCh38.txt"; ## combined text files of variants and participant genotypes
##
my $file5a="ALLunfiltered_Counts_perVariant_noMAF-filter_NonBronchRels_NCFB_BRR_combined_PCDx92_GRCh38_allCHROMS.txt"; ## InHouse MAFs GRCh38 combined files from R script
my $file5b="ALLunfiltered_Counts_perVariant_noMAF-filter_NonBronchRels_Bronc-probands_Bronc-relatives_combined_PCDx92_GRCh37_allCHROMS.txt"; ## InHouse MAFs GRCh37
##
my $file7a="ALLunfiltered_Counts_perVariant_noMAF-filter_NonBronchRels_Bronc-probands_Bronc-relatives_combined_PCDx92_GRCh38_extended-Info.txt"; ##annotation for variants GRCh38 combined files from R script
my $file7b="ALLunfiltered_Counts_perVariant_noMAF-filter_NonBronchRels_Bronc-probands_Bronc-relatives_combined_PCDx92_GRCh37_extended-Info.txt"; ##annotation for variants GRCh37
##################################################################

# my $panelonly="yes"; ###yes if only Panel genes required, otherwise ="no" or "anything else other than yes"!
my $panelonly="no";

my $omimfile=$path1."/info_files/Ensembl_OMIM_GOterms_Combined_rel94_GRCh38.txt"; ##GRCh38
my $omimfileb=$path1."/info_files/Ensembl_OMIM_GOterms_Combined_rel94_GRCh37.txt"; ##GRCh37
my $PAgenefile=$path1."/info_files/PCDx92_GeneSymbols.csv"; ##GRCh38
my $PAgenefileb=$path1."/info_files/PCDx92_GeneSymbols_GRCh37.csv"; ##GRCh37

my $infile=$path1."/info_files/".$file;
my $infile2=$path1."/".$file2;
my $infile5a=$path1."/".$file5a;
my $infile5b=$path1."/".$file5b;
my $infile7a=$path1."/".$file7a;
my $infile7b=$path1."/".$file7b;

my %FamVarsLK=();
my %Genes=();
my %GenoFams=();
my %GenoVarFams=();
my %GenoVarAllele=();
my %GeneInfo=();
my %InHouse5=();
my %PanApp=();
my %PanCounts=();
my %HPOterms=();
my %Solved=();
my %Anno7=();
my @PKs=();

## Read GRCh38 Genes of interest into %PanApp##
open INPUT6, $PAgenefile or die "Cannot open $PAgenefile\n";
#my $head6=<INPUT6>;
loop5: while (<INPUT6>){
	 my $Line6=$_;
	 chomp $Line6;
	 $Line6 =~s/\"//g;
		if(!exists $PanApp{$Line6}){$PanApp{$Line6}=0;}
		if(!exists $PanCounts{$Line6}){$PanCounts{$Line6}=0;}	
	}
 close INPUT6;
## Read GRCh37 Genes of interest into %PanApp##
open INPUT6b, $PAgenefileb or die "Cannot open $PAgenefileb\n";
#my $head6=<INPUT6>;
loop5b: while (<INPUT6b>){
	 my $Line6b=$_;
	 chomp $Line6b;
	 $Line6b =~s/\"//g;
		if(!exists $PanApp{$Line6b}){$PanApp{$Line6b}=0;}
		if(!exists $PanCounts{$Line6b}){$PanCounts{$Line6b}=0;}	
	}
 close INPUT6b;
 
## Read GRCh38 'case' and 'control' (Non-bronchiectasis) Rare disease relative variant counts and MAFs  into %InHouse##
open INPUT5, $infile5a or die "Cannot open $infile5a\n";
my $head5=<INPUT5>;
loop5: while (<INPUT5>){
	 my $Line5=$_;
	 chomp $Line5;
	 $Line5=~s/\"//g;
	 my @Linesplit5=split(/\t/,$Line5);
	 my $var=$Linesplit5[0]."_".$Linesplit5[1]."_".$Linesplit5[2]."_".$Linesplit5[3];
		if(!exists $InHouse5{$var}){
			my $maf=0;
			#NOT Bronch Relatives
			$InHouse5{$var}{'NBrR'}{'Hets'}=$Linesplit5[6];
			$InHouse5{$var}{'NBrR'}{'Homs'}=$Linesplit5[7];
			$maf = ($Linesplit5[6] + (2 * $Linesplit5[7])) / (2 * 33376); ## need to change number of 'controls' depending on cohort
			my $roundedN = sprintf("%.6f", $maf);
			$InHouse5{$var}{'NBrR'}{'MAF'}=$roundedN;
			#NonCF Bronch probands
			$InHouse5{$var}{'BrP'}{'Hets'}=$Linesplit5[9];
			$InHouse5{$var}{'BrP'}{'Homs'}=$Linesplit5[10];
			$maf = ($InHouse5{$var}{'BrP'}{'Hets'} + (2 * $InHouse5{$var}{'BrP'}{'Homs'})) / (2 * 123);  ## need to change number of probands depending on cohort
			my $roundedC = sprintf("%.6f", $maf);
			$InHouse5{$var}{'BrP'}{'MAF'}=$roundedC;
			#Bronch relatives
			$InHouse5{$var}{'BrR'}{'Hets'}=$Linesplit5[12];
			$InHouse5{$var}{'BrR'}{'Homs'}=$Linesplit5[13];
			$maf = ($Linesplit5[13] + (2 * $Linesplit5[13])) / (2 * 260);  ## need to change number of relatives depending on cohort
			my $roundedI = sprintf("%.6f", $maf);
			$InHouse5{$var}{'BrR'}{'MAF'}=$roundedI;
			}
 }
 close INPUT5;
## Read GRCh37 'case' and 'control' (Non-bronchiectasis) Rare disease relative variant counts and MAFs  into %InHouse##
open INPUT5b, $infile5b or die "Cannot open $infile5b\n";
my $head5b=<INPUT5b>;
loop5: while (<INPUT5b>){
	 my $Line5=$_;
	 chomp $Line5;
	 $Line5=~s/\"//g;
	 my @Linesplit5=split(/\t/,$Line5);
	 my $var=$Linesplit5[0]."_".$Linesplit5[1]."_".$Linesplit5[2]."_".$Linesplit5[3];
		if(!exists $InHouse5{$var}){
			my $maf=0;
			#NOT Bronch Relatives
			$InHouse5{$var}{'NBrR'}{'Hets'}=$Linesplit5[6];
			$InHouse5{$var}{'NBrR'}{'Homs'}=$Linesplit5[7];
			$maf = ($Linesplit5[6] + (2 * $Linesplit5[7])) / (2 * 5696);  ## need to change number of 'controls' depending on cohort
			my $roundedN = sprintf("%.6f", $maf);
			$InHouse5{$var}{'NBrR'}{'MAF'}=$roundedN;
			#All Bronch probands
			$InHouse5{$var}{'BrP'}{'Hets'}=$Linesplit5[9];
			$InHouse5{$var}{'BrP'}{'Homs'}=$Linesplit5[10];
			$maf = ($Linesplit5[9] + (2 * $Linesplit5[10])) / (2 * 53);  ## need to change number of probands depending on cohort
			my $roundedC = sprintf("%.6f", $maf);
			$InHouse5{$var}{'BrP'}{'MAF'}=$roundedC;
			#All Bronch relatives
			$InHouse5{$var}{'BrR'}{'Hets'}=$Linesplit5[12];
			$InHouse5{$var}{'BrR'}{'Homs'}=$Linesplit5[13];
			$maf = ($Linesplit5[12] + (2 * $Linesplit5[13])) / (2 * 59);  ## need to change number of relatives depending on cohort
			my $roundedI = sprintf("%.6f", $maf);
			$InHouse5{$var}{'BrR'}{'MAF'}=$roundedI;
			}
 }
 close INPUT5b;

## Read in GRCh38 variant annotation information##
open INPUT7, $infile7a or die "Cannot open $infile7a\n";
my $head7=<INPUT7>;
loop7: while (<INPUT7>){
	 my $Line7=$_;
	 chomp $Line7;
	 $Line7=~s/\"//g;
	 my @Linesplit7=split(/\t/,$Line7);
	 my $var7=$Linesplit7[0]."_".$Linesplit7[1]."_".$Linesplit7[2]."_".$Linesplit7[3];
		
		## Filter out common variants
		unless(!defined $InHouse5{$var7}{'NBrR'}{'MAF'}){unless($InHouse5{$var7}{'NBrR'}{'MAF'}!~/\d/){if($InHouse5{$var7}{'NBrR'}{'MAF'}>=0.001){next loop7;}}}#unless !defined InHouse MAF (Non-Immune relatives)
		unless(!defined $InHouse5{$var7}{'NBrR'}{'Homs'}){unless($InHouse5{$var7}{'NBrR'}{'Homs'}!~/\d/){if($InHouse5{$var7}{'NBrR'}{'Homs'}>0){next loop7;}}} ##Don't include variants that are HOMOzygous in 'controls'
		#unless(!defined $InHouse5{$var7}{'BrP'}{'MAF'}){unless($InHouse5{$var7}{'BrP'}{'MAF'}!~/\d/){if($InHouse5{$var7}{'BrP'}{'MAF'}>=0.1){next loop7;}}}#unless !defined InHouse MAF (Immune Probands)

		#unless($Linesplit7[14]!~/\d/){if($Linesplit7[14]>0.1){next loop7;}}#maxAF Gnomad

		if(!exists $Anno7{$var7}){
			$Anno7{$var7}{'symbol'}=$Linesplit7[4];
			$Anno7{$var7}{'consq'}=$Linesplit7[8];
			$Anno7{$var7}{'hgvsc'}=$Linesplit7[9];
			$Anno7{$var7}{'hgvsp'}=$Linesplit7[10];
				$Linesplit7[11]=~s/\&/, /;
				$Linesplit7[11]=~s/\_/ /;
			$Anno7{$var7}{'rsid'}=$Linesplit7[11];
			$Anno7{$var7}{'sift'}=$Linesplit7[12];
			$Anno7{$var7}{'polyphen'}=$Linesplit7[13];
			$Anno7{$var7}{'maxaf'}=$Linesplit7[14];
			$Anno7{$var7}{'maxpop'}=$Linesplit7[15];
				if($Anno7{$var7}{'maxpop'}=~/\S+\&\S+\&/){$Anno7{$var7}{'maxpop'}="gnomAD_ALL";}
				$Linesplit7[16]=~s/\&/, /;
				$Linesplit7[16]=~s/\_/ /;
			$Anno7{$var7}{'clinsig'}=$Linesplit7[16];
			$Anno7{$var7}{'UTR5cq'}=$Linesplit7[17];
			}
 }
 close INPUT7;

## Read in GRCh37 variant annotation information##
open INPUT7b, $infile7b or die "Cannot open $infile7b\n";
my $head7b=<INPUT7b>;
loop7: while (<INPUT7b>){
	 my $Line7=$_;
	 chomp $Line7;
	 $Line7=~s/\"//g;
	 my @Linesplit7=split(/\t/,$Line7);
	 my $var7=$Linesplit7[0]."_".$Linesplit7[1]."_".$Linesplit7[2]."_".$Linesplit7[3];
		
		## Filter out common variants
		unless(!defined $InHouse5{$var7}{'NBrR'}{'MAF'}){unless($InHouse5{$var7}{'NBrR'}{'MAF'}!~/\d/){if($InHouse5{$var7}{'NBrR'}{'MAF'}>=0.001){next loop7;}}}#unless !defined InHouse MAF (Non-Immune relatives)
		unless(!defined $InHouse5{$var7}{'NBrR'}{'Homs'}){unless($InHouse5{$var7}{'NBrR'}{'Homs'}!~/\d/){if($InHouse5{$var7}{'NBrR'}{'Homs'}>0){next loop7;}}} ##Don't include variants that are HOMOzygous in 'controls'
		#unless(!defined $InHouse5{$var7}{'BrP'}{'MAF'}){unless($InHouse5{$var7}{'BrP'}{'MAF'}!~/\d/){if($InHouse5{$var7}{'BrP'}{'MAF'}>=0.1){next loop7;}}}#unless !defined InHouse MAF (Immune Probands)

		#unless($Linesplit7[14]!~/\d/){if($Linesplit7[14]>0.1){next loop7;}}#maxAF Gnomad

		if(!exists $Anno7{$var7}){
			$Anno7{$var7}{'symbol'}=$Linesplit7[4];
			$Anno7{$var7}{'consq'}=$Linesplit7[8];
			$Anno7{$var7}{'hgvsc'}=$Linesplit7[9];
			$Anno7{$var7}{'hgvsp'}=$Linesplit7[10];
				$Linesplit7[11]=~s/\&/, /;
				$Linesplit7[11]=~s/\_/ /;
			$Anno7{$var7}{'rsid'}=$Linesplit7[11];
			$Anno7{$var7}{'sift'}=$Linesplit7[12];
			$Anno7{$var7}{'polyphen'}=$Linesplit7[13];
			$Anno7{$var7}{'maxaf'}=$Linesplit7[14];
			$Anno7{$var7}{'maxpop'}=$Linesplit7[15];
				if($Anno7{$var7}{'maxpop'}=~/\S+\&\S+\&/){$Anno7{$var7}{'maxpop'}="gnomAD_ALL";}
				$Linesplit7[16]=~s/\&/, /;
				$Linesplit7[16]=~s/\_/ /;
			$Anno7{$var7}{'clinsig'}=$Linesplit7[16];
			$Anno7{$var7}{'UTR5cq'}=$Linesplit7[17];
			}
 }
 close INPUT7b;

## Read txt information file and output per GENE html files ##
open INPUT1, $infile or die "Cannot open $infile\n";
my $head1=<INPUT1>;
loop1: while (<INPUT1>){
	 my $Line1=$_;
	 chomp $Line1;
	 $Line1=~s/\"//g;
	 my @Linesplit=split(/\t/,$Line1);
		## %FamVarsLK{platekeyID}{'family'}=FamilyID
		if(!exists $FamVarsLK{$Linesplit[6]}){
			$FamVarsLK{$Linesplit[6]}{'Family'}=$Linesplit[0];
			$FamVarsLK{$Linesplit[6]}{'RelativeID'}=$Linesplit[1];
			$FamVarsLK{$Linesplit[6]}{'RelativeType'}=$Linesplit[7];
			push(@PKs,$Linesplit[6]);
			}
		if(!exists $HPOterms{$Linesplit[0]}){
			$HPOterms{$Linesplit[0]}=$Linesplit[4];
			}
		if($Linesplit[7]=~/Proband/i){	
			if(!exists $Solved{$Linesplit[0]}){
				$Solved{$Linesplit[0]}=$Linesplit[20];
				}
			}#if proband	
}
close INPUT1;

##Read per Platekey ID Genes/variants into %Genes
open INPUT2, $infile2 or die "Cannot open $infile2\n";
my $head2=<INPUT2>;
loop2: while (<INPUT2>){
	 my $Line2=$_;
	 chomp $Line2;
	 $Line2=~s/\"//g;
	 my @Linesplit2=split(/\t/,$Line2);
		##ONLY INCLUDE NON-REF GENOTYPES
		my @varsplit=split(/\_/,$Linesplit2[1]);
		if($varsplit[2] =~ /^$varsplit[3]$/){next loop2;} #exclude homozygous reference genotypes
	
		##Filter common variants
		unless(!defined $InHouse5{$Linesplit2[1]}{'NBrR'}{'MAF'}){unless($InHouse5{$Linesplit2[1]}{'NBrR'}{'MAF'}!~/\d/){if($InHouse5{$Linesplit2[1]}{'NBrR'}{'MAF'}>=0.001){next loop2;}}}#unless !defined InHouse MAF (Non-Immune relatives)
		unless(!defined $InHouse5{$Linesplit2[1]}{'NBrR'}{'Homs'}){unless($InHouse5{$Linesplit2[1]}{'NBrR'}{'Homs'}!~/\d/){if($InHouse5{$Linesplit2[1]}{'NBrR'}{'Homs'}>0){next loop2;}}} ##Don't include variants that are HOMOzygous in 'controls'
		#unless(!defined $InHouse5{$Linesplit2[1]}{'BrP'}{'MAF'}){unless($InHouse5{$Linesplit2[1]}{'BrP'}{'MAF'}!~/\d/){if($InHouse5{$Linesplit2[1]}{'BrP'}{'MAF'}>=0.1){next loop2;}}}#unless !defined InHouse MAF (Immune Probands)
	
		#unless(!defined $Anno7{$Linesplit2[1]}{'maxaf'}){unless($Anno7{$Linesplit2[1]}{'maxaf'}!~/\d/){if($Anno7{$Linesplit2[1]}{'maxaf'}>=0.1){next loop2;}}}#maxAF Gnomad
	
		##if ONLY INCLUDE PANEL GENES option is selected!!
		if($panelonly=~/yes/i){
			my @SYMsplit=();
			my $pgen="no";
			if($Linesplit2[3]=~/\, /){
				@SYMsplit=split(/\, /,$Linesplit2[3]);
			}elsif($Linesplit2[3]!~/\, /){$SYMsplit[0]=$Linesplit2[3];}
		for(my $c=0;$c<scalar(@SYMsplit);$c++){
			if(exists $PanCounts{$SYMsplit[$c]}){$pgen="yes";}
			}
		if($pgen=~/no/){next loop2;}## ONLY include PANEL genes
			}##panelonly if
		
		##split multi-comma separated gene names
		my @csgenes=();
		if($Linesplit2[3]=~/\, /){
			@csgenes=split(/\, /,$Linesplit2[3]);
			}elsif($Linesplit2[3]!~/\, /){$csgenes[0]=$Linesplit2[3];}
	
	#initialise %Genes{Gene}{Variant}{maf/popn/clinsig}
	foreach my $csg (@csgenes){
			if ($csg =~/\//){$csg=~s/\//\-/g;} ##substitute gene names that include a \/ !!
			if(!exists $Genes{$csg}{$Linesplit2[1]}){
			 $Genes{$csg}{$Linesplit2[1]}{'maf'}=$Linesplit2[4];
			 if($Linesplit2[5]=~/gnomAD_AFR&gnomAD_AMR&gnomAD_ASJ&gnomAD_EAS&gnomAD_FIN&gnomAD_NFE&gnomAD_OTH&gnomAD_SAS/){$Linesplit2[5]="gnomAD_ALL";}
			 $Genes{$csg}{$Linesplit2[1]}{'popn'}=$Linesplit2[5];
			 $Genes{$csg}{$Linesplit2[1]}{'clinsig'}=$Linesplit2[6];
			 $GeneInfo{$csg}{'descr'}=".";
			 $GeneInfo{$csg}{'pheno'}=".";
			 $GeneInfo{$csg}{'ensg'}=".";
			}
		
		#Add family IDs to to %GenoVarFams{gene}{variant}{familyID}
		if(!exists $FamVarsLK{$Linesplit2[0]}){
			# print "Family info not available for PKid: $Linesplit2[0]\n"; 
			next loop2;
			}
		my $family=$FamVarsLK{$Linesplit2[0]}{'Family'};
		my $relType=$FamVarsLK{$Linesplit2[0]}{'RelativeType'};
		if(!exists $GenoVarFams{$csg}{$Linesplit2[1]}{$family}){
			$GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'proband'}=0;
			$GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'father'}=0;
			$GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'mother'}=0;
			$GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'Other'}=0;
			}
		#Check REF and ALT alleles from genotype, against those of variant call !!! indels !!! $GenoVarAllele{gene symbol}{variant}
		if(!exists $GenoVarAllele{$csg}{$Linesplit2[1]}){
		$GenoVarAllele{$csg}{$Linesplit2[1]}{'ref'}="$varsplit[2]";
		$GenoVarAllele{$csg}{$Linesplit2[1]}{'alt'}="$varsplit[3]";
		}
		#Add family IDs to to %GenoFams{gene}{familyID}
		if(!exists $GenoFams{$csg}{$family}){
		$GenoFams{$csg}{$family}{'proband'}{'het'}=".";
		$GenoFams{$csg}{$family}{'proband'}{'hem'}=".";
		$GenoFams{$csg}{$family}{'proband'}{'hom'}=".";
		$GenoFams{$csg}{$family}{'father'}{'het'}=".";
		$GenoFams{$csg}{$family}{'father'}{'hom'}=".";
		$GenoFams{$csg}{$family}{'father'}{'hem'}=".";
		$GenoFams{$csg}{$family}{'mother'}{'het'}=".";
		$GenoFams{$csg}{$family}{'mother'}{'hom'}=".";
		$GenoFams{$csg}{$family}{'mother'}{'hem'}=".";
		$GenoFams{$csg}{$family}{'Other'}{'het'}=".";
		$GenoFams{$csg}{$family}{'Other'}{'hom'}=".";
		$GenoFams{$csg}{$family}{'Other'}{'hem'}=".";
		}
		##Add per gene genotype counts for proband
		if($relType=~/Proband/i and $GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'proband'}<1){
			unless($Linesplit2[2]=~/^NA$/){
				if($GenoFams{$csg}{$family}{'proband'}{'het'}=~/\./){$GenoFams{$csg}{$family}{'proband'}{'het'}=0;$GenoFams{$csg}{$family}{'proband'}{'hom'}=0;$GenoFams{$csg}{$family}{'proband'}{'hem'}=0;}
				my @gtsplit=split(/\//,$Linesplit2[2]);
				unless($Linesplit2[2]!~/\// or $Linesplit2[1]=~/Mb/){if($gtsplit[0]!~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'proband'}{'het'}++;$GenoVarAllele{$csg}{$Linesplit2[1]}{'ref'}=$gtsplit[0];$GenoVarAllele{$csg}{$Linesplit2[1]}{'alt'}=$gtsplit[1];}}
				unless($Linesplit2[2]!~/\// or $Linesplit2[1]=~/Mb/){if($gtsplit[0]=~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'proband'}{'hom'}++;$GenoVarAllele{$csg}{$Linesplit2[1]}{'alt'}=$gtsplit[1];}}
				if($Linesplit2[2]=~/^[ACTG\-]$/){$GenoFams{$csg}{$family}{'proband'}{'hom'}++;}#ChrX male single allele #make sure CN2 LOH vars are not added here!!
				if($Linesplit2[2]=~/^CN[0456789]$/){$GenoFams{$csg}{$family}{'proband'}{'hom'}++;}
				if($Linesplit2[2]=~/^CN[13]$/){$GenoFams{$csg}{$family}{'proband'}{'het'}++;}
				unless($Linesplit2[2]!~/\//){if($Linesplit2[1]=~/Mb/ and $Linesplit2[2]!~/^CN/ and $gtsplit[0]!~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'proband'}{'het'}++;}}
				unless($Linesplit2[2]!~/\//){if($Linesplit2[1]=~/Mb/ and $Linesplit2[2]!~/^CN/ and $gtsplit[0]=~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'proband'}{'hom'}++;}}
				$GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'proband'}++;
				}#unless proband has NA genotype
			}#if proband and the particular variant has not already been counted/printed
		##Add per gene genotype counts for father
		if($relType=~/Father/i and $GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'father'}<1){
			unless($Linesplit2[2]=~/^NA$/){
				if($GenoFams{$csg}{$family}{'father'}{'het'}=~/\./){$GenoFams{$csg}{$family}{'father'}{'het'}=0;$GenoFams{$csg}{$family}{'father'}{'hom'}=0;$GenoFams{$csg}{$family}{'father'}{'hem'}=0;}
				my @gtsplit=split(/\//,$Linesplit2[2]);
				unless($Linesplit2[2]!~/\// or $Linesplit2[1]=~/Mb/){if($gtsplit[0]!~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'father'}{'het'}++;$GenoVarAllele{$csg}{$Linesplit2[1]}{'ref'}=$gtsplit[0];$GenoVarAllele{$csg}{$Linesplit2[1]}{'alt'}=$gtsplit[1];}}
				unless($Linesplit2[2]!~/\// or $Linesplit2[1]=~/Mb/){if($gtsplit[0]=~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'father'}{'hom'}++;$GenoVarAllele{$csg}{$Linesplit2[1]}{'alt'}=$gtsplit[1];}}
				if($Linesplit2[2]=~/^[ACTG\-]$/){$GenoFams{$csg}{$family}{'father'}{'hom'}++;}#ChrX male single allele #make sure CN2 LOH vars are not added here!!
				if($Linesplit2[2]=~/^CN[0456789]$/){$GenoFams{$csg}{$family}{'father'}{'hom'}++;}
				if($Linesplit2[2]=~/^CN[13]$/){$GenoFams{$csg}{$family}{'father'}{'het'}++;}
				unless($Linesplit2[2]!~/\//){if($Linesplit2[1]=~/Mb/ and $Linesplit2[2]!~/^CN/ and $gtsplit[0]!~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'father'}{'het'}++;}}
				unless($Linesplit2[2]!~/\//){if($Linesplit2[1]=~/Mb/ and $Linesplit2[2]!~/^CN/ and $gtsplit[0]=~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'father'}{'hom'}++;}}
				$GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'father'}++;
				}#unless NA genotype
			}#if participant and the particular variant has not already been counted/printed
		##Add per gene genotype counts for mother
		if($relType=~/Mother/i and $GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'mother'}<1){
			unless($Linesplit2[2]=~/^NA$/){
				my @gtsplit=split(/\//,$Linesplit2[2]);
				if($GenoFams{$csg}{$family}{'mother'}{'het'}=~/\./){$GenoFams{$csg}{$family}{'mother'}{'het'}=0;$GenoFams{$csg}{$family}{'mother'}{'hom'}=0;$GenoFams{$csg}{$family}{'mother'}{'hem'}=0;}				
				unless($Linesplit2[2]!~/\// or $Linesplit2[1]=~/Mb/){if($gtsplit[0]!~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'mother'}{'het'}++;$GenoVarAllele{$csg}{$Linesplit2[1]}{'ref'}=$gtsplit[0];$GenoVarAllele{$csg}{$Linesplit2[1]}{'alt'}=$gtsplit[1];}}
				unless($Linesplit2[2]!~/\// or $Linesplit2[1]=~/Mb/){if($gtsplit[0]=~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'mother'}{'hom'}++;$GenoVarAllele{$csg}{$Linesplit2[1]}{'alt'}=$gtsplit[1];}}
				if($Linesplit2[2]=~/^[ACTG\-]$/){$GenoFams{$csg}{$family}{'mother'}{'hom'}++;}#ChrX male single allele #make sure CN2 LOH vars are not added here!!
				if($Linesplit2[2]=~/^CN[0456789]$/){$GenoFams{$csg}{$family}{'mother'}{'hom'}++;}
				if($Linesplit2[2]=~/^CN[13]$/){$GenoFams{$csg}{$family}{'mother'}{'het'}++;}
				unless($Linesplit2[2]!~/\//){if($Linesplit2[1]=~/Mb/ and $Linesplit2[2]!~/^CN/ and $gtsplit[0]!~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'mother'}{'het'}++;}}
				unless($Linesplit2[2]!~/\//){if($Linesplit2[1]=~/Mb/ and $Linesplit2[2]!~/^CN/ and $gtsplit[0]=~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'mother'}{'hom'}++;}}
				$GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'mother'}++;
				}#unless NA genotype
			}#if participant and the particular variant has not already been counted/printed
		##Add per gene genotype counts for ANY Other Relative
		if($relType!~/Proband/i and $relType!~/Father/i and $relType!~/Mother/i and $GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'Other'}<1){
			unless($Linesplit2[2]=~/^NA$/){
				my @gtsplit=split(/\//,$Linesplit2[2]);
				if($GenoFams{$csg}{$family}{'Other'}{'het'}=~/\./){$GenoFams{$csg}{$family}{'Other'}{'het'}=0;$GenoFams{$csg}{$family}{'Other'}{'hom'}=0;$GenoFams{$csg}{$family}{'Other'}{'hem'}=0;}				
				unless($Linesplit2[2]!~/\// or $Linesplit2[1]=~/Mb/){if($gtsplit[0]!~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'Other'}{'het'}++;$GenoVarAllele{$csg}{$Linesplit2[1]}{'ref'}=$gtsplit[0];$GenoVarAllele{$csg}{$Linesplit2[1]}{'alt'}=$gtsplit[1];}}
				unless($Linesplit2[2]!~/\// or $Linesplit2[1]=~/Mb/){if($gtsplit[0]=~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'Other'}{'hom'}++;$GenoVarAllele{$csg}{$Linesplit2[1]}{'alt'}=$gtsplit[1];}}
				if($Linesplit2[2]=~/^[ACTG\-]$/){$GenoFams{$csg}{$family}{'Other'}{'hom'}++;}#ChrX male single allele #make sure CN2 LOH vars are not added here!!
				if($Linesplit2[2]=~/^CN[0456789]$/){$GenoFams{$csg}{$family}{'Other'}{'hom'}++;}
				if($Linesplit2[2]=~/^CN[13]$/){$GenoFams{$csg}{$family}{'Other'}{'het'}++;}
				unless($Linesplit2[2]!~/\//){if($Linesplit2[1]=~/Mb/ and $Linesplit2[2]!~/^CN/ and $gtsplit[0]!~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'other'}{'het'}++;}}
				unless($Linesplit2[2]!~/\//){if($Linesplit2[1]=~/Mb/ and $Linesplit2[2]!~/^CN/ and $gtsplit[0]=~/^$gtsplit[1]$/){$GenoFams{$csg}{$family}{'other'}{'hom'}++;}}
				$GenoVarFams{$csg}{$Linesplit2[1]}{$family}{'Other'}++;
				}#unless NA genotype
			}#if participant and the particular variant has not already been counted/printed	
		}#foreach individual GENE
	}#while input2 loop2
close INPUT2;

##read info for GRCh37 genes symbols into %GeneInfo (some older gene symbols only present in GRCh37 - read these in first and then overwrite any that exist in GRCh38 with newer info!!)
open INPUT3b, $omimfileb or die "Cannot open $omimfileb\n";
my $head3b=<INPUT3b>;
loop3: while (<INPUT3b>){
	my $Line3=$_;
	chomp $Line3;
	$Line3=~s/\"//g;
	my @ls=split(/\t/,$Line3);
	if(exists $GeneInfo{$ls[3]}){
		$GeneInfo{$ls[3]}{'ensg'}=$ls[4];
		$GeneInfo{$ls[3]}{'descr'}=$ls[5];
		$GeneInfo{$ls[3]}{'pheno'}=$ls[6];
		}
	}
close INPUT3b;
##read info for GRCh38 gene symbols into %GeneInfo
open INPUT3, $omimfile or die "Cannot open $omimfile\n";
my $head3=<INPUT3>;
loop3: while (<INPUT3>){
	my $Line3=$_;
	chomp $Line3;
	$Line3=~s/\"//g;
	my @ls=split(/\t/,$Line3);
	if(exists $GeneInfo{$ls[3]}){
		$GeneInfo{$ls[3]}{'ensg'}=$ls[4];
		$GeneInfo{$ls[3]}{'descr'}=$ls[5];
		$GeneInfo{$ls[3]}{'pheno'}=$ls[6];
		}
	}
close INPUT3;
##############################################################

## Print Summary file of variant information ##
open(OUT2, ">$outfile2") || die "Cannot open file \"$outfile2\" to write to!\n";
print OUT2 "FAMILY\tHPO-terms\tSOLVED-status\tGENE\tENSGID\tPHENO\tCHROM\tPOS\tREF\tALT\tBuild\tGeneLocation\tUTR5uORF-CONSQ\t";
print OUT2 "MaxMAF\tMaxMAFpopn\tInhouseMAF-NBrR\tHomozygotes-NBrR\tInhouseMAF-BrP\tHomozygotes-BrP\tInhouseMAF-BrR\tHomozygotes-BrR\tClinvar\n";

open(OUT3, ">$outfile3") || die "Cannot open file \"$outfile3\" to write to!\n";
print OUT3 "FAMILY\tHPO-terms\tSOLVED-STATUS\tGENE\tENSGID\tPHENO\tProband-hets\tProband-homs\tFather-Hets\tFather-Homs\tMother-hets\tMother-Homs\tOtherRelative-hets\tOtherRelative-homs\n";
###########FOR EACH GENE ########################################
foreach my $g (keys %Genes){
##print out each variant per gene
	loopg:foreach my $v (sort keys %{$Genes{$g}}){
		
		my $build="GRCh37";
		if($v=~/chr/){$build="GRCh38";}
		
		my @nvar=split(/\_/,$v);
	
		$Anno7{$v}{'consq'}=~s/\,/\, /g;
		
		unless(defined $InHouse5{$v}{'NBrR'}{'MAF'}){$InHouse5{$v}{'NBrR'}{'MAF'}="\.";}
		unless(defined $InHouse5{$v}{'NBrR'}{'Homs'}){$InHouse5{$v}{'NBrR'}{'Homs'}="\.";}
		unless(defined $InHouse5{$v}{'BrP'}{'MAF'}){$InHouse5{$v}{'BrP'}{'MAF'}="\.";}
		unless(defined $InHouse5{$v}{'BrR'}{'MAF'}){$InHouse5{$v}{'BrR'}{'MAF'}="\.";}
		
		foreach my $fm (sort keys %{$GenoVarFams{$g}{$v}}){
			if($GenoVarFams{$g}{$v}{$fm}{'proband'}>0){
			
				print OUT2 "$fm\t$HPOterms{$fm}\t$Solved{$fm}\t";
				print OUT2 "$g\t";
				print OUT2 "$GeneInfo{$g}{'ensg'}\t";
				print OUT2 "$GeneInfo{$g}{'pheno'}\t";
				print OUT2 "$nvar[0]\t$nvar[1]\t$nvar[2]\t$nvar[3]\t";	
				print OUT2 "$build\t";
				print OUT2 "$Anno7{$v}{'consq'}\t";
				print OUT2 "$Anno7{$v}{'UTR5cq'}\t";
				print OUT2 "$Anno7{$v}{'maxaf'}\t$Anno7{$v}{'maxpop'}\t";
				print OUT2 "$InHouse5{$v}{'NBrR'}{'MAF'}\t";
				print OUT2 "$InHouse5{$v}{'NBrR'}{'Homs'}\t";
				print OUT2 "$InHouse5{$v}{'BrP'}{'MAF'}\t";
				print OUT2 "$InHouse5{$v}{'BrP'}{'Homs'}\t";
				print OUT2 "$InHouse5{$v}{'BrR'}{'MAF'}\t";
				print OUT2 "$InHouse5{$v}{'BrR'}{'Homs'}\t";
				print OUT2 "$Anno7{$v}{'clinsig'}\n";
			}#if proband has a genotype
		}#foreach family
	}#foreach variant

##Per family potential HET and HOM check list
	foreach my $f (sort keys %{$GenoFams{$g}}){
		if($GenoFams{$g}{$f}{'proband'}{'het'}>=1 or $GenoFams{$g}{$f}{'proband'}{'hom'}>=1 or $GenoFams{$g}{$f}{'proband'}{'hem'}>=1){
		print OUT3 "$f\t$HPOterms{$f}\t$Solved{$f}\t$g\t";
		print OUT3 "$GeneInfo{$g}{'ensg'}\t";
		print OUT3 "$GeneInfo{$g}{'pheno'}\t";
		print OUT3 "$GenoFams{$g}{$f}{'proband'}{'het'}\t$GenoFams{$g}{$f}{'proband'}{'hom'}\t";
		print OUT3 "$GenoFams{$g}{$f}{'father'}{'het'}\t$GenoFams{$g}{$f}{'father'}{'hom'}\t";
		print OUT3 "$GenoFams{$g}{$f}{'mother'}{'het'}\t$GenoFams{$g}{$f}{'mother'}{'hom'}\t";
		print OUT3 "$GenoFams{$g}{$f}{'Other'}{'het'}\t$GenoFams{$g}{$f}{'Other'}{'hom'}\n";
		}#potential comp het or hom
	}#foreach family

}#foreach gene
close OUT2;
close OUT3;

exit;
