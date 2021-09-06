#!/usr/bin/perl

#Script to convert vcf file to tab delim text file
#UTR script

use strict;
use warnings;

my($path1,$file,$symbols,$number)=@ARGV;

#define variables
my $vcf_file = $path1."/".$file;
my $outfile = $path1."/".$file."_Variants.txt";
my $maxMAF=-1.0;
my %pidgenes=();

##open gene symbols file and read into %pidgenes
open INPUT3, $symbols or die "Cannot open $symbols\n";
        foreach my $gene (<INPUT3>){
                chomp $gene;
		$gene=~s/\"//g;
	if(!exists $pidgenes{$gene}){$pidgenes{$gene}=0;}
	}
close INPUT3;

## Read vcf file into %info and print to txt files
open(OUT, ">$outfile") || die "Cannot open file \"$outfile\" to write to!\n";

if($vcf_file =~/.gz$/){
	open(INPUT2, "gunzip -c $vcf_file |") or die "Cannot open compressed file $vcf_file\n";
		}
	else{open INPUT2, $vcf_file or die "Cannot open $vcf_file\n";}

	loop2: while (<INPUT2>){
		my $Line=$_;
                chomp $Line;
                if($Line=~/\#\#/){next loop2;}
		if($Line=~/\?/){$Line=~s/\?/\./g;}
                my @lsplit=split(/\t/,$Line);

                #extract header line info and print relevant column headers to out file
                if($Line=~/\#CHROM/){
                        #SAMBCF ##CHROM POS     ID      REF     ALT     QUAL    FILTER  INFO FORMAT(GT:GQ:GQX:DP:DPF:AD:PL) etc.
                        for(my $c=0;$c<7;$c++){print OUT "$lsplit[$c]\t";}
			my @headA=
			("TOTC","HETS","HOMS",
			"Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE",
			"EXON","INTRON","HGVSc","HGVSp","Existing_variation",
			"STRAND","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID","CANONICAL","TSL",
			"GENE_PHENO","SIFT","PolyPhen","DOMAINS",
			"MAX_AF","MAX_AF_POPS","CLIN_SIG","PHENO",
			"UTR_exist_IFoORF","UTR_exist_OFoORF","exist_uORF","UTR5_annot","UTR5_consq");

			for(my $c1=0;$c1<scalar(@headA);$c1++){print OUT "$headA[$c1]\t";}
                        print OUT "END\n";
			next loop2;
		}#end header line
	
	#Initialise %info keys
	my %info=(Allele=>".",Consequence=>".",IMPACT=>".",SYMBOL=>".",Gene=>".",
	Feature_type=>".",Feature=>".",BIOTYPE=>".",EXON=>".",INTRON=>".",
	HGVSc=>".",HGVSp=>".",Existing_variation=>".",STRAND=>".",VARIANT_CLASS=>".",SYMBOL_SOURCE=>".",HGNC_ID=>".",
	CANONICAL=>".",TSL=>".",GENE_PHENO=>".",SIFT=>".",PolyPhen=>".",
	DOMAINS=>".",MAX_AF=>".",MAX_AF_POPS=>".",CLIN_SIG=>".",PHENO=>".",OMIM=>".",
	TOTC=>".",HETS=>".",HOMS=>".",
	UTR_exist_IFoORF=>".",UTR_exist_OFoORF=>".",exist_uORF=>".",UTR5_annot=>".",UTR5_consq=>".");

	if($Line=~/TOTC\=(\d+)\;HETS\=(\d+)\;HOMS\=(\d+)\;/){
				$info{'TOTC'}=$1;
				$info{'HETS'}=$2;
				$info{'HOMS'}=$3;
			my $inhMAF = ($info{'HETS'} + (2 * $info{'HOMS'})) / (2 * $number);
			my $fracall = $info{'TOTC'} / $number;
			##if($inhMAF>0.01 or $fracall>0.1){$info{'FILTEROUT'}="YES";}
				}
	
	#split CSG entries, first by comma delimiter and then by | delim. Combine all unique info entries
	if($Line=~/\;CSQ\=(\S+?)\t/){
		my $CSQ=$1;
		my @CSQcomma=();
		if($CSQ!~/\,/){$CSQcomma[0]=$CSQ;}
		if($CSQ=~/\,/){@CSQcomma=split(/\,/,$CSQ);}
		$maxMAF=-1.0;
		
	geneloop: foreach my $cq (@CSQcomma){
			my @CSQsplit=split(/\|/,$cq);
			my %ENSIDs=();
			
			if(!exists $pidgenes{$CSQsplit[3]}){next geneloop;}	##only include genes from input list
			
		for(my $val=0;$val<82;$val++){if(!defined $CSQsplit[$val]){$CSQsplit[$val]=".";}}

		#extract annotation information
			if($info{'Allele'}!~/$CSQsplit[0]/){$info{'Allele'}=$info{'Allele'}."\,".$CSQsplit[0];}
			if($info{'Consequence'}!~/$CSQsplit[1]/){$info{'Consequence'}=$info{'Consequence'}."\,".$CSQsplit[1];}
			if($info{'IMPACT'}!~/$CSQsplit[2]/){$info{'IMPACT'}=$info{'IMPACT'}."\,".$CSQsplit[2];}
			if($info{'SYMBOL'}!~/$CSQsplit[3]/){$info{'SYMBOL'}=$info{'SYMBOL'}."\,".$CSQsplit[3];}
			if($info{'Gene'}!~/$CSQsplit[4]/){$info{'Gene'}=$info{'Gene'}."\,".$CSQsplit[4];}
				if(!exists $ENSIDs{$CSQsplit[4]}){$ENSIDs{$CSQsplit[4]}=0;}
			
			if($info{'Feature_type'}!~/$CSQsplit[5]/){$info{'Feature_type'}=$info{'Feature_type'}."\,".$CSQsplit[5];}
			if($info{'Feature'}!~/$CSQsplit[6]/){$info{'Feature'}=$info{'Feature'}."\,".$CSQsplit[6];}
			if($info{'BIOTYPE'}!~/$CSQsplit[7]/){$info{'BIOTYPE'}=$info{'BIOTYPE'}."\,".$CSQsplit[7];}
			if($info{'EXON'}!~/$CSQsplit[8]/){$info{'EXON'}=$info{'EXON'}."\,".$CSQsplit[8];}
			if($info{'INTRON'}!~/$CSQsplit[9]/){$info{'INTRON'}=$info{'INTRON'}."\,".$CSQsplit[9];}
			if($info{'HGVSc'}!~/$CSQsplit[10]/){$info{'HGVSc'}=$info{'HGVSc'}."\,".$CSQsplit[10];}
			if($info{'HGVSp'}!~/$CSQsplit[11]/){$info{'HGVSp'}=$info{'HGVSp'}."\,".$CSQsplit[11];}
			if($info{'Existing_variation'}!~/$CSQsplit[17]/){$info{'Existing_variation'}=$info{'Existing_variation'}."\,".$CSQsplit[17];}
			if($info{'STRAND'}!~/$CSQsplit[19]/){$info{'STRAND'}=$info{'STRAND'}."\,".$CSQsplit[19];}
			if($info{'VARIANT_CLASS'}!~/$CSQsplit[21]/){$info{'VARIANT_CLASS'}=$info{'VARIANT_CLASS'}."\,".$CSQsplit[21];}
			if($info{'SYMBOL_SOURCE'}!~/$CSQsplit[22]/){$info{'SYMBOL_SOURCE'}=$info{'SYMBOL_SOURCE'}."\,".$CSQsplit[22];}
			if($info{'HGNC_ID'}!~/$CSQsplit[23]/){$info{'HGNC_ID'}=$info{'HGNC_ID'}."\,".$CSQsplit[23];}
			if($info{'CANONICAL'}!~/$CSQsplit[24]/){$info{'CANONICAL'}=$info{'CANONICAL'}."\,".$CSQsplit[24];}
			if($info{'TSL'}!~/$CSQsplit[26]/){$info{'TSL'}=$info{'TSL'}."\,".$CSQsplit[26];}
			
			if($info{'GENE_PHENO'}!~/$CSQsplit[33]/){$info{'GENE_PHENO'}=$info{'GENE_PHENO'}."\,".$CSQsplit[33];}
			if($info{'SIFT'}!~/$CSQsplit[34]/){$info{'SIFT'}=$info{'SIFT'}."\,".$CSQsplit[34];}
			if($info{'PolyPhen'}!~/$CSQsplit[35]/){$info{'PolyPhen'}=$info{'PolyPhen'}."\,".$CSQsplit[35];}
			if($info{'DOMAINS'}!~/$CSQsplit[36]/){$info{'DOMAINS'}=$info{'DOMAINS'}."\,".$CSQsplit[36];}
				
			if($info{'MAX_AF'}!~/$CSQsplit[56]/){$info{'MAX_AF'}=$info{'MAX_AF'}."\,".$CSQsplit[56];}
				$info{'MAX_AF'}=~s/^\.\,//; $info{'MAX_AF'}=~s/^\,//;
				if($info{'MAX_AF'}=~/\d+/){if($info{'MAX_AF'}>$maxMAF){$maxMAF=$info{'MAX_AF'};}}
			if($info{'MAX_AF_POPS'}!~/$CSQsplit[57]/){$info{'MAX_AF_POPS'}=$info{'MAX_AF_POPS'}."\,".$CSQsplit[57];}
			if($info{'CLIN_SIG'}!~/$CSQsplit[58]/){$info{'CLIN_SIG'}=$info{'CLIN_SIG'}."\,".$CSQsplit[58];}
			if($info{'PHENO'}!~/$CSQsplit[60]/){$info{'PHENO'}=$info{'PHENO'}."\,".$CSQsplit[60];}

			if($info{'UTR_exist_IFoORF'}!~/$CSQsplit[66]/){$info{'UTR_exist_IFoORF'}=$info{'UTR_exist_IFoORF'}."\,".$CSQsplit[66];}
			if($info{'UTR_exist_OFoORF'}!~/$CSQsplit[67]/){$info{'UTR_exist_OFoORF'}=$info{'UTR_exist_OFoORF'}."\,".$CSQsplit[67];}
			if($info{'exist_uORF'}!~/$CSQsplit[68]/){$info{'exist_uORF'}=$info{'exist_uORF'}."\,".$CSQsplit[68];}
			if($info{'UTR5_annot'}!~/$CSQsplit[69]/){$info{'UTR5_annot'}=$info{'UTR5_annot'}."\,".$CSQsplit[69];}
			if($info{'UTR5_consq'}!~/$CSQsplit[70]/){$info{'UTR5_consq'}=$info{'UTR_consq'}."\,".$CSQsplit[70];}

		}#foreach , sep CSQ
	}#if CSQ

	foreach my $k (keys %info){
		$info{$k}=~s/^\.\,//;
		$info{$k}=~s/\,{2,}/\,/g;
		$info{$k}=~s/^\,//;
		$info{$k}=~s/\,$//;
		$info{$k}=~s/^\,$/\./;
		if($info{$k}=~/^\s+$/){$info{$k}=".";}
		if($info{$k}=~/^$/){$info{$k}=".";}
		}

	my @lineA=("$info{'TOTC'}","$info{'HETS'}","$info{'HOMS'}","$info{'Allele'}","$info{'Consequence'}","$info{'IMPACT'}","$info{'SYMBOL'}","$info{'Gene'}",
	"$info{'Feature_type'}","$info{'Feature'}","$info{'BIOTYPE'}","$info{'EXON'}","$info{'INTRON'}","$info{'HGVSc'}","$info{'HGVSp'}",
	"$info{'Existing_variation'}","$info{'STRAND'}","$info{'VARIANT_CLASS'}","$info{'SYMBOL_SOURCE'}","$info{'HGNC_ID'}",
	"$info{'CANONICAL'}","$info{'TSL'}","$info{'GENE_PHENO'}","$info{'SIFT'}","$info{'PolyPhen'}","$info{'DOMAINS'}",
	"$info{'MAX_AF'}","$info{'MAX_AF_POPS'}","$info{'CLIN_SIG'}","$info{'PHENO'}",
	"$info{'UTR_exist_IFoORF'}","$info{'UTR_exist_OFoORF'}","$info{'exist_uORF'}","$info{'UTR5_annot'}","$info{'UTR5_consq'}");
        
	### print line info to output txt file
	unless($info{'HETS'}==0 and $info{'HOMS'}==0){
		for(my $c=0;$c<7;$c++){print OUT "$lsplit[$c]\t";}       
		for(my $c1=0;$c1<scalar(@lineA);$c1++){print OUT "$lineA[$c1]\t";}
		print OUT "END\n";
		}##unless no het or hom genotype counts (Gel included some known pathogenic variants in vcf calling, even if not seen)
	@lineA=();
}# while file lines
close INPUT2;
close OUT;

exit;
