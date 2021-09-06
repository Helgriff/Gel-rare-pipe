##R-code used to combine per variant genotype counts across cohorts

setwd("/file/path/to/analysis/directory")
workDIR <-getwd()

numGa="PCDx92"
TOTparts <- 123+260+33376 ##NCFBr-probands, BronchPlus-relatives & NonBr-relatives

#chromNo <-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")
chromNo <-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","21","X")

for (c in 1:length(chromNo)){

##NCFB
NCFBvarID <- read.table(paste0(workDIR,"/NCFB_probands/perCHROM_norm_",numGa,"/preVariants_NCFB_probands_MAF-noFilt_",numGa,"_chr",chromNo[c],".txt_allCILgeneVars-extended"),comment.char="",header=TRUE,sep="\t",stringsAsFactors = FALSE,
                       colClasses = c("character","numeric","character","character","character","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character"))
NCFBvarIDshort <- NCFBvarID[,1:8]

##BRONCH relatives
BRRvarID <- read.table(paste0(workDIR,"/BRONCHplus_relatives/perCHROM_norm_",numGa,"/preVariants_BRONCHplus_relatives_MAF-noFilt_",numGa,"_chr",chromNo[c],".txt_allCILgeneVars-extended"),comment.char="",header=TRUE,sep="\t",stringsAsFactors = FALSE,
                        colClasses = c("character","numeric","character","character","character","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character"))
BRRvarIDshort <- BRRvarID[,1:8] 

##NON-BRONCH relatives
NBRvarID <- read.table(paste0(workDIR,"/NON-BRONCH_relatives/perCHROM_norm_",numGa,"/preVariants_NON-BRONCH_relatives_MAF-noFilt_",numGa,"_chr",chromNo[c],".txt_allCILgeneVars-extended"),comment.char="",header=TRUE,sep="\t",stringsAsFactors = FALSE,
                        colClasses = c("character","numeric","character","character","character","numeric","numeric","numeric","character","character","character","character","character","character","character","character","character","character"))
NBRvarIDshort <- NBRvarID[,1:8] 

#SHORT
jnt.vars1 <- full_join(NBRvarIDshort, NCFBvarIDshort, by=c("X.CHROM","POS","REF","ALT","SYMBOL"))
jnt.vars <- full_join(jnt.vars1, BRRvarIDshort, by=c("X.CHROM","POS","REF","ALT","SYMBOL"))
#LONG
jnt.vars1L <- full_join(NBRvarID, NCFBvarID, by=c("X.CHROM","POS","REF","ALT","SYMBOL","Consequence","HGVSc","HGVSp","Existing_variation","SIFT","PolyPhen","MAX_AF","MAX_AF_POPS","CLIN_SIG","UTR5_consq"))
jnt.varsL <- full_join(jnt.vars1L, BRRvarID, by=c("X.CHROM","POS","REF","ALT","SYMBOL","Consequence","HGVSc","HGVSp","Existing_variation","SIFT","PolyPhen","MAX_AF","MAX_AF_POPS","CLIN_SIG","UTR5_consq"))

jnt.vars[is.na(jnt.vars)] <- 0
jnt.varsL[is.na(jnt.varsL)] <- 0

colnames(jnt.vars) <- c("CHROM","POS","REF","ALT","ID","nbrTOT","nbrHET","nbrHOM","ncfbpTOT","ncfbpHET","ncfbpHOM","brrTOT","brrHET","brrHOM")
colnames(jnt.varsL) <- c("CHROM","POS","REF","ALT","ID","nbrTOT","nbrHET","nbrHOM","Consequence","HGVSc","HGVSp","Existing_variation","SIFT","PolyPhen","MAX_AF","MAX_AF_POPS","CLIN_SIG","UTR5_consq","ncfbpTOT","ncfbpHET","ncfbpHOM","brrTOT","brrHET","brrHOM")

HET.nci <- jnt.vars$nbrHET + jnt.vars$ncfbpHET + jnt.vars$brrHET
HOM.nci <- jnt.vars$nbrHOM + jnt.vars$ncfbpHOM + jnt.vars$brrHOM
MAF.nci <- ((2*HOM.nci) + HET.nci)/(2*TOTparts)
HET.nciL <- jnt.varsL$nbrHET + jnt.varsL$ncfbpHET + jnt.varsL$brrHET
HOM.nciL <- jnt.varsL$nbrHOM + jnt.varsL$ncfbpHOM + jnt.varsL$brrHOM
MAF.nciL <- ((2*HOM.nciL) + HET.nciL)/(2*TOTparts)
jnt.vars <- cbind(jnt.vars,HET.nci,HOM.nci,MAF.nci)
jnt.varsL <- cbind(jnt.varsL,HET.nciL,HOM.nciL,MAF.nciL)

write.table(jnt.vars, file=paste0(workDIR,"/ALLunfiltered_Counts_perVariant_noMAF-filter_NonBronchRels_NCFB_BRR_combined_",numGa,"_chr",chromNo[c],".txt"), sep="\t",row.names = F)
write.table(jnt.varsL, file=paste0(workDIR,"/ALLunfiltered_Counts_perVariant_noMAF-filter_NonBronchRels_NCFB_BRR_combined_",numGa,"_chr",chromNo[c],"_extended-Info.txt"), sep="\t",row.names = F)
}

#### END ####
