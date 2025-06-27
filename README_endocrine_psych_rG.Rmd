# Endocrine Hormone and Psychiatric Disorder Genetic Correlation Analysis

# Obtaining EUR allele frequencies to fill in frequencies for studies missing this annotation
```{r}
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# importing 1kg allele frequencies for files that have missing AF fields
freq_1kg <- fread("/gpfs/data/huo-lab/BCAC/james.li/xtrait/data/1kg_reference_bfile/EUR/chr_all.afreq") 
freq_1kg <- freq_1kg %>% select(ID,ALT,ALT_FREQS)
bim <- fread("/gpfs/data/huo-lab/BCAC/james.li/xtrait/data/1kg_reference_bfile/EUR/chr_all.bim")
bim <- bim %>% rename(ID=V2) %>% select(ID,V1,V4,V5,V6)
joined_1kg <- inner_join(freq_1kg,bim,by=c("ID"))
chr_all_ALT_FREQS <- joined_1kg %>% mutate(REF=ifelse(ALT==V5,V6,V5)) %>% rename(CHR=V1,BP=V4) %>% select(CHR,BP,ALT,REF,ALT_FREQS)
save(chr_all_ALT_FREQS,file="/gpfs/data/huo-lab/BCAC/james.li/xtrait/data/1kg_reference_bfile/EUR/chr_all_ALT_FREQS.RData")

# load in rsID dictionary
SNP_GRCh37_38_match_update <- readRDS("/gpfs/data/pierce-lab/rsID_dictionary/SNP_GRCh37_38_match_update.rds")
SNP_GRCh37_38_match_update <- SNP_GRCh37_38_match_update %>% rename(CHR=chr,BP=pos37,A1=allele1_37,A2=allele2_37,snpid=rsid)
SNP_GRCh37_38_match_update <- SNP_GRCh37_38_match_update %>% select(CHR,BP,A1,A2,snpid)

# creating a conversion dictionary for clumping prior to the MR analysis
colnames(bim) <- c("ID","CHR","BP","tmp_A1","tmp_A2")
bim <- bim %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) 
bim <- bim %>% mutate(new_id = paste(CHR,BP,A2,A1,sep=":")) %>% select(ID,new_id)
bim <- bim[!(duplicated(bim$new_id) | duplicated(bim$new_id, fromLast = TRUE)), ]

write.table(bim,file="/gpfs/data/huo-lab/BCAC/james.li/xtrait/data/1kg_reference_bfile/EUR/ID_CONVERSION.txt",quote=F,row.names=F,col.names=T,sep="\t")

system("module load plink/2.0; plink2 --bfile /gpfs/data/huo-lab/BCAC/james.li/xtrait/data/1kg_reference_bfile/EUR/chr_all --update-name /gpfs/data/huo-lab/BCAC/james.li/xtrait/data/1kg_reference_bfile/EUR/ID_CONVERSION.txt 2 1 --make-bed --out /gpfs/data/huo-lab/BCAC/james.li/xtrait/data/1kg_reference_bfile/EUR/converted_chr_all")
```


# Processing GWAS summary statistics
```{bash}
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
 
# import hm3 list
# hm3_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/preconfluence/data/hm3_snplist/merge_alleles.txt"))$SNP

# load in rsID dictionary
SNP_GRCh37_38_match_update <- readRDS("/gpfs/data/pierce-lab/rsID_dictionary/SNP_GRCh37_38_match_update.rds")
SNP_GRCh37_38_match_update <- SNP_GRCh37_38_match_update %>% rename(CHR=chr,BP=pos37,tmp_A1=allele1_37,tmp_A2=allele2_37,snpid=rsid)
SNP_GRCh37_38_match_update <- SNP_GRCh37_38_match_update %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(CHR,BP,A1,A2,snpid) 
simplified_info <- SNP_GRCh37_38_match_update %>% select(CHR,BP,snpid) 

# loading in 1kg allele frequencies and parsing them
load("/gpfs/data/huo-lab/BCAC/james.li/xtrait/data/1kg_reference_bfile/EUR/chr_all_ALT_FREQS.RData")
chr_all_ALT_FREQS <- chr_all_ALT_FREQS %>% mutate(A1=ifelse(ALT<REF,ALT,REF)) %>% mutate(A2=ifelse(ALT<REF,REF,ALT)) %>% mutate(MAF=ifelse(ALT_FREQS<0.5,ALT_FREQS,1-ALT_FREQS))
MAF_1kg <- chr_all_ALT_FREQS %>% select(CHR,BP,A1,A2,MAF)

#############################
#############################
# LIST OF STUDIES AND PMIDS #
#############################
#############################
# AN: 31308545
# TS: 30818990
# OCD: 38548983
# SCZ: 35396580
#	ADHD: 36702997
# BIP: 34002096
# MDD: 30718901
# ANX: 31712720

########################################
########################################
########################################
# Parsing summary statistics for traits
########################################
########################################
########################################

########################################
# importing and parsing ANX
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/raw_trait")
trait_list <- list.files()[grepl("panic",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% filter(`#CHROM` %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(CHR=`#CHROM`,BP=POS,beta=BETA,se=SE,P=PVAL,info=IMPINFO) %>% mutate(EAF=(FCAS*NCAS + FCON*NCON)/(NCAS+NCON), N=NCAS+NCON) %>% mutate(MAF=ifelse(EAF<0.5,EAF,1-EAF))
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=A1,tmp_A2=A2) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)
  output_trait_file <- str_split(trait_list[index],pattern="\\.")[[1]][1]
  write.table(processed_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_trait/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

########################################
# importing and parsing OCD
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/raw_trait")
trait_list <- list.files()[grepl("obsessive",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% filter(CHR %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(se=SE,info=INFO) %>% mutate(beta=log(OR),EAF=(FRQ_A_33943*Nca + FRQ_U_33943*Nco)/(Nca+Nco), N=Nca+Nco) %>% mutate(MAF=ifelse(EAF<0.5,EAF,1-EAF))
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=A1,tmp_A2=A2) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)
  output_trait_file <- str_split(trait_list[index],pattern="-")[[1]][1]
  write.table(processed_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_trait/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

########################################
# importing and parsing AN
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/raw_trait")
trait_list <- list.files()[grepl("AN",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% filter(CHROM %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(CHR=CHROM,BP=POS,beta=BETA,se=SE,info=IMPINFO,P=PVAL) %>% mutate(N=NCAS+NCON)
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=ALT,tmp_A2=REF) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- inner_join(tmp_trait,MAF_1kg,by=c("CHR","BP","A1","A2")) 
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)
  output_trait_file <- str_split(trait_list[index],pattern="\\.")[[1]][1]
  write.table(processed_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_trait/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

########################################
# importing and parsing TS
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/raw_trait")
trait_list <- list.files()[grepl("TS",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% filter(CHR %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(se=SE,info=INFO) %>% mutate(beta=log(OR), N=14307)
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=A1,tmp_A2=A2) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- inner_join(tmp_trait,MAF_1kg,by=c("CHR","BP","A1","A2")) 
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)
  output_trait_file <- str_split(trait_list[index],pattern="_")[[1]][1]
  write.table(processed_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_trait/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

########################################
# importing and parsing MDD
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/raw_trait")
trait_list <- list.files()[grepl("depression",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) %>% rename(snpid=MarkerName)
  tmp_trait <- tmp_trait %>% mutate(A1=toupper(A1),A2=toupper(A2))
  tmp_trait <- inner_join(tmp_trait,simplified_info,by=c("snpid"))
  tmp_trait <- tmp_trait %>% filter(CHR %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(se=StdErrLogOR,beta=LogOR) %>% mutate(info=1,N=500199) %>% mutate(MAF=ifelse(Freq<0.5,Freq,1-Freq))
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=A1,tmp_A2=A2) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)
  output_trait_file <- str_split(trait_list[index],pattern="_")[[1]][3]
  write.table(processed_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_trait/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

########################################
# importing and parsing ADHD
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/raw_trait")
trait_list <- list.files()[grepl("ADHD",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% filter(CHR %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(CHR=CHR,BP=BP,se=SE,P=P,info=INFO) %>% mutate(beta=log(OR),EAF=(FRQ_A_38691*Nca + FRQ_U_186843*Nco)/(Nca+Nco), N=Nca+Nco) %>% mutate(MAF=ifelse(EAF<0.5,EAF,1-EAF))
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=A1,tmp_A2=A2) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)
  output_trait_file <- str_split(trait_list[index],pattern="_")[[1]][1]
  write.table(processed_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_trait/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

########################################
# importing and parsing SCZ
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/raw_trait")
trait_list <- list.files()[grepl("SCZ",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% filter(CHROM %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(CHR=CHROM,BP=POS,beta=BETA,se=SE,P=PVAL,info=IMPINFO) %>% mutate(EAF=(FCAS*NCAS + FCON*NCON)/(NCAS+NCON), N=NCAS+NCON) %>% mutate(MAF=ifelse(EAF<0.5,EAF,1-EAF))
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=A1,tmp_A2=A2) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)
  output_trait_file <- str_split(trait_list[index],pattern="_")[[1]][2]
  write.table(processed_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_trait/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

########################################
# importing and parsing BIP
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/raw_trait")
trait_list <- list.files()[grepl("pgc-bip",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% filter(`#CHROM` %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(CHR=`#CHROM`,BP=POS,beta=BETA,se=SE,P=PVAL,info=IMPINFO) %>% mutate(EAF=(FCAS*NCAS + FCON*NCON)/(NCAS+NCON), N=NCAS+NCON) %>% mutate(MAF=ifelse(EAF<0.5,EAF,1-EAF))
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=A1,tmp_A2=A2) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)
  output_trait_file <- str_split(trait_list[index],pattern="\\.")[[1]][1]
  write.table(processed_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_trait/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

##############################################
##############################################
# importing and parsing non-thyroid hormones #
##############################################
##############################################

########################################
# importing and parsing SHBG
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/data/UKBB/")
trait_list <- list.files()[grepl("SHBG",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% rename(CHR=chr,BP=pos)
  tmp_trait <- tmp_trait %>% filter(CHR %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(se=se_meta,beta=beta_meta) %>% mutate(info=1,N=381526,MAF=ifelse(af_meta>0.5,1-af_meta,af_meta),P=10^(-neglog10_pval_meta))
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=alt,tmp_A2=ref) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)

  # remove rows with NA in beta, se, P, or MAF, and then duplicate snpids
  cleaned_trait <- processed_trait[!is.na(beta) & !is.na(se) & !is.na(P) & !is.na(MAF)]
  duplicate_snpids <- cleaned_trait[, .N, by = snpid][N > 1, snpid]
  cleaned_trait <- cleaned_trait[!(snpid %in% duplicate_snpids)]
  output_trait_file <- str_split(trait_list[index],pattern="_")[[1]][1]
  write.table(cleaned_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_thormone/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

########################################
# importing and parsing estradiol
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/data/UKBB/")
trait_list <- list.files()[grepl("estradiol",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% rename(CHR=chr,BP=pos)
  tmp_trait <- tmp_trait %>% filter(CHR %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(se=se_meta,beta=beta_meta) %>% mutate(info=1,N=67623,MAF=ifelse(af_meta>0.5,1-af_meta,af_meta),P=10^(-neglog10_pval_meta))
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=alt,tmp_A2=ref) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)

  # remove rows with NA in beta, se, P, or MAF, and then duplicate snpids
  cleaned_trait <- processed_trait[!is.na(beta) & !is.na(se) & !is.na(P) & !is.na(MAF)]
  duplicate_snpids <- cleaned_trait[, .N, by = snpid][N > 1, snpid]
  cleaned_trait <- cleaned_trait[!(snpid %in% duplicate_snpids)]
  output_trait_file <- str_split(trait_list[index],pattern="_")[[1]][1]
  write.table(cleaned_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_thormone/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}

########################################
# importing and parsing testosterone
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/data/UKBB/")
trait_list <- list.files()[grepl("testosterone",list.files())]
for (index in 1:length(trait_list)) {
  print(paste("Processing summary statistics for:",trait_list[index]))
  tmp_trait <- fread(trait_list[index]) 
  tmp_trait <- tmp_trait %>% rename(CHR=chr,BP=pos)
  tmp_trait <- tmp_trait %>% filter(CHR %in% c(1:22)) 
  tmp_trait <- tmp_trait %>% rename(se=se_meta,beta=beta_meta) %>% mutate(info=1,N=381081,MAF=ifelse(af_meta>0.5,1-af_meta,af_meta),P=10^(-neglog10_pval_meta))
  # standardizing A1 and A2 alleles
  tmp_trait <- tmp_trait %>% rename(tmp_A1=alt,tmp_A2=ref) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2)
  tmp_trait <- tmp_trait %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  tmp_trait <- tmp_trait %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  processed_trait <- left_join(tmp_trait,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_trait <- processed_trait %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)

  # remove rows with NA in beta, se, P, or MAF, and then duplicate snpids
  cleaned_trait <- processed_trait[!is.na(beta) & !is.na(se) & !is.na(P) & !is.na(MAF)]
  duplicate_snpids <- cleaned_trait[, .N, by = snpid][N > 1, snpid]
  cleaned_trait <- cleaned_trait[!(snpid %in% duplicate_snpids)]
  output_trait_file <- str_split(trait_list[index],pattern="_")[[1]][1]
  write.table(cleaned_trait,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_thormone/",output_trait_file),quote=F,sep="\t",row.names=F,col.names=T)
}


########################################
########################################
# importing and parsing thyroid hormones
########################################
########################################
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/raw_thormone")
thormone_list <- list.files()
for (i in 1:length(thormone_list)) {
    print(paste("Processing summary statistics for:",thormone_list[i]))
  tmp_thormone <- fread(thormone_list[i]) 
  tmp_thormone <- tmp_thormone %>% separate(MarkerName,sep="\\:",into=c("CHR","BP","VAR_TYPE"))
  tmp_thormone <- tmp_thormone %>% mutate(A1 = toupper(Allele1),A2=toupper(Allele2),beta=Effect,se=StdErr,P=P.value,info=1,MAF = ifelse(Freq1>0.5,1-Freq1, Freq1)) 
  tmp_thormone <- tmp_thormone %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  tmp_thormone <- tmp_thormone %>% mutate(CHR=as.integer(CHR),BP=as.integer(BP))
  # standardizing A1 and A2 alleles
  tmp_thormone <- tmp_thormone %>% rename(tmp_A1=A1,tmp_A2=A2) %>% mutate(beta=ifelse(tmp_A1<tmp_A2,beta,-beta)) %>% mutate(A1=ifelse(tmp_A1<tmp_A2,tmp_A1,tmp_A2)) %>% mutate(A2=ifelse(tmp_A1<tmp_A2,tmp_A2,tmp_A1)) %>% select(-tmp_A1,-tmp_A2) 
  processed_thormone <- tmp_thormone %>% select(CHR,BP,A1,A2,beta,se,P,info,MAF,N)
  processed_thormone <- left_join(tmp_thormone,SNP_GRCh37_38_match_update,by=c("CHR","BP","A1","A2"))
  processed_thormone <- processed_thormone %>% rename(bp=BP) %>% select(snpid,CHR,bp,A1,A2,beta,se,P,info,MAF,N)
  output_thormone_file <- gsub("_","",gsub("_overall_130421_invvar1.txt-QCfiltered_GC.txt.gz","",gsub("formatted_","",thormone_list[i])))
  write.table(processed_thormone,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_thormone/",output_thormone_file),quote=F,sep="\t",row.names=F,col.names=T)
}
```

# Munging summary statistics for LDSC
```{bash}
# munge sumstats
for o in thormone trait
do
  cd /gpfs/data/huo-lab/BCAC/james.li/xtrait/input/processed_${o}
  for j in *
  do
    sbatch --export=ARGS1=${o},ARGS2=${j} /gpfs/data/huo-lab/BCAC/james.li/xtrait/code/munge_sumstats_compute_h2.sh
  done
done
```

# Computing genetic correlations
```{bash}
# obtaining thormone munge sumstats list 
cd /gpfs/data/huo-lab/BCAC/james.li/xtrait/input/munge_thormone
thormone_list=`ls *.sumstats.gz`
# obtaining gwas_trait munge sumstats list 
cd /gpfs/data/huo-lab/BCAC/james.li/xtrait/input/munge_trait
trait_list=`ls *.sumstats.gz`

# computing genetic correlations
for thormone in $thormone_list
do
  for trait in $trait_list
  do
    sbatch --export=ARGS1=${thormone},ARGS2=${trait} /gpfs/data/huo-lab/BCAC/james.li/xtrait/code/compute_rg.sh
  done
done

# parsing output from LDSC files into a table
cd /gpfs/data/huo-lab/BCAC/james.li/xtrait/output/rG
for thormone in invnormFT3 invnormFT4 invnormTSH invnormTT3 lnFT3FT4 lnTT3FT4 estradiol testosterone SHBG
do
  header_file=`ls ${thormone}*.log | head -1`
  tail -5 ${header_file} | head -1 > ./parsed_output/${thormone}

  for i in `ls ${thormone}*.log`
  do
    tail -4 $i | head -1 >> ./parsed_output/${thormone}
  done
done
```

# Computing rG only amongst hormones 
```{bash}
# Obtaining thormone munge sumstats list
cd /gpfs/data/huo-lab/BCAC/james.li/xtrait/input/munge_thormone
thormone_list=($(ls *.sumstats.gz))
  
# Computing genetic correlations with all possible pairs (including same files)
num_files=${#thormone_list[@]}
    
for ((i=0; i<num_files; i++)); do
for ((j=0; j<num_files; j++)); do  # Start from 0 to include same-file pairs
  sbatch --export=ARGS1=${thormone_list[i]},ARGS2=${thormone_list[j]} /gpfs/data/huo-lab/BCAC/james.li/xtrait/code/compute_rg_only_hormone.sh
done
done

# parsing output from LDSC files into a table
cd /gpfs/data/huo-lab/BCAC/james.li/xtrait/output/rG_only_hormone
for thormone in invnormFT3 invnormFT4 invnormTSH invnormTT3 lnFT3FT4 lnTT3FT4 estradiol SHBG testosterone
do
  header_file=`ls ${thormone}*.log | head -1`
  tail -5 ${header_file} | head -1 > ./parsed_output/${thormone}

  for i in `ls ${thormone}*.log`
  do
    tail -4 $i | head -1 >> ./parsed_output/${thormone}
  done
done
```


# Analying and plotting cross-trait genetic correlation results in R 
```{r}
library(data.table)
library(dplyr)
library(ggplot2)

#######################
# CROSS-TRAIT RESULTS #
#######################
# set working directory
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/output/rG/parsed_output")
thormone_list <- c(
  "SHBG","invnormFT3","invnormTSH","lnFT3FT4","testosterone","estradiol","invnormFT4","invnormTT3","lnTT3FT4"
)

# obtaining all genetic correlation results
result_df <- data.frame()
for (thormone in thormone_list) {
  cor_df<-fread(thormone)
  cor_df<-cor_df %>% mutate(p1=gsub("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/munge_thormone/","",p1),p2=gsub("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/munge_trait/","",p2)) %>% mutate(p1=gsub(".sumstats.gz","",p1),p2=gsub(".sumstats.gz","",p2)) %>% mutate(rg=as.numeric(rg),se=as.numeric(se)) %>% mutate(Hormone=thormone)
  # computing p-values
  print(thormone)
  cor_df <- cor_df %>% select(p2,rg,se,p,Hormone)
  result_df <- rbind(result_df,cor_df)
}

# Define mappings
p2_mapping <- c(
  "ADHD2022" = "ADHD",
  "SCZ" = "Schizophrenia",
  "TS" = "Tourette's Syndrome",
  "ocs2024obsessive" = "Obsessive-Compulsive Symptoms",
  "pgc-bip2021-BDI" = "Bipolar Disorder (Type 1)",
  "pgc-bip2021-BDII" = "Bipolar Disorder (Type 2)",
  "pgc-bip2021-all" = "Bipolar Disorder (All)",
  "pgc-panic2019" = "Panic Disorder",
  "pgcAN2" = "Anorexia Nervosa",
  "depression" = "Major Depressive Disorder"
)

hormone_mapping <- c(
  "SHBG" = "SHBG",
  "estradiol" = "Estradiol",
  "invnormFT3" = "FT3",
  "invnormFT4" = "FT4",
  "invnormTSH" = "TSH",
  "invnormTT3" = "TT3",
  "lnFT3FT4" = "FT3/FT4",
  "lnTT3FT4" = "TT3/FT4",
  "testosterone" = "Testosterone"
)

# Transform dataframe
result_df <- result_df %>%
  mutate(
    p2 = recode(p2, !!!p2_mapping),
    Hormone = recode(Hormone, !!!hormone_mapping)
  )

# Print transformed dataframe
print(result_df)

hormone_order <- c("TSH", "FT3", "TT3", "FT4", "FT3/FT4", "TT3/FT4", "SHBG", "Estradiol", "Testosterone")
hormone_order <- rev(hormone_order)
result_df$Hormone <- factor(result_df$Hormone, levels = hormone_order)


# Load necessary libraries
library(ggplot2)
library(viridis)
library(dplyr)
library(scales)

# Assuming `result_df` is already loaded

# Convert p-values into a color scale
result_df <- result_df %>%
  mutate(color = case_when(
    rg > 0 ~ "#FA8072",
    rg < 0 ~ "#6495ED",
    rg == 0 ~ "white"
  )) %>%
  mutate(darkness = case_when(
    p < 0.05/90 ~ 1,
    p < 0.001 ~ 0.8,
    p < 0.05 ~ 0.6,
    TRUE ~ 0.1
  ))

# Format p-values with fixed two decimal places in scientific notation
result_df$p_formatted <- sapply(result_df$p, function(x) {
  if (x < 0.1) {
    exponent <- floor(log10(x))
    base <- formatC(x / 10^exponent, format = "f", digits = 2)
    bquote(.(base) %*% 10^.(exponent))
  } else {
    formatC(x, format = "f", digits = 2)
  }
})

# Create plot
ggplot(result_df, aes(x = p2, y = Hormone)) +
  geom_tile(aes(fill = color, alpha = darkness), color = "black") +
  scale_fill_identity() +
  scale_alpha_continuous(range = c(0.1, 1)) +
  geom_text(aes(label = p_formatted), parse = TRUE, size = 5, color = "black") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(fill = "white", color = "black", size = 1)
  )

# Save the plot
ggsave("/gpfs/data/huo-lab/BCAC/james.li/xtrait/output/rG/parsed_output/corrplot.png", width = 15, height = 13.5, dpi = 1200)

# save output file
TABLE_S3<-result_df %>% select(Hormone,p2,rg,se,p) %>% rename(`Psychiatric Disorder`=p2,`Genetic Correlation (rG)`=rg,`Standard Error of rG`=se,`P-value of rG`=p)
write.table(TABLE_S3,file="/gpfs/data/huo-lab/BCAC/james.li/xtrait/output/rG/parsed_output/TABLE_S3.tsv",quote=F,row.names=F,col.names=T,sep="\t")
```

# Analying and plotting cross-hormone genetic correlation results in R 
```{r}
library(data.table)
library(dplyr)
library(ggplot2)

########################
# HORMONE-ONLY RESULTS #
########################
# set working directory
setwd("/gpfs/data/huo-lab/BCAC/james.li/xtrait/output/rG_only_hormone/parsed_output")
thormone_list <- c(
  "SHBG","invnormFT3","invnormTSH","lnFT3FT4","testosterone","estradiol","invnormFT4","invnormTT3","lnTT3FT4"
)

# obtaining all genetic correlation results
result_df <- data.frame()
for (thormone in thormone_list) {
  cor_df<-fread(thormone)
  cor_df<-cor_df %>% mutate(p1=gsub("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/munge_thormone/","",p1),p2=gsub("/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/munge_thormone/","",p2)) %>% mutate(p1=gsub(".sumstats.gz","",p1),p2=gsub(".sumstats.gz","",p2)) %>% mutate(rg=as.numeric(rg),se=as.numeric(se)) %>% mutate(Hormone=thormone)
  # computing p-values
  print(thormone)
  cor_df <- cor_df %>% select(p1,p2,rg,se,p,Hormone)
  result_df <- rbind(result_df,cor_df)
}

hormone_mapping <- c(
  "SHBG" = "SHBG",
  "estradiol" = "Estradiol",
  "invnormFT3" = "FT3",
  "invnormFT4" = "FT4",
  "invnormTSH" = "TSH",
  "invnormTT3" = "TT3",
  "lnFT3FT4" = "FT3/FT4",
  "lnTT3FT4" = "TT3/FT4",
  "testosterone" = "Testosterone"
)

# Transform dataframe
result_df <- result_df %>%
  mutate(
    p2 = recode(p2, !!!hormone_mapping),
    p1 = recode(p1, !!!hormone_mapping)
  )

# Print transformed dataframe
print(result_df)

hormone_order <- c("TSH", "FT3", "TT3", "FT4", "FT3/FT4", "TT3/FT4", "SHBG", "Estradiol", "Testosterone")
#hormone_order <- rev(hormone_order)
result_df$p1 <- factor(result_df$p1, levels = rev(hormone_order))
result_df$p2 <- factor(result_df$p2, levels = hormone_order)

# Load necessary libraries
library(ggplot2)
library(viridis)
library(dplyr)
library(scales)

# Assuming `result_df` is already loaded

# Convert p-values into a color scale
result_df <- result_df %>%
  mutate(color = case_when(
    p < 1e-90 ~ "black",
    rg > 0 ~ "#FA8072",
    rg < 0 ~ "#6495ED",
    rg == 0 ~ "white"
  )) %>%
  mutate(darkness = case_when(
    p < 0.05/36 ~ 1,
    p < 0.001 ~ 0.8,
    p < 0.05 ~ 0.6,
    TRUE ~ 0.1
  ))

# Format p-values with superscript scientific notation and two decimal places
result_df$p_formatted <- sapply(result_df$p, function(x) {
  if (x < 0.1) {
    exponent <- floor(log10(x))
    base <- formatC(x / 10^exponent, format = "f", digits = 2)
    bquote(.(base) %*% 10^.(exponent))
  } else {
    formatC(x, format = "f", digits = 2)
  }
})

# filtering out unnecessary panels
result_df_filtered <- result_df %>% 
  filter(!(p1=="TSH" & p2=="FT3")) %>% 
  filter(!(p1=="TSH" & p2=="TT3")) %>% 
  filter(!(p1=="FT3" & p2=="TT3")) %>% 
  filter(!(p1=="TSH" & p2=="FT4")) %>% 
  filter(!(p1=="FT3" & p2=="FT4")) %>% 
  filter(!(p1=="TT3" & p2=="FT4")) %>% 
  filter(!(p1=="TSH" & p2=="FT3/FT4")) %>% 
  filter(!(p1=="FT3" & p2=="FT3/FT4")) %>% 
  filter(!(p1=="TT3" & p2=="FT3/FT4")) %>% 
  filter(!(p1=="FT4" & p2=="FT3/FT4")) %>% 
  filter(!(p1=="FT3" & p2=="TT3/FT4")) %>% 
  filter(!(p1=="TT3" & p2=="TT3/FT4")) %>% 
  filter(!(p1=="TSH" & p2=="TT3/FT4")) %>% 
  filter(!(p1=="FT4" & p2=="TT3/FT4")) %>% 
  filter(!(p1=="FT3/FT4" & p2=="TT3/FT4")) %>% 
  filter(!(p1=="TSH" & p2=="SHBG")) %>% 
  filter(!(p1=="FT3" & p2=="SHBG")) %>% 
  filter(!(p1=="TT3" & p2=="SHBG")) %>% 
  filter(!(p1=="FT4" & p2=="SHBG")) %>% 
  filter(!(p1=="FT3/FT4" & p2=="SHBG")) %>% 
  filter(!(p1=="TT3/FT4" & p2=="SHBG")) %>% 
  filter(!(p1=="TSH" & p2=="Estradiol")) %>% 
  filter(!(p1=="FT3" & p2=="Estradiol")) %>% 
  filter(!(p1=="TT3" & p2=="Estradiol")) %>% 
  filter(!(p1=="FT4" & p2=="Estradiol")) %>% 
  filter(!(p1=="FT3/FT4" & p2=="Estradiol")) %>% 
  filter(!(p1=="TT3/FT4" & p2=="Estradiol")) %>% 
  filter(!(p1=="SHBG" & p2=="Estradiol")) %>% 
  filter(!(p1=="TSH" & p2=="Testosterone")) %>% 
  filter(!(p1=="FT3" & p2=="Testosterone")) %>% 
  filter(!(p1=="TT3" & p2=="Testosterone")) %>% 
  filter(!(p1=="FT4" & p2=="Testosterone")) %>% 
  filter(!(p1=="FT3/FT4" & p2=="Testosterone")) %>% 
  filter(!(p1=="TT3/FT4" & p2=="Testosterone")) %>% 
  filter(!(p1=="SHBG" & p2=="Testosterone")) %>% 
  filter(!(p1=="Estradiol" & p2=="Testosterone"))

# Create plot
ggplot(result_df_filtered, aes(x = p2, y = p1)) +
  geom_tile(aes(fill = color, alpha = darkness), color = "black") +
  scale_fill_identity() +
  scale_alpha_continuous(range = c(0.1, 1)) +
  geom_text(aes(label = p_formatted), parse = TRUE, size = 5, color = "black") +  # Enable math parsing
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(fill = "white", color = "black", size = 1)
  )

# Save the plot
ggsave("/gpfs/data/huo-lab/BCAC/james.li/xtrait/output/rG_only_hormone/parsed_output/corrplot.png", width = 15, height = 13, dpi = 1200)

# save output file
TABLE_S2<-result_df_filtered %>% select(p1,p2,rg,se,p) %>% filter(p1!=p2) %>% rename(`Hormone #1`=p1,`Hormone #2`=p2,`Genetic Correlation (rG)`=rg,`Standard Error of rG`=se,`P-value of rG`=p)
write.table(TABLE_S2,file="/gpfs/data/huo-lab/BCAC/james.li/xtrait/output/rG_only_hormone/parsed_output/TABLE_S2.tsv",quote=F,row.names=F,col.names=T,sep="\t")
```

