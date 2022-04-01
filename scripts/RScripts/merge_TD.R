#!/bin/env Rscript
library(argparser)
library(tidyverse)

p<- arg_parser('Merge the RESULTS of the TajimasD per gene, synonymous Tajimas D, nonsynonymous Tajimas D',
               'merge Tajimas D')
p<- add_argument(p,short = '-g', arg = '--genes', default = 'Genes.txt', 
                 help = 'the name or identifier of the individual that was sequenced')
p<- add_argument(p,short = '-e', arg = '--length', default = 'length.txt', 
                 help = 'a list with the lenght of each gene, I got this from the Combined.vcf using grep and should coincide in order with Genes.txt')
p<- add_argument(p, short = '-l', arg = '--list', default = 'TDlist',
                 help = 'the list with the name of the files in the RESULTS directory that contain Tajimas D for each gene')
p<- add_argument(p, short = '-o', arg = '--output', default = 'Cgodm',
                 help = 'the name that you whant in the output <output>_TajimasD.csv')
p<- add_argument(p, short = '-s', arg = '--op_system', default = 'mac',
                 help = 'the operative system between linux of mac os, this is only to selec the proper separator of new line for the list')
argu<- parse_args(p)

#### def functions



### the header for each of the TajimasD 
headers<-c('CHROM','BIN_START','N_SNPS','TajimaD')
identifiers<- scan(paste(argu$list), character(), quote = '')

if (argu$op_system == 'mac'){
  genes<-scan(paste(argu$genes), character(), quote = '',sep = '\r')
}
if (argu$op_system == 'linux'){
  genes<-scan(paste(argu$genes), character(), quote = '',sep = '\n')
}
if (argu$op_system == 'windows'){
  genes<-scan(paste(argu$genes), character(), quote = '',sep = '\r\n')
}

len<-scan(paste(argu$length), character(), quote = '')


for (i in 1:length(identifiers)){
  if (length(grep('perGene',identifiers[i])) != 0){
    TDperGene<- read.table(paste0(identifiers[i]),header = F,sep = '\t')
    colnames(TDperGene)<-headers
  }
  if (length(grep('Nonsynonymous',identifiers[i])) != 0){
    TDNonsynonymous<- read.table(paste0(identifiers[i]),header = F,sep = '\t')
    colnames(TDNonsynonymous)<-headers
  }
  if (length(grep('Synonymous',identifiers[i])) != 0){
    TDSynonymous<- read.table(paste0(identifiers[i]),header = F,sep = '\t')
    colnames(TDSynonymous)<-headers
  }
}

T_NT<-c()

for (i in 1:length(genes)){
  if( length(grep('TOXIN',genes[i])) != 0) {
    T_NT<-c(T_NT,'Toxin')
  }
  else{
    T_NT<-c(T_NT,'Nontoxin')
  }
}

TPM_df2 <- data.frame(gene.id = genes)
TPM_df2 <- TPM_df2 %>% mutate( class = T_NT)
d <-c("3FTx","AChe","AChE","BPP","CNP","CRISP","CTL","Ficolin","FusedToxin","GluCyc",
      "Goannatyrotoxin","HYAL","KUN","LAAO","MP","MYO","NGF","NUC","PDE","PLA2","PLB",
      "SVMMP","SVMPI","SVMPII","SVMPIII","SVSP","TruncHYAL","VDP4","VEGF",
      "VESP","Vespryn","VF","Waprin")

ToxFam1<-rep('NT',length(TPM_df2$gene.id))

for (i in 1:length(TPM_df2$gene.id)){
  if (TPM_df2$class[i] == 'Toxin'){
    ToxFam1[i]<- 'uncharacterised'
    for (j in d){
      if( length(grep(paste(j),TPM_df2$gene.id[i])) != 0) {
        ToxFam1[i]<- paste(j)
      }
      else{
      }
    }
  }
  else{
  }
}

TPM_df2 <- TPM_df2 %>%
mutate(toxin_family = ToxFam1, length = len )
       


blank<- rep(0,length(genes))
output<-data.frame(gene.id = genes,
                   class = TPM_df2$class,
                   toxin_family = TPM_df2$toxin_family,
                   length = TPM_df2$length,
                   N_SNP = blank,
                   TajimasD = blank,
                   Ns_N_SNP = blank,
                   Ns_TajimasD = blank,
                   S_N_SNP = blank,
                   S_TajimasD = blank )

for (i in 1:nrow(TDperGene)){
  output$N_SNP[which(output$gene.id == TDperGene$CHROM[i])]<- TDperGene$N_SNPS[i]
  output$TajimasD[which(output$gene.id == TDperGene$CHROM[i])]<- TDperGene$TajimaD[i]
}
for (i in 1:nrow(TDNonsynonymous)){
  output$Ns_N_SNP[which(output$gene.id == TDNonsynonymous$CHROM[i])]<- TDNonsynonymous$N_SNPS[i]
output$Ns_TajimasD[which(output$gene.id == TDNonsynonymous$CHROM[i])]<- TDNonsynonymous$TajimaD[i]
}
for (i in 1:nrow(TDSynonymous)){
  output$S_N_SNP[which(output$gene.id == TDSynonymous$CHROM[i])]<- TDSynonymous$N_SNPS[i]
output$S_TajimasD[which(output$gene.id == TDSynonymous$CHROM[i])]<- TDSynonymous$TajimaD[i]
}

write.csv(output,paste0(argu$output,"_TajimasD.csv"), row.names = F)





  
