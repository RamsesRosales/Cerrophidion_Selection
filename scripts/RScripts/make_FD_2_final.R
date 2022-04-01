#!/bin/env Rscript
library(argparser)

p<- arg_parser('Merge the <name>_TajimasD.csv, the <name>_consense_df.csv, the <name>_Dif_Exp_tab.csv',
               'merge Tajimas D')
p<- add_argument(p,short = '-t', arg = '--TD_df', default = 'Cgodm_TajimasD.csv', 
                 help = 'The Tajimas de data frame path')
p<- add_argument(p, short = '-r', arg = '--tpm_df', default = 'Cgodm_consense_df.csv',
                 help = 'the list with the name of the files in the RESULTS directory that contain Tajimas D for each gene')
p<- add_argument(p, short = '-d', arg = '--Dif_Exp', default = 'Cgodm_DifExp_tab.csv',
                 help = 'the list with the name of the files in the RESULTS directory that contain Tajimas D for each gene')
p<- add_argument(p, short = '-o', arg = '--output', default = 'Cgodm',
                 help = 'the name that you whant in the output <output>_finaldata.csv')
argu<- parse_args(p)


#### load libraries

library(tidyverse)

#### define functions

count_names_len<-function(data){
  longer_name<-c(0)
  smaller_name<-c(1000)
  nids<-c()
  nidl<-c()
  for (i in 1:length(data)){
    if (nchar(data[i]) > longer_name){
      longer_name <- nchar(data[i])
      nidl<- i
    }
    if (nchar(data[i]) < smaller_name){
      smaller_name <- nchar(data[i])
      nids<-i
    }
  }
  result<-data.frame('min/max' = c("min", "max"),"position" = c(nids,nidl), "length" = c(smaller_name, longer_name))
  print(result)
}

check_equality<-function(col1,col2){
  x<- col1
  y<- col2
  positions<-c(0)
  for (i in 1:length(x)){
    if(x[i] != y[i]){
      positions<-c(positions,i)
    }
  }
  print(positions)
}

##### load data frames
TDdf<-read.csv(paste(argu$TD_df),header = T)
Tpmdf<-read.csv(paste(argu$tpm_df), header = T)
Dif_E<-read.csv(paste(argu$Dif_Exp), header = T)

if (nrow(TDdf) != nrow(Tpmdf)){
  TDdf1<-TDdf
  print("modify TDdf to add the missing genes witout values")
  TDdf2<-data.frame(gene.id = Tpmdf$gene_id,class = Tpmdf$class, toxin_family = Tpmdf$toxin_family)
  TDdf2[,4:26]<-0
  for (i in 1:nrow(TDdf2)){
    if(length(which(TDdf2$gene.id[i] == TDdf$gene.id) !=0 )){
      TDdf2[i,2:26]<-TDdf[which(TDdf2$gene.id[i] == TDdf$gene.id),2:26]
    }
  }
  colnames(TDdf2)<- colnames(TDdf)
  TDdf<-TDdf2
}


x<-count_names_len(TDdf$gene.id)
y<-count_names_len(Tpmdf$gene_id)
w<-count_names_len(Dif_E$X)

if (x$length[1] != y$length[1] & x$length[2] != y$length[2]){
  print("check that names of gene.id are the same")
  q(save = F)
}

if (w$length[1] != y$length[1] & w$length[2] != y$length[2]){
  print("check that names of gene.id are the same")
  #q(save = F)
}

#colnames(Tpmdf)[1]<- 'genes.id'
colnames(Tpmdf)[2]<- 'class1'
colnames(Tpmdf)[3]<- 'toxinfamily1'


Tpmdf<- Tpmdf %>% arrange(gene_id)
TDdf<- TDdf %>% arrange(gene.id)
Dif_E<- Dif_E %>% arrange(X)
output<-cbind(TDdf,Tpmdf,Dif_E)

a<-check_equality(col1 = output$gene.id, col2 = output$gene_id)
b<-check_equality(output$class,output$class1)
c<-check_equality(output$toxin_family,output$toxinfamily1)
d<-check_equality(output$gene.id,output$X)

if (length(a) > 1 | length(b) > 1 |length(c) > 1 | length(d) > 1){
  print('check that the class and toxin family coincide in both data frames')
  q(save = F)
}

sdev_tpm<-c()
for (i in 1:nrow(output)){
  tmp<-sd(output[i,grep(paste0(argu$output,'.CLP'),colnames(output))])
  sdev_tpm[i]<- tmp
}



output$sdev_tpm <- sdev_tpm
output$bcv_tpm <- output$sdev_tpm/output$Average

ind<- colnames(output) == 'Pop'
e<-length(colnames(output)[ind])

if (e > 0){
  ind<-colnames(output) == 'Pop'
  colnames(output)[ind]<- 'Dif_Exp'
}
if (e == 0){
  ind<-colnames(output) == 'SVL'
  colnames(output)[ind]<- 'Dif_Exp'
}



ind<-grep(paste0(argu$output,'.CLP'),colnames(output))
x<-output[,ind]

output<- select(output, gene_id, class, toxin_family,
                Average,sdev_tpm,bcv_tpm,Pi,
                Dif_Exp,length,Variant_rate,
                VI_High,VI_Low,VI_Moderate,
                VE_Missense, VE_Nonsense,
                Ns_N_SNP,S_N_SNP,N_SNP,
                TajimasD,
                Ns_TajimasD,S_TajimasD,SPvNP,
                LRT,p_value,,CW1,CW2,CW3,UCW1,UCW2,UCW3)



output<-data.frame(output[1:3],x,output[4:ncol(output)])


z<-colnames(output)[grep(paste0(argu$output,'.CLP'),colnames(output))]
y<-c('ID',"class","toxin_family",z,
     "average_tpm","stdev_tpm","bcv_tpm",'pi',"diff_exprs","length","Variant_rate",'High_impact','Low_impact','Moderate_impact',
     'Missense','Nonsense','Total_Nonsynonymous','Synonymous','Total_variants','tajimasD','tajimasD_nonsyn','tajimasD_synon',
     'SPvNP',"LRT","p_value","CW1","CW2","CW3","UCW1","UCW2", "UCW3")

colnames(output)<- y

ind<- output$diff_exprs == 0
output$diff_exprs[ind]<- 'FALSE'
ind<- output$diff_exprs == 1
output$diff_exprs[ind]<- 'TRUE'


write.csv(output,paste0(argu$output,'_finaldata.csv'), row.names = F)
