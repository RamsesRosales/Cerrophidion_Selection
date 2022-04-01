#!/bin/env Rscript
#TranCompPlots for every conbination of individuals
library(zCompositions)
library(argparser)
library(tidyverse)
source('~/Documents/bin/Venomancer/Guide/PlottingFunctions1.R')

p<- arg_parser('Use the functions from the script PlottingFunctions1 in ToxCodAn to get the plots of comparison between species')
p<- add_argument(p,short = '-i', arg = '--input', default = 'Cgodm_consense_df.csv', 
                 help = 'the name of the input file')
p<- add_argument(p,short = '-a', arg = '--average', default = 'NO', 
                 help = 'if NO, then it returns the comparison between individuals of the species, if YES, it also return each individual against the average')
argu<- parse_args(p)
Cerro_df<-read.csv(paste(argu$input),header = T)
Cerro_df<-Cerro_df[which(Cerro_df$Average != 0),]
b<-colnames(Cerro_df)
d<-length(b)
TPM_df <- Cerro_df %>% mutate_if(is.character,as.factor)
TPM_df <- df_clean(Cerro_df, class="class", toxin_family="toxin_family", colors=toxin_colors)
TPM_df2 <- t(cmultRepl(t(TPM_df[,4:d]),output = "p-counts"))
TPM_df2 <- cbind(TPM_df[,1:3], TPM_df2)
rownames(TPM_df2) <- TPM_df2$gene_id

a<-colnames(TPM_df2)

if (argu$average == 'NO'){
  for (i in 4:(length(a)-2)){
    print(i)
    for (j in (i+1):(length(a)-1)){
      pdf(paste0(a[i],'vs',a[j],"_Comparison.pdf"),width = 16, height=10)
      TransCompPlot(id1 = paste(a[i]),id2 = paste(a[j]))
      dev.off()
      x=c(i,j)
      print(x)
    }
  }
}else{
  for (i in 4:(length(a)-2)){
    print(i)
    for (j in (i+1):(length(a)-1)){
      pdf(paste0(a[i],'vs',a[j],"_Comparison.pdf"),width = 16, height=10)
      TransCompPlot(id1 = paste(a[i]),id2 = paste(a[j]))
      dev.off()
      x=c(i,j)
      print(x)
    }
  }
  
  for (i in 4:(length(a)-1)){
    pdf(paste0(a[i],'vs','Average',"_Comparison.pdf"),width = 16, height=10)
    TransCompPlot(id1 = paste(a[i]),id2 = paste(a[length(a)]))
    dev.off()
  }
}






