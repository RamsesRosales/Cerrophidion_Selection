#!/bin/env Rscript

library(argparser)
library(tidyverse)

p<- arg_parser('Use the functions from the script PlottingFunctions1 in ToxCodAn to get the plots from a csv with stacked colums withowt a average')
p<- add_argument(p,short = '-i', arg = '--input', default = 'No', 
                 help = 'the name of the input file')
p<- add_argument(p, short = '-n', arg = '--names', default = 'No',
                 help = 'a list with the name of the individuals')
p<- add_argument(p, short = '-p', arg = '--plot_function', default = '/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/PlottingFunctions1.R',
                 help = 'the path to the ploting functions in toxcodan, note that I modify that file
                 to work with my data frame, changes might be neded to run this code in the normal ploting function')

argu<- parse_args(p)

Cerro_df<-read.csv(paste(argu$input),header = T)

Cerro_df<-Cerro_df[which(Cerro_df$TPM != 0),]

#TPM_df <- df_clean(TPM_df, class="class", toxin_family="toxin_family", colors=toxin_colors)

#library(zCompositions)

h<-(read.table(paste(argu$names), header = F))
dirlist<-c(h$V1)

source(paste(argu$plot_function))


for (i in 1:length(dirlist)){
  pdf(paste0(dirlist[i],"_transcriptome.pdf"),width = 16, height=10)
  FancyFigure(df = Cerro_df[which(Cerro_df$ID == paste(dirlist[i])),], id = 'TPM')
  dev.off()
}
