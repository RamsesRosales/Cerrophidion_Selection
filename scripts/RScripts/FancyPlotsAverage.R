#!/bin/env Rscript
#fancy plots for the species Average

library(argparser)
library(tidyverse)

p<- arg_parser('Use the functions from the script PlottingFunctions1 in ToxCodAn to get the plots for the species average from a csv with colums for each individual')
p<- add_argument(p,short = '-i', arg = '--input', default = 'No', 
                 help = 'the name of the input file')
p<- add_argument(p,short = '-n', arg = '--name', default = 'Cgpdm', 
                 help = 'the name of the output, <name>_transcriptome.pdf')
p<- add_argument(p, short = '-p', arg = '--plot_function', default = '/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/PlottingFunctions1.R',
                 help = 'the path to the ploting functions in toxcodan, note that I modify that file
                 to work with my data frame, changes might be neded to run this code in the normal ploting function')

argu<- parse_args(p)

Cerro_df<-read.csv(paste(argu$input),header = T)

Cerro_df<-Cerro_df[which(Cerro_df$Average != 0),]
Cerro_df1<-data.frame('col1' = rep('relleno',length(Cerro_df$Average)), Cerro_df)

#TPM_df <- df_clean(TPM_df, class="class", toxin_family="toxin_family", colors=toxin_colors)

source(paste(argu$plot_function))

pdf(paste0(argu$name,"_transcriptome.pdf"),width = 16, height=10)
FancyFigure(df = Cerro_df1, id = 'Average')
dev.off()  
