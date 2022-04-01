#!/bin/env Rscript
#fancy plots for the species Average North and South

library(argparser)
library(tidyverse)

p<- arg_parser('Use the functions from the script PlottingFunctions1 in ToxCodAn to get the plots for the species average from a csv with colums for each individual')
p<- add_argument(p,short = '-i', arg = '--input', default = 'No', 
                 help = 'the name of the input file')
p<- add_argument(p,short = '-l', arg = '--list', default = 'No', 
                 help = 'list with the names of the North population in the first line and the Name of the South population in the second line
                 separated by a ","; the names must be the same that colum names for example in Cgodm.CLP2359, must like that and not Cgodm-CLP2359')
p<- add_argument(p,short = '-n', arg = '--name', default = 'Cgodm', 
                 help = 'the name of the output, <name>_transcriptome.pdf')

argu<- parse_args(p)

if(argu$list == "No"){
  North<- c("Cgodm.CLP2377","Cgodm.CLP2378")
  South<- c("Cgodm.CLP2359","Cgodm.CLP2360","Cgodm.CLP2362","Cgodm.CLP2388")
}else{
  North<-scan(argu$list,nlines = 1,what = "character",sep = ",")
  South<-scan(argu$list,skip = 1, what = "character", sep = ",")
}

print("North individuals:")
print(North)
print("South individuals:")
print(South)


Cerro_df<-read.csv(paste(argu$input),header = T)

Average_North<-c()
Average_South<-c()
for (i in 1:nrow(Cerro_df)){
  Average_North<-c(Average_North,mean(as.numeric(Cerro_df[i,North])))
  Average_South<-c(Average_South,mean(as.numeric(Cerro_df[i,South])))
}

Cerro_df<-mutate(Cerro_df,Average_North = Average_North, Average_South = Average_South)


Cerro_dfN<-Cerro_df[which(Cerro_df$Average_North != 0),]
Cerro_dfS<-Cerro_df[which(Cerro_df$Average_South != 0),]
Cerro_dfN<-data.frame('col1' = rep('relleno',length(Cerro_dfN$Average_North)), Cerro_dfN)
Cerro_dfS<-data.frame('col1' = rep('relleno',length(Cerro_dfS$Average_South)), Cerro_dfS)

#TPM_df <- df_clean(TPM_df, class="class", toxin_family="toxin_family", colors=toxin_colors)

#source('~/Documents/bin/Venomancer/Guide/PlottingFunctions1.R')
source('/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/PlottingFunctions1.R')


pdf(paste0(argu$name,"AverageNorth_transcriptome.pdf"),width = 16, height=10)
FancyFigure(df = Cerro_dfN, id = 'Average_North')
dev.off()  

pdf(paste0(argu$name,"AverageSouth_transcriptome.pdf"),width = 16, height=10)
FancyFigure(df = Cerro_dfS, id = 'Average_South')
dev.off()


