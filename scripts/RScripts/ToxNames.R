#!/bin/env Rscript


library(argparser)

ar<-arg_parser(description = "program to create a csv with one colum with original name of the toxins and other colum with the name made of the toxin family and the rank of the toxin by average expression in deacreasing order, it will also take the file of the expected count for each inidividual and will modify the names based on the old phyle",
           name = "ToxNames")
ar<-add_argument(ar,arg = "--input",short = "-i",
                 help = "name of a csv file with the average expression the gene.id and the toxin family in colums, generated with consense.R",
                 default = "No")
ar<-add_argument(ar,arg = "--modify",short = "-m",
                 help = "if you want to create a new version of a file with the new names give the path of the file here", 
                 default = "No" )

ar<-add_argument(ar,arg = "--species",short = "-s",
                 help = "name to attach to the output, and to the start of the names, <name>_ToxNames.csv",
                 default = "Cerrophidion")

ars<-parse_args(ar)

#if (ars$input == "No"){
#  print("No file selected, give a input file")
#  quit()
#}

library(readr)
library(dplyr)

input<-read_csv(ars$input)

input<-input[which(input$class == "Toxin"),]

input<-input[order(input$Average,decreasing = T),]

x<-rep(" ",nrow(input))
ToxNam<-data_frame(gene_id = x, ToxNam = x)

for (i in 1:nrow(input)){
  ToxNam[i,1]<-input$gene_id[i]
  ToxNam[i,2]<-paste0(ars$species,'_',input$toxin_family[i],"_",i)
}

write_csv(ToxNam,paste0(ars$species,"_ToxNames.csv"))

if (ars$modify == "No"){
  print("table with the names relashipship generated, if you want to add the new names to a file give the path using option -m")
}
if (ars$modify != "No"){
  original<-read_csv(ars$modify)
  write_csv(original,paste0(gsub(".csv","",ars$modify),"_original.csv"))
  modifyed<-original
  for (i in 1:nrow(ToxNam)){
    modifyed$gene_id[which(ToxNam$gene_id[i] == modifyed$gene_id)]<-ToxNam$ToxNam[i]
  }
  write_csv(modifyed,ars$modify)
}



