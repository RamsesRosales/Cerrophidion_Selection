#!/bin/env Rscript

#Script to rename fasta headers with a pattern or a reference table


library(argparser)
argu<-arg_parser(description = "This is a small program to rename fasta headers, it takes the fasta file and a pattern
                 and rename the headers with consecutive number with that patter, the ouput is the renamed fasta and
                 a csv with the olname and the new names",name = "Rename Fasta Headers ")

argu<-add_argument(argu,arg = "--input",short = "-i",help = "the name of the fasta file to read",default = "test.fasta")

argu<-add_argument(argu,arg = "--pattern", short = "-p", help = "the pattern you want to use, e.g. if you use header it will be header1, header2 etc.", default = "No")

argu<-add_argument(argu,arg = "--ref_table", short = "-t", help = "in case you already have a reference table, the name of the file to read, it should be a csv file
                   with the first colum as old name and second colum as new name, the code assume the table has headers", default = "No")
argu<-add_argument(argu,arg = "--ordered", short = "-o", 
                   help = "set as True in case you have a file with the new names in order, but the old names dont match exactly", 
                   default = F)
argu<-add_argument(argu,arg = "--unordered", short = "-u", 
                   help = "set as True in case you have a file with the new names unordered and with weird characters as '|,[]'", 
                   default = F)

argum<-parse_args(argu)

###
print("load or install and load libraries phylotools and seqinr")

if ("phylotools" %in% installed.packages()[,1] == F){
  install.packages("phylotools", repo = "https://cran.microsoft.com/")
}

if ("seqinr" %in% installed.packages()[,1] == F){
  install.packages("seqinr", repo = "https://cran.microsoft.com/")
}

library(phylotools)
library(seqinr)


print("read fasta file")
data<-read.fasta(file = argum$input,seqtype = "DNA")

fasta_name<-argum$input
fasta_name<-gsub(".fasta","",fasta_name)

if (argum$pattern != "No" & argum$ref_table != "No"){
  print("Chose pattern or ref_table, the program cannot run both")
  quit()
}

if (argum$pattern != "No"){
  print("make referenced headers table")
  ref_table<-data.frame("old_names" = getName(data))
  new_names<- c()
  for (i in 1:nrow(ref_table)){
    new_names<-c(new_names,paste0(argum$pattern,'_',i))
  }
  ref_table$new_names<- new_names
}

if (argum$ref_table != "No"){
  print("read referenced headers table")
  ref_table<-read.csv(paste0(argum$ref_table),header = T)
}

if (argum$ordered){
  for (i in 1:nrow(ref_table)){
    attr(data[[i]],which = "name") <- ref_table$new_names[i]
  }
  print("rename fasta headers")
  write.fasta(data,file = paste0(fasta_name,"_renamed.fasta"),names = getName.list(data),nbchar = max(getLength.list(data)))
  quit()
}
if (argum$unordered){
  print("old headers unordered and with weird characters, will try to rename")
  print("this only work if colnames are old_names,new_names")
  checked<-(1:length(getName.list(data)))
  for (i in 1:length(getName.list(data))){
    if(length(which(ref_table$old_names == getName(data[[i]]))) != 0){
      ind<-which(ref_table$old_names == getName(data[[i]]))
      attr(data[[i]],which = "name")<-ref_table$new_names[ind]
      checked<-checked[which(checked != i)]
    }
  }
  for (i in checked){
    if(length(grep(getName(data[[i]]),ref_table$old_names)) == 1){
      ind<-grep(getName(data[[i]]),ref_table$old_names)
      attr(data[[i]],which = "name")<-ref_table$new_names[ind]
      checked<-checked[which(checked != i)]
    }
  }
  for (i in checked){
    ind<-grep(getName(data[[i]]),ref_table$old_names)
#    if(length(ind) < length(getName.list(data))){
      print("************************************************")
      print(getName(data[[i]]))
      print(paste0(length(ind)," matches for name ",i))
      print("name unmodified")
#    }
    print("************************************************")
  }
  print(paste0("rename ", length(getName.list(data))-length(checked) ," fasta headers, ",length(checked)," names unchanged"))
  write.fasta(data,file = paste0(fasta_name,"_renamed.fasta"),names = getName.list(data),nbchar = max(getLength.list(data)))
  quit()
}


print("rename fasta headers")
rename.fasta(paste(argum$input),ref_table,outfile = paste0(fasta_name,"_renamed.fasta"))

if (argum$pattern != "No" & argum$ref_table == "No"){
  print("print reference table")
  write.csv(ref_table,paste0(fasta_name,"ref_table.csv"),row.names = F)
}
