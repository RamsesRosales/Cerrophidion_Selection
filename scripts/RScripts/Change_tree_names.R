#!/bin/env Rscript

library(argparser)

ar<-arg_parser(description = "change the names of a phylogenetic tree",name = "Rename_tree")

ar<-add_argument(ar,arg = "--input",short = "-i",help = "input tree in newick format",
                 default = "Combined_PLA2s_trim.aln.trim.fasta.contree")

ar<-add_argument(ar,arg = "--names",short = "-n",
                 help= "csv with the current names and the new names, the colums should be named, <tree> for the old names and <hm> for the new names ",
                 default = "PLA2_names_tree_n_hm_1.csv")
ar<-add_argument(ar,arg = "--output",short = "-o",
                 help= "the name to add to the output, if not specified it will use the name of the input+ _new_names.tree",
                 default = "No")
ar<-add_argument(ar,arg = "--get_names",short = "-g",
                 help= "it will get you a list of the names as output, is easier to edit than edit the tree file",
                 default = F)
ar<-add_argument(ar,arg = "--ordered_names",short = "-r",
                 help= "use a list with the names in the same order than the names in the tree to change, you can get the name list order with -g T, it assumes the list have a header and use the first colum as it reads as a csv",
                 default = F)
ars<-parse_args(ar)

#script

library(readr)
library(dplyr)
library(ape)

print(paste0("reading tree ",ars$input))
tree<-read.tree(ars$input)

if (ars$get_names){
  print("writing csv with names and quiting")
  x<-tree[["tip.label"]]
  output<-data.frame("names" = x)
  write.csv(output,"tree_names.csv",row.names = F)
  q()
}

if (ars$ordered_names){
  print("using list with ordered names as reference")
  data_n<-read.csv(ars$names)
  tree[["tip.label"]]<-data_n[,1]
  
  print(paste0("writing ouput file",ars$names))
  
  if(ars$output == "No"){
    print(paste0("writing ouput file ",ars$input,".new_names.tree"))
    write.tree(phy=tree,paste0(ars$input,".new_names.tree"))
  }
  if(ars$output != "No"){
    print(paste0("writing ouput file ",ars$output,"_new_names.tree"))
    write.tree(tree,paste0(ars$output,"_new_names.tree"))
  }
  q()
}

print(paste0("reading names data base ",ars$names))

data_n<-read_csv(ars$names)

#for (i in 1:nrow(data_n)){
#  print(paste0(data_n$old_hm[i]," , ",data_n$old_tree [i]))
#  print(grep(data_n$old_hm,data_n$old_tree [i]))
#}

print(paste0("replacing the names",ars$names))

x<-tree[["tip.label"]]
for (i in 1:length(x)){
  if(length(which(x[i] == data_n$tree)) != 0){
    x[i]<-data_n$hm[which(x[i] == data_n$tree)]
  }
}

tree[["tip.label"]]<-x

print(paste0("writing ouput file",ars$names))

if(ars$output == "No"){
  print(paste0("writing ouput file ",ars$input,".new_names.tree"))
  write.tree(phy=tree,paste0(ars$input,".new_names.tree"))
}
if(ars$output != "No"){
  print(paste0("writing ouput file ",ars$output,"_new_names.tree"))
  write.tree(tree,paste0(ars$output,"_new_names.tree"))
}


