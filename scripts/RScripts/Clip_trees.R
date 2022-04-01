#!/bin/env Rscript

#Script to keep or clip branches of a tree, in newick format

if ("ape" %in% installed.packages()[,1] == F){
  install.packages("ape", repo = "https://cran.microsoft.com/")
}

library(ape)
library(argparser)

argu<-arg_parser(description = "This is a small program to clip phylogenetic trees in base of the labesls
           the program uses library ape and argparser",name = "tree cliper ")

argu<-add_argument(argu,arg = "--input",short = "-i",help = "the name of the tree to cut in newik format",default = "test.tree")

argu<-add_argument(argu,arg = "--label_list", short = "-l",help = "list with the labels to keep or clip it has to ve in a txt file in one colum list", default = "test.list")

argu<-add_argument(argu,arg = "--mode", short = "-m", help = "keep to keep the tips in the list, clip to cut the tips in the list", default = "keep")

argu<-add_argument(argu,arg = "--output", short = "-o", help = "name of the resulting newick tree, the extention .tree will be added", default = "test_keep")

argu<-add_argument(argu,arg = "--pattern", short = "-p", help = "in case there wass a easy to fix substitution between your names and the tip labels put the pattern and sub in a list, 
                   with the pattern as the first element and the substitution in the element the change will be made in the tips, not in your list", default = "No sub")

argum<-parse_args(argu)



tree<-read.tree(argum$input)
clip_tips<-scan(argum$label_list,what = "character",sep = "\n")

if (argum$pattern != "No sub"){
  pat_sub<-scan(argum$pattern,what = "character",sep = "\n")
  tree[["tip.label"]]<-gsub(paste(pat_sub[1]),paste(pat_sub[2]),tree[["tip.label"]])  
}


if(argum$mode == "keep"){
  cut_tree<-keep.tip(tree,clip_tips)
}
if(argum$mode == "clip"){
  cut_tree<- drop.tip(tree,clip_tips)
}

write.tree(cut_tree,paste0(argum$output,".tree"))


