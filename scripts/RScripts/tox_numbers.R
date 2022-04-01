#!/bin/env Rscript

library(argparser)


args<-arg_parser(description = "Script to make a stack plot of toxins families, alternatively you can import the function",
                 name = "",
                 hide.opts = F)
args<-add_argument(args,arg = "--input",
                   short = "-i",
                   help = "if run from shell, path to the input file, the colum names must be class(Toxin, Nontoxin,), toxin_family",
                   default = "Cgodm_consense_df.csv")
args<-add_argument(args,arg = "--number_of_individuals",
                   short = "-n",
                   help = "if run from shell, path to the input file, the colum names must be class(Toxin, Nontoxin,), toxin_family",
                   default = 6)
args<-add_argument(args,arg = "--north",
                   short = "-r",
                   help = "if run from shell, path to the input file, the colum names must be class(Toxin, Nontoxin,), toxin_family",
                   nargs = "?",
                   default = c(4,5))
argmts<-parse_args(args)


mult_vector<-function(vec){
  cmd = paste0(vec[1])
  for (i in 2:(length(vec))){
    cmd<- paste0(vec[i],'*',cmd)
  }
#  print(cmd)
  return(eval(parse(tex=cmd)))
}


Cgodm<-read.csv(argmts$input, header = T)
Cgodm<-Cgodm[,1:(4 + argmts$number_of_individuals)]
Cgodm<-Cgodm[which(Cgodm$class == "Toxin"),]
AvN<-c()
AvS<-c()
for (i in 1:nrow(Cgodm)){
  AvN<-c(AvN,mean(as.numeric(Cgodm[i,c(argmts$north + 3)])))
  AvS<-c(AvS,mean(as.numeric(Cgodm[i,c(3+c(1:argmts$number_of_individuals)[-c(argmts$north)])])))
}
print(paste0(rep("*",50)))
print("missing in North population")
print(length(which(AvN == 0)))
print("missing in South population")
print(length(which(AvS == 0)))

print(paste0(rep("*",50)))
print("Toxin missing in all the individuals")
length(which(Cgodm[,ncol(Cgodm)] == 0))

print(paste0(rep("*",50)))
print("Toxins missing by individual")
for (i in 4:c(3 + argmts$number_of_individuals)){
  print(colnames(Cgodm)[i])
  print(length(which(Cgodm[,i] == 0)))
}

core<-0
for (i in 1:nrow(Cgodm)){
  if(mult_vector(as.numeric(Cgodm[i,4:c(3 + argmts$number_of_individuals)])) != 0){
    core<-core +1
  }
}
print(paste0(rep("*",50)))
print("core genes")
print(core)

print(paste0(rep("*",50)))
print("genes of North population")
print(nrow(Cgodm)+length(which(AvS == 0)))

print(paste0(rep("*",50)))
print("genes of South population")
print(nrow(Cgodm)-length(which(AvN == 0)))

