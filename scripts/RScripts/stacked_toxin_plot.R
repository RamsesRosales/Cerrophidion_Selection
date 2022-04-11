library(argparser)


args<-arg_parser(description = "Script to make a stack plot of toxins families, alternatively you can import the function\n
                 by sourcing the script, to run it from the shell you have to set -f F",
           name = "",
           hide.opts = F)
args<-add_argument(args,arg = "--functions",
                   short = "-f",
                   help = "By default set as true, if sourced the script it only import the function and the color file",
                   default = T)
args<-add_argument(args,arg = "--input",
                   short = "-i",
                   help = "if run from shell, path to the input file, the colum names must be class(Toxin, Nontoxin,), toxin_family",
                   default = "Cgodm_consense_df.csv")
args<-add_argument(args,arg = "--output",
                   short = "-o",
                   help = "if run from shell, path and name of the output file, <path/name>_stacked_barplot.svg",
                   default = "Cgodm")
args<-add_argument(args,arg = "--number_of_individuals",
                   short = "-n",
                   help = "number of individuals",
                   default = 6)
args<-add_argument(args,arg = "--tf_substitution",
                   short = "-t",
                   help = "list with a substitution of toxin_families with the name of the gene firts and the new value separated by a ',', you can have several substitutions",
                   default = "No")
args<-add_argument(args,arg = "--color_class",
                   short = "-c",
                   help = "if toxin family substitution, you can set the color for your new class, <color>,<toxin_familr>",
                   default = "No")
args<-add_argument(args,arg = "--height",
                   short = "-e",
                   default = 15,
                   help = "hight of the output plot default is 15>",
                   type = "numeric")
args<-add_argument(args,arg = "--width",
                   short = "-w",
                   default = 25,
                   help = "width of the ouput the default is 25",
                   type = "numeric")
argmts<-parse_args(args)


if(argmts$functions == F){
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("load libraries")
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(RColorBrewer)
  
  
  #####
  # Toxin Colors
  toxin_colors<-c("#d30b94","3FTx",
                  "#201923","BPP",
                  "#772b9d","CRISP",
                  "#0ec434","CTL",
                  "#a4a4a4","Ficolin",
                  "#ffffff","FusedToxin",
                  "#29bdab","HYAL",
                  "#ffdab9","KUN",#e68f66
                  "#ffc413","LAAO",
                  "#fcff5d","Lacta",
                  "#f47a22","MYO",
                  "#fffac8","NGF",#ffcba5
                  "#5d4c86","NUC",
                  "#7cd676","PDE",
                  "#f52727","PLA2",
                  "#991919","PLA2_neuro",
                  "#c3a5b4","PLB",
                  "#3998f5","SVMMP",
                  "#277da7","SVMP",
                  "#277da7","SVMPI",
                  "#3750db","SVMPII",
                  "#2f2aa0","SVMPIII",
                  "#228c68","SVSP",
                  "#aaffc3","VEGF",#946aa2
                  "#aa6e28","Vespryn",#c56133
                  "#cccccc","VF",
                  "#f07cab","Waprin",
                  "#8B5A2B","AChE",
                  "#235b54","uncharacterised",
                  "#37294f","VDP4",
                  "#cd0000","PLA2_pseudoneuro",
                  "#632819","zOpenColor5",
                  "#b732cc","zOpenColor6",
                  "#8ad8e8","zOpenColor7")
  
  ##### ADD extra colors
  if (argmts$color_class != "No"){
    x<-scan(argmts$color_class,what = "character",sep = "\n")
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("adding colors from list:")
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print(x)
    for (i in 1:length(x)){
      y<-strsplit(x,split = ",")
      toxin_colors<-c(toxin_colors,y[[1]][1],y[[1]][2])
    }
    print(toxin_colors)
  }
  #####
  toxin_colors_df<-matrix(toxin_colors,nrow=length(toxin_colors)/2,ncol=2,byrow=T)
  toxin_colors_df<-as.data.frame(cbind(toxin_colors_df,rep(1,nrow(toxin_colors_df))))
  toxin_colors_df$V2 <- factor(toxin_colors_df$V2, levels = toxin_colors_df$V2)
  toxin_colors_df$V3 <- as.numeric(toxin_colors_df$V3)
  #ggbarplot(toxin_colors_df, "V2","V3", fill=toxin_colors_df$V1, width = 1, xlab="",ylab="", main="Toxin Colors") + rotate_x_text(angle = 45)
  toxin_colors<-palette(as.vector(toxin_colors_df$V1))
  toxin_colors<-palette(as.vector(toxin_colors_df$V1))
  names(toxin_colors)<-toxin_colors_df$V2
  
  
  ###### set functions
  
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("loading functions:")
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  
  percentage<-function(datum,data){
    datum <- as.numeric(datum)
    data<- as.vector(data)
    if(datum == 0){
      print("percentage of 0")
      return(0)
    }
    if(datum < 0){
      print("this function is only for positive numbers")
    }
    if(datum > 0 & sum(data) > 0){
      msg<-paste0("percentage of",datum," from ",sum(data),": ",(datum/sum(data))*100,"%")
      print(msg)
      return(datum/sum(data))
    }
  }
  stacked_barplot<-function(data = data,toxin_family = 'toxin_family',class = 'class',n_samples = 6,print=F){
    co<- n_samples+1
    print("pass1")
    data<- data %>% filter(!!as.name(class) == "Toxin") %>%
      mutate_if(is.character,as.factor) %>%
      group_by(!!as.name(toxin_family)) %>%
      summarise_if(is.numeric,sum) %>%
      ungroup()
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("control stop:")
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%") 
    data <- data %>% gather(key = 'individual',value = "TPM", 2:!!as.numeric(co)) %>%
      filter(TPM != 0) 
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("File correctly formatted")
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%") 
    
    a<-ggplot(data, aes(fill=toxin_family, y=TPM, x=individual)) + 
      geom_bar(position="fill", stat="identity", width = 1, colour = "black", size = 0.5 )+
      theme_classic()+
      scale_fill_manual(values = toxin_colors)+
      ylab(label = "Toxin Family Expression %")
    
    if(print==TRUE){
      print(a)
    }
    if(print==FALSE){
      return(a)
    }
  }
  
  #######
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("reading data from ")
  print(argmts$input)
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  
  
  data<-read_csv(argmts$input)
#  data<-read_csv("~/Desktop/Cerrophidion/Cgodm_1/17rsem/Cgodm_consense_df.csv")
  
  if (argmts$tf_substitution != "No"){
    x<-scan(argmts$tf_substitution,what = "character",sep = "\n")
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("substititution of toxin_families from list:")
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%") 
    for (i in 1:length(x)){
      y<-strsplit(x[i],split = ",")
      print(y[[1]])
      data$toxin_family[which(data$gene_id == y[[1]][1])] <- y[[1]][2]
      y<-c()
    }
  }
  

  a<- unique(data$toxin_family)[which(unique(data$toxin_family) != "NT")] 
  tmp<- c()
  for (i in a){
    tmp<-c(tmp, toxin_colors[which(names(toxin_colors) == i)])
  }
  tmp<-tmp[order(names(tmp))]
  toxin_colors<-tmp
  
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print("Running Function:")
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%") 
  a<-stacked_barplot(data = data,toxin_family = 'toxin_family',class = 'class',n_samples = argmts$number_of_individuals,print=F)
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print(paste0("Saving File: ",argmts$output,"_stacked_barplot.svg"))
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%") 

  ggsave(paste0(argmts$output,"_stacked_barplot.svg"),plot = a,height = argmts$height,width = argmts$width,units = "cm")
}

if (argmts$functions == T){
  rm(args)
  rm(argmts)
  toxin_colors<-c("#d30b94","3FTx",
                  "#201923","BPP",
                  "#772b9d","CRISP",
                  "#0ec434","CTL",
                  "#a4a4a4","Ficolin",
                  "#ffffff","FusedToxin",
                  "#29bdab","HYAL",
                  "#ffdab9","KUN",#e68f66
                  "#ffc413","LAAO",
                  "#fcff5d","Lacta",
                  "#f47a22","MYO",
                  "#fffac8","NGF",#ffcba5
                  "#5d4c86","NUC",
                  "#7cd676","PDE",
                  "#f52727","PLA2",
                  "#991919","PLA2_neuro",
                  "#c3a5b4","PLB",
                  "#3998f5","SVMMP",
                  "#277da7","SVMP",
                  "#277da7","SVMPI",
                  "#3750db","SVMPII",
                  "#2f2aa0","SVMPIII",
                  "#228c68","SVSP",
                  "#aaffc3","VEGF",#946aa2
                  "#aa6e28","Vespryn",#c56133
                  "#cccccc","VF",
                  "#f07cab","Waprin",
                  "#8B5A2B","AChE",
                  "#235b54","uncharacterised",
                  "#37294f","VDP4",
                  "#cd0000","PLA2_pseudoneuro",
                  "#632819","zOpenColor5",
                  "#b732cc","zOpenColor6",
                  "#8ad8e8","zOpenColor7")
  #####
  toxin_colors_df<-matrix(toxin_colors,nrow=length(toxin_colors)/2,ncol=2,byrow=T)
  toxin_colors_df<-as.data.frame(cbind(toxin_colors_df,rep(1,nrow(toxin_colors_df))))
  toxin_colors_df$V2 <- factor(toxin_colors_df$V2, levels = toxin_colors_df$V2)
  toxin_colors_df$V3 <- as.numeric(toxin_colors_df$V3)
  #ggbarplot(toxin_colors_df, "V2","V3", fill=toxin_colors_df$V1, width = 1, xlab="",ylab="", main="Toxin Colors") + rotate_x_text(angle = 45)
  toxin_colors<-palette(as.vector(toxin_colors_df$V1))
  toxin_colors<-palette(as.vector(toxin_colors_df$V1))
  names(toxin_colors)<-toxin_colors_df$V2
  print("only source stacked_barplot_function")
  stacked_barplot<-function(data = data,toxin_family = 'toxin_family',class = 'class',n_samples = 6,print=TRUE){
    co<- n_samples+1
    data<- data %>% filter(!!as.name(class) == "Toxin") %>%
      mutate_if(is.character,as.factor) %>%
      group_by(!!as.name(toxin_family)) %>%
      summarise_if(is.numeric,sum) %>%
      ungroup()
    data <- data %>% gather(key = 'individual',value = "TPM", 2:co) %>%
      filter(TPM != 0) 
    
    a<-ggplot(data, aes(fill=toxin_family, y=TPM, x=individual)) + 
      geom_bar(position="fill", stat="identity", width = 1, colour = "black", size = 0.5 )+
      theme_classic()+
      scale_fill_manual(values = toxin_colors)
    if(print==TRUE){
      print(a)
    }
    if(print==FALSE){
      return(a)
    }
  }
  
}



#data$toxin_family[which(data$gene_id == "TOXIN_extenderContig120||Cgodm-CLP2360trinity_PLA2-3a_trinityContig14327_PLA2")] <- "PLA2_neuro"
#data$toxin_family[which(data$gene_id == "TOXIN_extenderContig449||Cgodm-CLP2360trinity_PLA2-2a_trinityContig14324_PLA2")] <- "PLA2_neuro"
#data$toxin_family[which(data$gene_id == "TOXIN_extenderContig167||Cgodm-CLP2377ext_PLA2-4_extContig372_PLA2")] <- "PLA2_pseudoneuro"

#to get the data base with percentage run these

#data<-data_original

#for(i in 1:nrow(data)){
#  for(j in 2:ncol(data)){
#    print(data$toxin_family[i]) 
#    print(colnames(data)[j])
#    print(c(i,j))
#    data[i,j]<- percentage(data_original[i,j],data_original[,j])
#  }
#}

