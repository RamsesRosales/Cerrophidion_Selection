#!/bin/env Rscript
######

library(argparser)

argu<-arg_parser(description = 'Take the TPM values of the rsem to do a heat map, it can take a table with 0 and 1 to show the diferential expression', name = 'DifExp_HeatMap')

argu<-add_argument(parser = argu, arg = "--input", short = '-i',help = 'input file, should be a csv with the expected counts for each individual from rsem',default = '../Cgodm_consense_df.csv')

argu<-add_argument(parser = argu, arg = '--metadata', short = '-m', help = 'metadata file, the file with the info for each individual in this case we need Pop, and SVL', default = '~/Desktop/Cerrophidion_Testing/Cerrophidion_specimens.csv' )

argu<-add_argument(parser = argu, arg = "--DifExp", short = '-d',help = 'Set as True if you have want to show DifExp in the rows',default = 'Cgodm_DifExp_tab.csv')

argu<-add_argument(argu, arg = '--species', short = '-s', help = 'name of the species in the tabla, example Cgodm, will be used to name the output', default = 'Cgodm')

argu<-add_argument(argu, arg = '--Pop', short = '-p', help = 'set as True if you want nave Pop colum', default = T)

argu<-add_argument(argu, arg = '--ToxFam', short = '-T', help = 'set as True if you want to run analysis of Toxin Families too', default = T)

argu<-add_argument(argu, arg = '--TF_input', short = '-t', help = 'name of the file to use as input if ToxFam is T, example Cgodm_DifExp_TF_tab.csv', default = 'Cgodm_DifExp_TF_tab.csv')

argu<-add_argument(argu, arg = '--colors', short = '-c', help = 'set as true to use a color file', default = T)

argu<-add_argument(argu, arg = '--Width', short = '-w', help = 'width of the output', default = 1200)

argu<-add_argument(argu, arg = '--Height', short = '-e', help = 'height of the output', default = 1000)

argu<-add_argument(argu, arg = '--WidthTF', short = '-z', help = 'width of the output', default = 1200)

argu<-add_argument(argu, arg = '--HeightTF', short = '-y', help = 'height of the output', default = 1000)

argmts<- parse_args(argu)


#argmts$input<-'../Ctzot_consense_df.csv'
#argmts$DifExp<-'Ctzot_DifExp_tab.csv'
#argmts$species<-'Ctzot'
#argmts$TF_input<-'Ctzot_DifExp_TF_tab.csv'
#argmts$Pop<-F
#argmts$colors<-T

source('/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/PlottingFunctions1.R')
source('/zfs/venom/Ramses/bin/RamsesScripts/Colors_HM_Cerrophidion_1.R')
library(pheatmap)
library(dplyr)
library(readr)  
library(tibble)
library(zCompositions)
library(ggplot2)
library(dendextend)
library(RColorBrewer)



#####def functions

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


if (argmts$colors){
  my_color <- get(eval(paste0('my_color_',argmts$species)))
  my_colorTF <-get(eval(paste0('my_colorTF_',argmts$species)))
}
 


##### prepare data


#### metadata for colums
metadata <- read_csv(paste0(argmts$metadata))
metadata <- metadata %>% filter(Species == paste(argmts$species)) %>% 
  mutate_if(is.character, as.factor) %>% 
  column_to_rownames("Individuals") #%>% 
#  select(-Collector, -DataLocation)

metadata<-metadata[,-c(12,13)]

##### data for rows 
D_E<-read.csv(paste0(argmts$DifExp), header = T)
rownames(D_E)<- D_E$X
D_E<- D_E[,-1]
ind<-grep('TOXIN',rownames(D_E))
D_E<- D_E[ind,]


#### data values
lnToxins<- read_csv(paste0(argmts$input))
torm<-c()
for (i in grep('TOXIN',lnToxins$gene_id)){
    if (lnToxins$Average[i] == 0){
    torm<-c(torm,lnToxins$gene_id[i])
  }
}


lnToxins<-lnToxins[which(lnToxins$Average != 0),]
TPM_df <- df_clean(lnToxins, class="class", toxin_family="toxin_family", colors=toxin_colors)
a<-ncol(lnToxins)
TPM_df2 <- t(cmultRepl(t(TPM_df[,4:a]),output = "p-counts"))
TPM_df2 <- cbind(TPM_df[,1:3], TPM_df2)
rownames(TPM_df2) <- TPM_df2$gene_id
lnToxins<- TPM_df2 %>% filter(class == 'Toxin') %>%
  mutate_if(is.numeric,log) #%>%
#  arrange(desc(Average)) 
rn<-c()
for (i in 1:nrow(lnToxins)){
  rn<-c(rn,paste0('toxin',i))
}
rownames(lnToxins)<-rn
x<-colnames(lnToxins)
x<-gsub(paste0(argmts$species,'-'),'',x)
colnames(lnToxins)<-x
y<-as.matrix(lnToxins[,4:(a-1)])

##### clean D_E dt from toxins with average 0


if (length(torm) >= 1){
  for (i in 1:length(torm)){
    assign(paste0('ind',i),rownames(D_E) == torm[i])
  }
  ind<-c(rep(F,length(rownames(D_E))))
  for (i in 1:length(torm)){
    for (j in 1:length(get(eval(paste0('ind',i)))))
      if (get(eval(paste0('ind',i)))[j] == T){
        ind[j]<- T
      }
  }
  D_E<- D_E[!ind,]
}



##### sort the data to show in the map
D_E1<-data.frame(ifelse(D_E$SVL== 1,'p<0.5_SVL','NS'))
colnames(D_E1)[1]<-c('SVL_DE')
rownames(D_E1)<-rownames(lnToxins)

if (argmts$Pop){
  b<-data.frame(ifelse(D_E$Pop== 1,'p<0.5_Pop','NS'))
  D_E1<-data.frame(D_E1,b)
  colnames(D_E1)[2]<-c('Pop_DE')
}

D_E1<-data.frame(D_E1,lnToxins$toxin_family)
colnames(D_E1)[ncol(D_E1)]<-c('ToxFam')


c<-cbind(y,D_E1)
c<-arrange(c,ToxFam)

data<-c[,1:ncol(y)]
HM_RN<-c[,(ncol(y)+1):ncol(c)]


##### code run

if (argmts$colors){
  z<-pheatmap(data,cluster_rows = F,cluster_cols = T,annotation_row = HM_RN, annotation_col = metadata[,3:6], annotation_legend = T, show_rownames = F,  annotation_colors = my_color, border_color = NA)
}
if (!argmts$colors){
  z<-pheatmap(data,cluster_rows = F,cluster_cols = T,annotation_row = HM_RN, annotation_col = metadata[,3:6], annotation_legend = T, show_rownames = F, border_color = NA)
}


save_pheatmap_png(z,paste0(argmts$species,'_heatmap.png'),height=as.numeric(argmts$Height),width=as.numeric(argmts$Width))


##### Here start the proces for ToxFam

if (argmts$ToxFam){
  lnToxins<- read_csv(paste0(argmts$input))
  lnToxins<-lnToxins[which(lnToxins$Average != 0),]
  TPM_df <- df_clean(lnToxins, class="class", toxin_family="toxin_family", colors=toxin_colors)
  a<-ncol(lnToxins)
  TPM_df2 <- t(cmultRepl(t(TPM_df[,4:a]),output = "p-counts"))
  TPM_df2 <- cbind(TPM_df[,1:3], TPM_df2)
  rownames(TPM_df2) <- TPM_df2$gene_id
  lnToxins<- TPM_df2 %>% filter(class == 'Toxin') %>%
    group_by(toxin_family,class) %>% 
    summarize_if(is.numeric,sum) 
  
  TPM_family_df <- as.data.frame(lnToxins)
  rownames(TPM_family_df) <- TPM_family_df$toxin_family
  
  lnClasses <- as.data.frame(TPM_family_df) %>%
    mutate_if(is.numeric,log)
  
  rownames(lnClasses) <- lnClasses$toxin_family
  
  x<-colnames(lnClasses)
  x<-gsub(paste0(argmts$species,'-'),'',x)
  colnames(lnClasses)<-x
  b<-length(lnClasses)-1
  y<-as.matrix(lnClasses[,3:b])
  
  D_E<-read.csv(paste0(argmts$TF_input), header = T)
  rownames(D_E)<- D_E$X
  D_E<- D_E[,-1]
  
  D_E1<-data.frame(ifelse(D_E$SVL== 1,'p<0.5_SVL','NS'))
  colnames(D_E1)[1]<-c('SVL_DE')
  rownames(D_E1)<-rownames(lnToxins)
  
  if (argmts$Pop){
    b<-data.frame(ifelse(D_E$Pop== 1,'p<0.5_Pop','NS'))
    D_E1<-data.frame(D_E1,b)
    colnames(D_E1)[2]<-c('Pop_DE')
  }
  
  
  c<-cbind(y,D_E1)
  
  data<-c[,1:ncol(y)]
  HM_RN<-c[,(ncol(y)+1):ncol(c)]
  HM_RN<-data.frame(HM_RN)
  rownames(HM_RN)<- rownames(c)
  if(ncol(HM_RN) == 1){
    colnames(HM_RN)<-'SVL_DE'
  }
  
  if (argmts$colors){
    z<-pheatmap(data,cluster_rows = F,cluster_cols = T,annotation_row = HM_RN,annotation_col = metadata[,3:6], annotation_legend = T, show_rownames = T, annotation_colors = my_colorTF, border_color = NA )
  }
  if (!argmts$colors){
    z<-pheatmap(data,cluster_rows = F,cluster_cols = T,annotation_row = HM_RN,annotation_col = metadata[,3:6], annotation_legend = T, show_rownames = T, border_color = NA )
  }
  
  save_pheatmap_png(z,paste0(argmts$species,'_TF_heatmap.png'),height=as.numeric(argmts$HeightTF),width=as.numeric(argmts$WidthTF))
}





