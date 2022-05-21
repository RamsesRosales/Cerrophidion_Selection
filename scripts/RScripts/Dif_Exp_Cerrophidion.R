#!/bin/env Rscript
##### program to do a DEanalysis
library(argparser)

argu<-arg_parser(description = 'Differential expression with DeSeq2 and EDGER, written for Cerrophidion', name = 'DifExp')

argu<-add_argument(parser = argu, arg = "--input", short = '-i',help = 'input file, should be a csv with the expected counts for each individual from rsem',default = 'Cgodm_consense_ExpCount_df.csv')

argu<-add_argument(parser = argu, arg = '--metadata', short = '-m', help = 'metadata file, the file with the info for each individual in this case we need Pop, and SVL', default = '~/Desktop/Cerrophidion_Testing/Cerrophidion_specimens.csv' )

argu<-add_argument(argu, arg = '--species', short = '-s', help = 'name of the species in the tabla, example Cgodm', default = 'Cgodm')

argu<-add_argument(argu, arg = '--treatment', short = '-t', help = 'set the name of the groups to compare as they are in the table', default = c('Pop','SVL'),nargs = '?')

argu<-add_argument(argu, arg = '--ToxFam', short = '-TF', help = 'set as True if you want to run analysis of Toxin Families too', default = T)

argmts<- parse_args(argu)

print('loading libraries')

# Read Libraries and Functions.
library(readr)
library(edgeR)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(ashr)

print('libraries loaded')

print('loading metadata')

# Prepare Metadata, because we don't have replicates on other variables, we only will test for north and south population on Cgodm and SVL in all Cerrophidion.
metadata <- read_csv(paste(argmts$metadata))
metadata <- metadata %>% filter(Species == paste(argmts$species)) %>% 
  mutate_if(is.character, as.factor) %>% 
  column_to_rownames("Individuals") %>% 
  select(-Collector, -DataLocation)

metadata<-metadata[sort(rownames(metadata)),]

print('metadata loaded')

print('import expected count data')
# Prepare Expected Counts
expcounts <- read_csv(paste(argmts$input))
#expcounts <- expcounts %>% filter(class=="Toxin") we did not run this as we get different results if we use all the genes and only toxins.
x<-ncol(expcounts)-1
tmp_df<-as.matrix(expcounts[,4:x])
rownames(tmp_df)<-expcounts$gene_id

colnames(tmp_df)<-gsub(paste0(argmts$species,"-"),"",colnames(tmp_df))


#DESeq2 only takes intergers so we have to do this.z
for (i in 1:ncol(tmp_df)){
  tmp_df[,i]<-as.integer(tmp_df[,i])
}

#prepare the variables to run the loop.
varnames<-argmts$treatment
fndf<- c()
for (i in 1:length(varnames)){
  fndf<-c(fndf,paste0('DESeq2_',varnames[i]))
}

#prepare the matrix with the proper desing with the Pop and SVL variables
for (i in 1:length(varnames)){
  cmd = paste0(fndf[i],"<-DESeqDataSetFromMatrix(tmp_df,metadata, design = ~",varnames[i],")")
  eval(parse(text=cmd))
}

dir.create('Dif_Exp')

#run the loop to get the DE of genes considering all the genes and the normal DESeq2 not the Likelihod Ratio Test (LRT)
for (i in 1:length(fndf)){
  DESeq2<-get(eval(paste0(fndf[i])))
  keep <- rowSums(counts(DESeq2)) >= 10
  DESeq2 <- DESeq2[keep,]
  rm(keep)
  #  DESeq2<-DESeq(DESeq2,fitType = "local",test = 'LRT', reduced = ~1)
  DESeq2<-DESeq(DESeq2,fitType = "local")
  DESeq2_res <- results(DESeq2,alpha=0.05)
  #  print(head(DESeq2_res))
  print(summary(DESeq2_res))
  DESeq2_res1<-lfcShrink(DESeq2, coef = paste(resultsNames(DESeq2)[2]), type = 'ashr')
  plotMA(DESeq2_res1)
  title(main = paste(resultsNames(DESeq2)[2]))
  DESeq2_res <- as.data.frame(DESeq2_res)
  table<-DESeq2_res[DESeq2_res$padj<0.05,]
  table<-table[complete.cases(table),c(2,6)]
  #  table<-head(table[complete.cases(table),c(2,6)])
  #  print(table)
  assign(paste0('tableDESeq2',varnames[i]), table)
  write.csv(table,paste0('Dif_Exp/',argmts$species,'_tableDESeq2',varnames[i],'.csv'),col.names = T,row.names = T)
  rm(DESeq2)
}


# Prepare Expected Counts
expcounts <- read_csv(paste(argmts$input))
#expcounts <- expcounts %>% filter(class=="Toxin") we did not run this as we get different results if we use all the genes and only toxins.
tmp_df<-as.matrix(expcounts[,4:x])
rownames(tmp_df)<-expcounts$gene_id

EDGER_Pop<-function(){
  group<-factor(c(metadata$Pop), labels = c('North','South'))
  # Prepare edgeR Dataframe
  edgeR_df <- DGEList(counts=tmp_df,group=group)
  keep <- rowSums(cpm(edgeR_df)>1) >= 2
  #keep1 <- filterByExpr(edgeR_df, group = group)
  #print(keep)
  edgeR_df <- edgeR_df[keep, , keep.lib.sizes=FALSE]
  rm(keep)
  edgeR_df <- calcNormFactors(edgeR_df)
  design1 <- model.matrix(~group)
  edgeR_df <- estimateDisp(edgeR_df,design1)
  fit1 <- glmFit(edgeR_df,design1)
  S_vs_N <- glmLRT(fit1,coef = 2)
  summary(decideTests(S_vs_N))
  plotMD(S_vs_N)
  resSN <- as.data.frame(topTags(S_vs_N, n=Inf, p.value=0.05))
  table1 <- resSN[,c(1,5)]
  write.csv(table1,paste0('Dif_Exp/',argmts$species,'_tableEDGER_Pop.csv'),col.names = T,row.names = T)
}

EDGER_SVL<-function(){
  edgeR_df <- DGEList(counts=tmp_df,samples = metadata)
  keep <- rowSums(cpm(edgeR_df)>1) >= 2
  edgeR_df <- edgeR_df[keep, , keep.lib.sizes=FALSE]
  rm(keep)
  edgeR_df <- calcNormFactors(edgeR_df)
  
  ###########
  design <- model.matrix(~SVL, data=edgeR_df$samples)
  dispersion <- estimateDisp(edgeR_df,design)
  fit <- glmFit(dispersion,design)
  SVLcomp<-glmLRT(fit)
  summary(decideTests(SVLcomp))
  plotMD(SVLcomp)
  resSVL<-as.data.frame(topTags(SVLcomp,n=Inf,p.value = 0.05))
  resSVL
  table2 <- resSVL[,c(1,5)]
  write.csv(table2,paste0('Dif_Exp/',argmts$species,'_tableEDGER_SVL.csv'),col.names = T,row.names = T)
}


if(length(grep('Pop',argmts$treatment)) == 1){
  EDGER_Pop()
}

if(length(grep('SVL',argmts$treatment)) == 1){
  EDGER_SVL()
}

EDGER_PopTF<-function(){
  group<-factor(c(metadata$Pop), labels = c('North','South'))
  # Prepare edgeR Dataframe
  edgeR_df <- DGEList(counts=tmp_Tox,group=group)
  keep <- rowSums(cpm(edgeR_df)>1) >= 2
  #keep1 <- filterByExpr(edgeR_df, group = group)
  #print(keep)
  edgeR_df <- edgeR_df[keep, , keep.lib.sizes=FALSE]
  rm(keep)
  edgeR_df <- calcNormFactors(edgeR_df)
  design1 <- model.matrix(~group)
  edgeR_df <- estimateDisp(edgeR_df,design1)
  fit1 <- glmFit(edgeR_df,design1)
  S_vs_N <- glmLRT(fit1,coef = 2)
  summary(decideTests(S_vs_N))
  plotMD(S_vs_N)
  resSN <- as.data.frame(topTags(S_vs_N, n=Inf, p.value=0.05))
  table1 <- resSN[,c(1,5)]
  write.csv(table1,paste0('Dif_Exp/',argmts$species,'_tableEDGER_TF_Pop.csv'),col.names = T,row.names = T)
}

EDGER_SVLTF<-function(){
  edgeR_df <- DGEList(counts=tmp_Tox,samples = metadata)
  keep <- rowSums(cpm(edgeR_df)>1) >= 2
  edgeR_df <- edgeR_df[keep, , keep.lib.sizes=FALSE]
  rm(keep)
  edgeR_df <- calcNormFactors(edgeR_df)
  
  ###########
  design <- model.matrix(~SVL, data=edgeR_df$samples)
  dispersion <- estimateDisp(edgeR_df,design)
  fit <- glmFit(dispersion,design)
  SVLcomp<-glmLRT(fit)
  summary(decideTests(SVLcomp))
  plotMD(SVLcomp)
  resSVL<-as.data.frame(topTags(SVLcomp,n=Inf,p.value = 0.05))
  resSVL
  table2 <- resSVL[,c(1,5)]
  write.csv(table2,paste0('Dif_Exp/',argmts$species,'_tableEDGER_TF_SVL.csv'),col.names = T,row.names = T)
}


expcountsTF<-as.data.frame(expcounts %>% group_by(toxin_family,class) %>%
                             summarize_if(is_numeric,sum))
y<-ncol(expcountsTF)-1
expcountsTF <- expcountsTF %>% filter(class=="Toxin")
tmp_df_TF<-as.matrix(expcountsTF[,3:y])
rownames(tmp_df_TF)<-expcountsTF$toxin_family

tmp_NT<-expcounts <- expcounts %>% filter(class=="Nontoxin")
NT_rownames<- tmp_NT$gene_id
tmp_NT<-as.matrix(expcounts[,4:x])
rownames(tmp_NT)<- NT_rownames

tmp_Tox<- as.matrix(rbind(tmp_df_TF,tmp_NT)) #I rbind the data set of the other genes to the toxin families sum as it change the results as the fit use all the genes.

colnames(tmp_Tox)<-gsub(paste0(argmts$species,"-"),"",colnames(tmp_Tox))

ToxinFamilies<- function(){
  #########repetition of all for toxin families.
  #DESeq2 only takes intergers so we have to do this.z
  for (i in 1:ncol(tmp_Tox)){
    tmp_Tox[,i]<-as.integer(tmp_Tox[,i])
  }
  
  #prepare the matrix with the proper desing with the Pop and SVL variables
  for (i in 1:length(varnames)){
    cmd = paste0(fndf[i],'_TF',"<-DESeqDataSetFromMatrix(tmp_Tox,metadata, design = ~",varnames[i],")")
    eval(parse(text=cmd))
  }
  
  for (i in 1:length(fndf)){
    DESeq2<-get(eval(paste0(fndf[i],'_TF')))
    keep <- rowSums(counts(DESeq2)) >= 10
    DESeq2 <- DESeq2[keep,]
    rm(keep)
    #  DESeq2<-DESeq(DESeq2,fitType = "local",test = 'LRT', reduced = ~1)
    DESeq2<-DESeq(DESeq2,fitType = "local")
    DESeq2_res <- results(DESeq2,alpha=0.05)
    #  print(head(DESeq2_res))
    print(summary(DESeq2_res))
    DESeq2_res1<-lfcShrink(DESeq2, coef = paste(resultsNames(DESeq2)[2]), type = 'ashr')
    plotMA(DESeq2_res1)
    title(main = paste(resultsNames(DESeq2)[2]))
    DESeq2_res <- as.data.frame(DESeq2_res)
    table<-DESeq2_res[DESeq2_res$padj<0.05,]
    table<-table[complete.cases(table),c(2,6)]
    #  table<-head(table[complete.cases(table),c(2,6)])
    #  print(table)
    assign(paste0('tableDESeq2_TF',varnames[i]), table)
    rm(DESeq2)
    write.csv(table,paste0('Dif_Exp/',argmts$species,'_tableDESeq2_TF_',varnames[i],'.csv'),col.names = T,row.names = T)
  }
  
  
  tmp_Tox<- as.matrix(rbind(tmp_df_TF,tmp_NT))
  
  
  if(length(grep('Pop',argmts$treatment)) == 1){
    EDGER_PopTF()
  }
  
  if(length(grep('SVL',argmts$treatment)) == 1){
    EDGER_SVLTF()
  }
  
  
}


if (argmts$ToxFam){
  ToxinFamilies()
}


















