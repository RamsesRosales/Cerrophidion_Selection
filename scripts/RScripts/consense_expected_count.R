#!/bin/env Rscript
library(argparser)
library(tidyverse)

p<- arg_parser('take the data frame of the previous program and make a data frame as the one in the example of the ToxCodAn program')
p<- add_argument(p, short = '-i', arg = '--input', default = 'No file selected',
                 help = 'the file with the data frame done with the other script')
p<- add_argument(p, short = '-s', arg = '--Species', default = 'Species_consense',
                 help = 'give the name of the output <name>_consense_df.csv')

argu<- parse_args(p)

print(paste('input data frame ',argu$input))

a<- read.csv(paste(argu$input), header = T)

b<-unique(a$ID)
aa<- a[which(a$ID == b[1]),]
tpm_df <- select(aa,gene_id,class,toxin_family)

for (i in 1:length(b)){
  x <- c(a$expected_count[which(a$ID == b[i])])
  x <- x[1:nrow(tpm_df)]
  tpm_df<- tpm_df %>% mutate(!!b[i] := x)
}
i<-1
z<- c()
for (i in 1:nrow(tpm_df)){
  x <- tpm_df[i,4:ncol(tpm_df)]
#  y <- levels(as.factor(x))
  z <- c(z,round(mean(as.numeric(x)),4))
}

tpm_df<-tpm_df %>% mutate(Average = z)

write.csv(tpm_df,paste0(argu$Species,'_consense_ExpCount_df.csv'), row.names = F)

#tpm_df<- tpm_df %>% select(-Average)
