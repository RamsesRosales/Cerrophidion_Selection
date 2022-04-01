R\_SCRIPTS:Venom variation and evolution in populations of montane
pitvipers (*Viperidae* :*Cerrophidion*).
================
Ramses A. Rosales-Garcia, Rhett M. Rautsaw , Erich P. Hofmann, Christoph
I. Grunwald, Jason M. Jones, Hector Franz-Chavez, Ivan T.
Ahumada-Carrillo, Ricardo Ramirez-Chaparro, Miguel Angel De la
Torre-Loranca, Jason L. Strickland, Andrew J. Mason, Matthew L. Holding,
Miguel Borja, Gamaliel Castaneda-Gaytan, Darin R. Rokyta, Tristan D.
Schramer, N. Jade Mellor, Edward A. Myers, Christopher Parkinson
2022 April 01

-   [R analysis](#r-analysis)
    -   [rename toxins by average expression and toxin
        families](#rename-toxins-by-average-expression-and-toxin-families)
    -   [snpEFF Data](#snpeff-data)
        -   [Is there a significant association between the type of
            mutation (nonsynonymous vs. synonymous) and the transcript
            class (toxin
            vs. nontoxin)?](#is-there-a-significant-association-between-the-type-of-mutation-nonsynonymous-vs-synonymous-and-the-transcript-class-toxin-vs-nontoxin)
        -   [Is there a significant difference in the number of SNPs per
            kilobase between toxin and
            nontoxins?](#is-there-a-significant-difference-in-the-number-of-snps-per-kilobase-between-toxin-and-nontoxins)
        -   [Is there a significant difference in the number of non
            synonymous SNPs per kilobase between toxin and
            nontoxins?](#is-there-a-significant-difference-in-the-number-of-non-synonymous-snps-per-kilobase-between-toxin-and-nontoxins)
        -   [Is there a significant difference in the number of
            synonymous SNPs per kilobase between toxin and
            nontoxins?](#is-there-a-significant-difference-in-the-number-of-synonymous-snps-per-kilobase-between-toxin-and-nontoxins)
    -   [Nucleotide diversity](#nucleotide-diversity)
        -   [Is there a significant difference in nucleotide diversity
            between toxins and
            nontoxins?](#is-there-a-significant-difference-in-nucleotide-diversity-between-toxins-and-nontoxins)
        -   [Does nucleotide diversity correlate with average expression
            levels?](#does-nucleotide-diversity-correlate-with-average-expression-levels)
        -   [Does nucleotide diversity correlate with average expression
            levels **in
            toxins**?](#does-nucleotide-diversity-correlate-with-average-expression-levels-in-toxins)
        -   [Does nucleotide diversity correlate with average expression
            levels and toxin family **in
            toxins**?](#does-nucleotide-diversity-correlate-with-average-expression-levels-and-toxin-family-in-toxins)
    -   [Tajima’s D](#tajimas-d)
        -   [Are estimates of Tajima’s D significantly different than 0
            (suggesting selection) for toxins or
            nontoxins?](#are-estimates-of-tajimas-d-significantly-different-than-0-suggesting-selection-for-toxins-or-nontoxins)
        -   [Is there a significant difference in Tajima’s D between
            toxins and
            nontoxins?](#is-there-a-significant-difference-in-tajimas-d-between-toxins-and-nontoxins)
        -   [Figure 2A](#figure-2a)
        -   [Does Tajima’s D (All, Positive, and Negative values)
            correlate with average expression
            levels?](#does-tajimas-d-all-positive-and-negative-values-correlate-with-average-expression-levels)
        -   [Does Tajima’s D (All, Positive, and Negative values)
            correlate with average expression levels **in
            toxins**?](#does-tajimas-d-all-positive-and-negative-values-correlate-with-average-expression-levels-in-toxins)
        -   [Does Tajima’s D (All, Positive, and Negative values)
            correlate with average expression levels and toxin families
            **in
            toxins**?](#does-tajimas-d-all-positive-and-negative-values-correlate-with-average-expression-levels-and-toxin-families-in-toxins)
        -   [Figure 3A](#figure-3a)
    -   [Global Fst](#global-fst)
        -   [Is there a significant difference in Global Fst between
            toxins and
            nontoxins?](#is-there-a-significant-difference-in-global-fst-between-toxins-and-nontoxins)
        -   [Figure 2B](#figure-2b)
        -   [Does Global Fst correlate with average expression
            levels?](#does-global-fst-correlate-with-average-expression-levels)
        -   [Does Global Fst correlate with average expression levels
            **in
            toxins**?](#does-global-fst-correlate-with-average-expression-levels-in-toxins)
        -   [Does Global Fst correlate with average expression levels
            and toxin family **in
            toxins**?](#does-global-fst-correlate-with-average-expression-levels-and-toxin-family-in-toxins)
        -   [Figure 3B](#figure-3b)
    -   [Tajima’s D Synonymous
        Mutations](#tajimas-d-synonymous-mutations)
        -   [Are estimates of Tajima’s D for synonymous mutations
            significantly different than 0 (suggesting selection) for
            toxins or
            nontoxins?](#are-estimates-of-tajimas-d-for-synonymous-mutations-significantly-different-than-0-suggesting-selection-for-toxins-or-nontoxins)
        -   [Is there a significant difference in Tajima’s D for
            synonymous mutations between toxins and
            nontoxins?](#is-there-a-significant-difference-in-tajimas-d-for-synonymous-mutations-between-toxins-and-nontoxins)
    -   [Tajima’s D Nonsynonymous
        Mutations](#tajimas-d-nonsynonymous-mutations)
        -   [Are estimates of Tajima’s D for nonsynonymous mutations
            significantly different than 0 (suggesting selection) for
            toxins or
            nontoxins?](#are-estimates-of-tajimas-d-for-nonsynonymous-mutations-significantly-different-than-0-suggesting-selection-for-toxins-or-nontoxins)
        -   [Is there a significant difference in Tajima’s D for
            nonsynonymous mutations between toxins and
            nontoxins?](#is-there-a-significant-difference-in-tajimas-d-for-nonsynonymous-mutations-between-toxins-and-nontoxins)
        -   [Does nonsynonymous Tajima’s D (All, Positive, and Negative
            values) correlate with average expression levels and toxin
            families **in
            toxins**?](#does-nonsynonymous-tajimas-d-all-positive-and-negative-values-correlate-with-average-expression-levels-and-toxin-families-in-toxins)
    -   [HyPhy Busted](#hyphy-busted)
        -   [Is there a significant difference in BUSTED between toxins
            and
            nontoxins?](#is-there-a-significant-difference-in-busted-between-toxins-and-nontoxins)
            -   [Non Parametric test](#non-parametric-test)
            -   [SLM](#slm)
        -   [Figure 2C](#figure-2c)
        -   [Do BUSTED Likelihood Ratios correlate with average
            expression
            levels?](#do-busted-likelihood-ratios-correlate-with-average-expression-levels)
            -   [Non parametric](#non-parametric)
            -   [Non parametric test](#non-parametric-test-1)
            -   [SLM](#slm-1)
            -   [Non parametric](#non-parametric-1)
            -   [SLM](#slm-2)
        -   [Do BUSTED Likelihood Ratios correlate with average
            expression levels **in
            toxins**?](#do-busted-likelihood-ratios-correlate-with-average-expression-levels-in-toxins)
            -   [Non parametric](#non-parametric-2)
            -   [SLM](#slm-3)
            -   [SLM](#slm-4)
        -   [Figure 3C](#figure-3c)
    -   [COWPLOTS](#cowplots)

# R analysis

Read the data

``` r
### install required packages
libraries<- c('readr','ggplot2','dplyr','cowplot','kableExtra','tidyr',"Hmisc")

for (i in libraries){
  if (!(i %in% installed.packages()[,1])){
    install.packages(i)
  }
}


#load libraries
library(readr)
#library(readxl)
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(cowplot)
library(kableExtra)
```

    ## 
    ## Attaching package: 'kableExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

``` r
library(tidyr)
library(compositions)
```

    ## Welcome to compositions, a package for compositional data analysis.
    ## Find an intro with "? compositions"

    ## 
    ## Attaching package: 'compositions'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     anova, cor, cov, dist, var

    ## The following objects are masked from 'package:base':
    ## 
    ##     %*%, norm, scale, scale.default

``` r
#read data
Data<-read_csv("Cgodm_finaldata.csv")
```

    ## Rows: 5541 Columns: 36

    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (3): ID, class, toxin_family
    ## dbl (31): Cgodm.CLP2359, Cgodm.CLP2360, Cgodm.CLP2362, Cgodm.CLP2377, Cgodm....
    ## lgl  (1): diff_exprs

    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Data$ID <- as.factor(Data$ID)
Data$class <- as.factor(Data$class)
Data$toxin_family <- as.factor(Data$toxin_family)

#import list of genes filter by the coverage
genelist<-scan('Genes.txt', character(), quote = '')
#t<-read.table("Genes.txt", header = F)
#genelist<-c(t$V1)

Data_Filter<-data.frame()
for (i in 1:length(genelist)){
  Data_Filter[i,1:ncol(Data)]<-Data[which(genelist[i] == Data$ID), ]
}
Data_Original<-Data
Data<-Data_Filter
rm(Data_Filter)

Data<-Data[!duplicated(Data),]
```

``` r
source("tab_functions.R")
```

## rename toxins by average expression and toxin families

``` r
Data1<-Data_Original[order(as.vector(Data_Original$class),as.vector(Data_Original$average_tpm),decreasing = T),]
ToxNam<-data.frame(gene_id = NA,Tox_Name = NA)
for (i in 1:table(Data1$class)[2]){
  ToxNam[i,1]<-paste(Data1$ID[i])
  ToxNam[i,2]<-paste0(Data1$toxin_family[i],"_",i)
}

for (i in 1:nrow(ToxNam)){
  levels(Data$ID)[which(ToxNam[i,1] == levels(Data$ID))] <- ToxNam[i,2]
}


#Data %>% mutate_if(is.character,as.vector)
```

## snpEFF Data

``` r
##separate part of the data
snpEff<-Data[,c("ID","class","average_tpm","stdev_tpm","bcv_tpm","diff_exprs","length",
                "Total_Nonsynonymous","Synonymous","Total_variants")]



###filter data to drom missing values
snpEff<-snpEff[which(snpEff$Total_variants > 0),]
snpEff<-snpEff[complete.cases(snpEff), ]

#### add colums with snps per 1000 bp 
snpEff$total_snpKB<-(snpEff$Total_variants/snpEff$length)*1000
snpEff$syn_snpKB<-(snpEff$Synonymous/snpEff$length)*1000
snpEff$nonsyn_snpKB<-(snpEff$Total_Nonsynonymous/snpEff$length)*1000

###subset tables in Toxins and Nontoxins
snpEff_tox<-droplevels.data.frame(subset(snpEff,snpEff$class=="Toxin"))
snpEff_nontox<-droplevels.data.frame(subset(snpEff,snpEff$class=="Nontoxin"))

###summarize table by class and calculating number of genes by class mean smp per kb standar deviation and sum of snps and proprtion of synonymous and nonsynonymous
x<-snpEff %>% group_by(class) %>% summarize(n=n(),
                                            perKB_all=mean(total_snpKB),
                                            perKB_all_sd=sd(total_snpKB),
                                            Total_Nonsynonymous=sum(Total_Nonsynonymous),
                                            Total_Synonymous=sum(Synonymous),
                                            Total_Variants=sum(Total_variants),
                                            Prop_nonsyn=(Total_Nonsynonymous/Total_Variants),
                                            Prop_syn=(Total_Synonymous/Total_Variants))


kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
class
</th>
<th style="text-align:right;">
n
</th>
<th style="text-align:right;">
perKB\_all
</th>
<th style="text-align:right;">
perKB\_all\_sd
</th>
<th style="text-align:right;">
Total\_Nonsynonymous
</th>
<th style="text-align:right;">
Total\_Synonymous
</th>
<th style="text-align:right;">
Total\_Variants
</th>
<th style="text-align:right;">
Prop\_nonsyn
</th>
<th style="text-align:right;">
Prop\_syn
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Nontoxin
</td>
<td style="text-align:right;">
4508
</td>
<td style="text-align:right;">
6.765957
</td>
<td style="text-align:right;">
4.307966
</td>
<td style="text-align:right;">
8563
</td>
<td style="text-align:right;">
24094
</td>
<td style="text-align:right;">
32657
</td>
<td style="text-align:right;">
0.2622102
</td>
<td style="text-align:right;">
0.7377898
</td>
</tr>
<tr>
<td style="text-align:left;">
Toxin
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
6.423349
</td>
<td style="text-align:right;">
3.580440
</td>
<td style="text-align:right;">
236
</td>
<td style="text-align:right;">
154
</td>
<td style="text-align:right;">
390
</td>
<td style="text-align:right;">
0.6051282
</td>
<td style="text-align:right;">
0.3948718
</td>
</tr>
</tbody>
</table>

### Is there a significant association between the type of mutation (nonsynonymous vs. synonymous) and the transcript class (toxin vs. nontoxin)?

``` r
chisq.test(x[,5:6])
```

    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  x[, 5:6]
    ## X-squared = 230.22, df = 1, p-value < 2.2e-16

``` r
models_tab<-chi2_table(model = chisq.test(x[,5:6]),id = "Chi2 type of mutation ~ class")
```

### Is there a significant difference in the number of SNPs per kilobase between toxin and nontoxins?

``` r
##transform variable to improve normality of the data
hist(sqrt(snpEff$total_snpKB))
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
##linear model to compare snpKB for toxins vs nontoxins
(m<-summary(lm(sqrt(snpEff$total_snpKB)~snpEff$class)))
```

    ## 
    ## Call:
    ## lm(formula = sqrt(snpEff$total_snpKB) ~ snpEff$class)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0124 -0.6010 -0.0225  0.5570  3.3970 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        2.46742    0.01225 201.487   <2e-16 ***
    ## snpEff$classToxin -0.04208    0.10046  -0.419    0.675    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8222 on 4574 degrees of freedom
    ## Multiple R-squared:  3.836e-05,  Adjusted R-squared:  -0.0001803 
    ## F-statistic: 0.1754 on 1 and 4574 DF,  p-value: 0.6753

``` r
####to check that the number is not afecting the result Rhett did 1000 repetitions of the linear model taking randon subsets of the nontoxins of the size of the toxins
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="SNPs per Kilobase ~ class")
results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(snpEff_nontox$total_snpKB,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","total_snpKB")
  data<-snpEff_tox[,c("class","total_snpKB")]
  data<-rbind(data,samp)
  z<-summary(lm(sqrt(as.numeric(data$total_snpKB))~data$class))
  results[i,1]<-round(z$r.squared,2)
  results[i,2]<-round(z$coefficients[2,4],2)
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="SNP's per kb")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

### Is there a significant difference in the number of non synonymous SNPs per kilobase between toxin and nontoxins?

``` r
##transform variable to improve normality of the data
hist(as.vector(clr(snpEff$nonsyn_snpKB)))
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
##linear model to compare snpKB for toxins vs nontoxins
(m<-summary(lm(as.vector(clr(snpEff$nonsyn))~snpEff$class)))
```

    ## 
    ## Call:
    ## lm(formula = as.vector(clr(snpEff$nonsyn)) ~ snpEff$class)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.31975 -0.23470  0.00729  0.27744  2.04871 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       -0.007289   0.008859  -0.823    0.411    
    ## snpEff$classToxin  0.490502   0.072672   6.750 1.67e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5948 on 4574 degrees of freedom
    ## Multiple R-squared:  0.009862,   Adjusted R-squared:  0.009645 
    ## F-statistic: 45.56 on 1 and 4574 DF,  p-value: 1.668e-11

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="Nonsynonymous SNPs per Kilobase ~ class")
####to check that the number is not afecting the result Rhett did 1000 repetitions of the linear model taking randon subsets of the nontoxins of the size of the toxins
results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(snpEff$nonsyn_snpKB,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","nonsyn_snpKB")
  data<-snpEff_tox[,c("class","nonsyn_snpKB")]
  data<-rbind(data,samp)
  z<-summary(lm(as.vector(clr(as.numeric(data$nonsyn_snpKB)))~data$class))
  results[i,1]<-round(z$r.squared,2)
  results[i,2]<-round(z$coefficients[2,4],2)
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="nonsynonymous SNP's per kb")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

### Is there a significant difference in the number of synonymous SNPs per kilobase between toxin and nontoxins?

``` r
##transform variable to improve normality of the data
hist(sqrt(snpEff$syn_snpKB))
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
##linear model to compare snpKB for toxins vs nontoxins
(m<-summary(lm(sqrt(snpEff$syn_snpKB)~snpEff$class)))
```

    ## 
    ## Call:
    ## lm(formula = sqrt(snpEff$syn_snpKB) ~ snpEff$class)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0582 -0.5261  0.0002  0.5609  3.3363 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        2.05822    0.01283 160.480  < 2e-16 ***
    ## snpEff$classToxin -0.67558    0.10521  -6.421 1.49e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8611 on 4574 degrees of freedom
    ## Multiple R-squared:  0.008934,   Adjusted R-squared:  0.008717 
    ## F-statistic: 41.23 on 1 and 4574 DF,  p-value: 1.489e-10

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="Synonymous SNPs per Kilobase ~ class")
####to check that the number is not afecting the result Rhett did 1000 repetitions of the linear model taking randon subsets of the nontoxins of the size of the toxins
results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(snpEff_nontox$syn_snpKB,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","syn_snpKB")
  data<-snpEff_tox[,c("class","syn_snpKB")]
  data<-rbind(data,samp)
  z<-summary(lm(sqrt(as.numeric(data$syn_snpKB))~data$class))
  results[i,1]<-round(z$r.squared,2)
  results[i,2]<-round(z$coefficients[2,4],2)
}
table(results[,2] < 0.05)
```

    ## 
    ## TRUE 
    ## 1000

``` r
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("TRUE"))[2])
lbls<-c("True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("green"),main="synonymous SNP's per kb")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

## Nucleotide diversity

Read the data

``` r
NDiv<-Data[,c("ID","class","average_tpm","stdev_tpm","bcv_tpm","diff_exprs","pi","Total_variants","toxin_family")]
NDiv<-NDiv[which(NDiv$Total_variants > 0),]
NDiv<-NDiv[complete.cases(NDiv), ]
NDiv_tox<-droplevels.data.frame(subset(NDiv,NDiv$class=="Toxin"))
NDiv_nontox<-droplevels.data.frame(subset(NDiv,NDiv$class=="Nontoxin"))

x<-NDiv %>% group_by(class) %>% summarize(mean=mean(pi),sd=sd(pi),n=length(pi))

kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
class
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Nontoxin
</td>
<td style="text-align:right;">
0.3149234
</td>
<td style="text-align:right;">
0.0812495
</td>
<td style="text-align:right;">
4508
</td>
</tr>
<tr>
<td style="text-align:left;">
Toxin
</td>
<td style="text-align:right;">
0.3252790
</td>
<td style="text-align:right;">
0.0843843
</td>
<td style="text-align:right;">
68
</td>
</tr>
</tbody>
</table>

### Is there a significant difference in nucleotide diversity between toxins and nontoxins?

``` r
hist(NDiv$pi)
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
(m<-summary(lm(NDiv$pi~NDiv$class)))
```

    ## 
    ## Call:
    ## lm(formula = NDiv$pi ~ NDiv$class)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.158612 -0.048526 -0.005706  0.044384  0.268410 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     0.314923   0.001211 260.091   <2e-16 ***
    ## NDiv$classToxin 0.010356   0.009933   1.043    0.297    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0813 on 4574 degrees of freedom
    ## Multiple R-squared:  0.0002376,  Adjusted R-squared:  1.901e-05 
    ## F-statistic: 1.087 on 1 and 4574 DF,  p-value: 0.2972

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="Nucleotide diversity ~ class")
results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(NDiv_nontox$pi,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","pi")
  data<-NDiv_tox[,c("class","pi")]
  data<-rbind(data,samp)
  z<-summary(lm(data$pi~data$class))
  results[i,1]<-round(z$r.squared,2)
  results[i,2]<-round(z$coefficients[2,4],2)
}

pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="Nucleotide Diversity")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

Save the results for plotting and subset toxins which fall outside the
95th percentile of nontoxin estimates.

``` r
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
CI<-quantile(NDiv_nontox$pi,c(0.025,0.975))
NDiv_tox$sig<-NA
for (i in 1:nrow(NDiv_tox)){
  if (NDiv_tox$pi[i] < CI[1] | NDiv_tox$pi[i] > CI[2]){NDiv_tox$sig[i]<-TRUE}
  else {NDiv_tox$sig[i]<-FALSE}
}
sig<-droplevels.data.frame(subset(NDiv_tox,NDiv_tox$sig=="TRUE"))
print(sig$ID[which(sig$sig == T)])
```

    ## [1] VEGF_30 SVSP_46
    ## Levels: VEGF_30 SVSP_46

``` r
ggplot(NDiv, aes(x=class, y=pi)) + 
  geom_violin(aes(fill=class),trim=F) + 
  ylim(NA,(1.25*max(NDiv$pi))) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  geom_point(aes(x=class,y=pi),data=sig,fill="white",colour="black",size=3,pch=23) +
  geom_hline(yintercept=CI,linetype="longdash",color="darkgray")  +
  geom_label(data=sig,aes(x=2.25,label=ID,y=pi),label.size = 0,size=5,inherit.aes = F) +
  draw_label(paste("Nontoxins",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.2*max(NDiv$pi)),size=10) +
  draw_label(paste("Toxins",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.2*max(NDiv$pi)),size=10) +
  draw_label(paste("R2=",r2,"\n p =",p,"\n b =",b),x=1.5,y=(1*max(NDiv$pi)),size=10) +
  guides(fill= "none")+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### Does nucleotide diversity correlate with average expression levels?

``` r
x<-NDiv %>% group_by(diff_exprs) %>% summarize(mean=mean(pi),sd=sd(pi),n=length(pi))
kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
diff\_exprs
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
0.3150844
</td>
<td style="text-align:right;">
0.0811962
</td>
<td style="text-align:right;">
4522
</td>
</tr>
<tr>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
0.3144824
</td>
<td style="text-align:right;">
0.0901802
</td>
<td style="text-align:right;">
54
</td>
</tr>
</tbody>
</table>

``` r
(m<-summary(lm(NDiv$pi~NDiv$diff_exprs)))
```

    ## 
    ## Call:
    ## lm(formula = NDiv$pi ~ NDiv$diff_exprs)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.14842 -0.04842 -0.00556  0.04427  0.26825 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          0.3150844  0.0012091 260.597   <2e-16 ***
    ## NDiv$diff_exprsTRUE -0.0006019  0.0111302  -0.054    0.957    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08131 on 4574 degrees of freedom
    ## Multiple R-squared:  6.394e-07,  Adjusted R-squared:  -0.000218 
    ## F-statistic: 0.002925 on 1 and 4574 DF,  p-value: 0.9569

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="NDiv ~ Dif_Exp")
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
ggplot(NDiv, aes(x=diff_exprs, y=pi)) + 
  geom_violin(aes(fill=diff_exprs),trim=F) +
  xlab("Differentially Expressed") +
  ylim(NA,(1.25*max(NDiv$pi))) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  draw_label(paste("FALSE",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.2*max(NDiv$pi)),size=10) +
  draw_label(paste("TRUE",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.2*max(NDiv$pi)),size=10) +
  draw_label(paste("R2 =",r2,"\n p =",p),x=1.5,y=(1*max(NDiv$pi)),size=10) +
  guides(fill="none")+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
all_lm<-lm(NDiv$pi~log(NDiv$average_tpm))
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = NDiv$pi ~ log(NDiv$average_tpm))
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.148715 -0.048473 -0.005594  0.044438  0.268140 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            0.3153570  0.0027405 115.074   <2e-16 ***
    ## log(NDiv$average_tpm) -0.0001004  0.0008843  -0.114     0.91    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08131 on 4574 degrees of freedom
    ## Multiple R-squared:  2.821e-06,  Adjusted R-squared:  -0.0002158 
    ## F-statistic: 0.0129 on 1 and 4574 DF,  p-value: 0.9096

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="NDiv ~ expression")
r2<-round(m$r.squared,3)
p<-round(m$coefficients[2,4],3)
ggplot(NDiv,aes(x=log(average_tpm),y=pi)) + xlab("Average TPM") + ylab("Pi") +
  geom_point() + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
  draw_label(paste("n =",nrow(NDiv),"; R2=",r2,"; p =",p),x=7.5,y=0.55,size=10)+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

### Does nucleotide diversity correlate with average expression levels **in toxins**?

``` r
x<-NDiv_tox %>% group_by(diff_exprs) %>% summarize(mean=mean(pi),sd=sd(pi),n=length(pi))
kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
diff\_exprs
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
0.3267182
</td>
<td style="text-align:right;">
0.0840496
</td>
<td style="text-align:right;">
59
</td>
</tr>
<tr>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
0.3158437
</td>
<td style="text-align:right;">
0.0911268
</td>
<td style="text-align:right;">
9
</td>
</tr>
</tbody>
</table>

``` r
(m<-summary(lm(NDiv_tox$pi~NDiv_tox$diff_exprs)))
```

    ## 
    ## Call:
    ## lm(formula = NDiv_tox$pi ~ NDiv_tox$diff_exprs)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.16005 -0.05257 -0.01169  0.06722  0.21874 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              0.32672    0.01106  29.546   <2e-16 ***
    ## NDiv_tox$diff_exprsTRUE -0.01087    0.03040  -0.358    0.722    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08494 on 66 degrees of freedom
    ## Multiple R-squared:  0.001936,   Adjusted R-squared:  -0.01319 
    ## F-statistic: 0.128 on 1 and 66 DF,  p-value: 0.7217

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="Toxins nucleotide diversity ~ class")
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
ggplot(NDiv_tox, aes(x=diff_exprs, y=pi)) + 
  geom_violin(aes(fill=diff_exprs),trim=F) +
  xlab("Differentially Expressed") +
  ylim(NA,(1.25*max(NDiv$pi))) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  draw_label(paste("FALSE",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.2*max(NDiv$pi)),size=10) +
  draw_label(paste("TRUE",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.2*max(NDiv$pi)),size=10) +
  draw_label(paste("R2 =",r2,"\n p =",p),x=1.5,y=(1*max(NDiv$pi)),size=10) +
  guides(fill="none")+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
all_lm<-lm(NDiv_tox$pi~log(NDiv_tox$average_tpm))
AIC(all_lm)
```

    ## [1] -140.5587

``` r
models_tab<-lm_table(all_lm,smr=F,add=T,tab = models_tab,id="Toxins nucleotide diversity ~ expression",aic = T)
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = NDiv_tox$pi ~ log(NDiv_tox$average_tpm))
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.161785 -0.048371 -0.009289  0.048849  0.217880 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               0.230069   0.064216   3.583 0.000645 ***
    ## log(NDiv_tox$average_tpm) 0.012009   0.007998   1.501 0.137997    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08361 on 66 degrees of freedom
    ## Multiple R-squared:  0.03303,    Adjusted R-squared:  0.01838 
    ## F-statistic: 2.254 on 1 and 66 DF,  p-value: 0.138

``` r
r2<-round(m$r.squared,3)
p<-round(m$coefficients[2,4],3)
ggplot(NDiv_tox,aes(x=log(average_tpm),y=pi)) + xlab("Average TPM") + ylab("Pi") +
  geom_point() + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
  draw_label(paste("n =",nrow(NDiv_tox),"; R2=",r2,"; p =",p),x=7.5,y=0.4,size=10)+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

### Does nucleotide diversity correlate with average expression levels and toxin family **in toxins**?

``` r
all_lm<-lm(NDiv_tox$pi~log(NDiv_tox$average_tpm)+NDiv_tox$toxin_family)
models_tab<-lm_table(all_lm,tf=T,smr=F,add=T,tab = models_tab,id="Toxins nucleotide diversity ~ expression and toxin family",aic = T)
AIC(all_lm)
```

    ## [1] -126.0214

``` r
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = NDiv_tox$pi ~ log(NDiv_tox$average_tpm) + NDiv_tox$toxin_family)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.17249 -0.03778  0.00000  0.04301  0.18935 
    ## 
    ## Coefficients:
    ##                                       Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                           0.196737   0.115315   1.706   0.0941 .
    ## log(NDiv_tox$average_tpm)             0.014598   0.011351   1.286   0.2042  
    ## NDiv_tox$toxin_familyCRISP           -0.044300   0.098065  -0.452   0.6534  
    ## NDiv_tox$toxin_familyCTL             -0.003779   0.059044  -0.064   0.9492  
    ## NDiv_tox$toxin_familyHYAL             0.036796   0.107028   0.344   0.7324  
    ## NDiv_tox$toxin_familyKUN             -0.059207   0.107817  -0.549   0.5853  
    ## NDiv_tox$toxin_familyLAAO            -0.077698   0.098043  -0.792   0.4317  
    ## NDiv_tox$toxin_familyNUC              0.018734   0.101533   0.185   0.8543  
    ## NDiv_tox$toxin_familyPDE              0.059338   0.104975   0.565   0.5744  
    ## NDiv_tox$toxin_familyPLA2            -0.022347   0.059930  -0.373   0.7108  
    ## NDiv_tox$toxin_familyPLB              0.052422   0.082356   0.637   0.5273  
    ## NDiv_tox$toxin_familySVMPI           -0.038287   0.069344  -0.552   0.5833  
    ## NDiv_tox$toxin_familySVMPII          -0.005924   0.060058  -0.099   0.9218  
    ## NDiv_tox$toxin_familySVMPIII          0.030325   0.056437   0.537   0.5934  
    ## NDiv_tox$toxin_familySVSP             0.040830   0.055698   0.733   0.4669  
    ## NDiv_tox$toxin_familyuncharacterised  0.036414   0.101534   0.359   0.7213  
    ## NDiv_tox$toxin_familyVEGF             0.194917   0.098119   1.987   0.0524 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08489 on 51 degrees of freedom
    ## Multiple R-squared:  0.2297, Adjusted R-squared:  -0.01196 
    ## F-statistic: 0.9505 on 16 and 51 DF,  p-value: 0.521

``` r
r2<-round(m$r.squared,3)
p<-round(m$coefficients[2,4],3)
ggplot(NDiv_tox,aes(x=log(average_tpm),y=pi)) + xlab("Average TPM") + ylab("Pi") +
  geom_point() + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
  draw_label(paste("n =",nrow(NDiv_tox),"; R2=",r2,"; p =",p),x=7.5,y=0.4,size=10)+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## Tajima’s D

get the data

``` r
Tajima<-Data[,c("ID","class","average_tpm","stdev_tpm","bcv_tpm","diff_exprs","tajimasD","Total_variants","toxin_family")]
Tajima<-Tajima[which(Tajima$Total_variants > 0),]
Tajima<-Tajima[complete.cases(Tajima), ]

Tajima_tox<-droplevels.data.frame(subset(Tajima,Tajima$class=="Toxin"))
Tajima_nontox<-droplevels.data.frame(subset(Tajima,Tajima$class=="Nontoxin"))

x<-Tajima %>% group_by(class) %>% summarize(mean=mean(tajimasD),sd=sd(tajimasD),n=length(tajimasD))

kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
class
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Nontoxin
</td>
<td style="text-align:right;">
-0.1875167
</td>
<td style="text-align:right;">
0.7775041
</td>
<td style="text-align:right;">
4508
</td>
</tr>
<tr>
<td style="text-align:left;">
Toxin
</td>
<td style="text-align:right;">
-0.0346594
</td>
<td style="text-align:right;">
0.8240109
</td>
<td style="text-align:right;">
68
</td>
</tr>
</tbody>
</table>

### Are estimates of Tajima’s D significantly different than 0 (suggesting selection) for toxins or nontoxins?

``` r
t.test(Tajima_nontox$tajimasD)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  Tajima_nontox$tajimasD
    ## t = -16.193, df = 4507, p-value < 2.2e-16
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2102193 -0.1648141
    ## sample estimates:
    ##  mean of x 
    ## -0.1875167

``` r
sd(Tajima_nontox$tajimasD)
```

    ## [1] 0.7775041

``` r
models_tab <-t_table(model = t.test(Tajima_nontox$tajimasD), add=T,tab = models_tab,id="NonToxins TD t.test")
t.test(Tajima_tox$tajimasD)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  Tajima_tox$tajimasD
    ## t = -0.34685, df = 67, p-value = 0.7298
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2341125  0.1647937
    ## sample estimates:
    ##  mean of x 
    ## -0.0346594

``` r
sd(Tajima_tox$tajimasD)
```

    ## [1] 0.8240109

``` r
models_tab <-t_table(model = t.test(Tajima_tox$tajimasD), add=T,tab = models_tab,id="Toxins TD t.test")
t.test(Tajima_nontox$tajimasD,Tajima_tox$tajimasD)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  Tajima_nontox$tajimasD and Tajima_tox$tajimasD
    ## t = -1.5195, df = 68.811, p-value = 0.1332
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.35354819  0.04783351
    ## sample estimates:
    ##  mean of x  mean of y 
    ## -0.1875167 -0.0346594

``` r
models_tab <-t_table(model = t.test(Tajima_nontox$tajimasD,Tajima_tox$tajimasD), add=T,tab = models_tab,id="Toxins vs Nontoxins TD t.test")
var.test(Tajima$tajimasD ~ Tajima$class)
```

    ## 
    ##  F test to compare two variances
    ## 
    ## data:  Tajima$tajimasD by Tajima$class
    ## F = 0.89031, num df = 4507, denom df = 67, p-value = 0.4617
    ## alternative hypothesis: true ratio of variances is not equal to 1
    ## 95 percent confidence interval:
    ##  0.6134523 1.2196957
    ## sample estimates:
    ## ratio of variances 
    ##          0.8903061

``` r
models_tab <-f_table(model = var.test(Tajima$tajimasD ~ Tajima$class), add=T,tab = models_tab,id="Toxins vs Nontoxins TD F.test")
```

``` r
results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(Tajima_nontox$tajimasD,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","tajimasD")
  data<-samp
  z<-t.test(as.numeric(data[,2]))
  results[i,1]<-round(z$estimate,2)
  results[i,2]<-round(z$p.value,2)
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="Tajima's D")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

### Is there a significant difference in Tajima’s D between toxins and nontoxins?

``` r
hist(Tajima$tajimasD)
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
(m<-summary(lm(Tajima$tajimasD~Tajima$class)))
```

    ## 
    ## Call:
    ## lm(formula = Tajima$tajimasD ~ Tajima$class)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7959 -0.5776 -0.0089  0.5094  2.3028 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       -0.18752    0.01159 -16.178   <2e-16 ***
    ## Tajima$classToxin  0.15286    0.09508   1.608    0.108    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7782 on 4574 degrees of freedom
    ## Multiple R-squared:  0.0005647,  Adjusted R-squared:  0.0003462 
    ## F-statistic: 2.585 on 1 and 4574 DF,  p-value: 0.108

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="TD ~ class")

results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(Tajima_nontox$tajimasD,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","tajimasD")
  data<-Tajima_tox[,c("class","tajimasD")]
  data<-rbind(data,samp)
  z<-summary(lm(data$tajimasD~data$class))
  results[i,1]<-round(z$r.squared,2)
  results[i,2]<-round(z$coefficients[2,4],2)
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="Tajima's D")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

``` r
#Save the results for plotting and subset toxins which fall outside the 95th percentile of nontoxin estimates.
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
CI<-quantile(Tajima_nontox$tajimasD,c(0.025,0.975))
Tajima_tox$sig<-NA
for (i in 1:nrow(Tajima_tox)){
  if (Tajima_tox$tajimasD[i] <= CI[1] | Tajima_tox$tajimasD[i] >= CI[2]){Tajima_tox$sig[i]<-TRUE}
  else {Tajima_tox$sig[i]<-FALSE}
}
sig<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$sig=="TRUE"))
print(sig[which(sig$sig == T),c("ID","tajimasD")])
```

    ##           ID tajimasD
    ## 819  VEGF_30  1.87288
    ## 829  SVSP_52  1.44125
    ## 851  SVSP_46  1.48617
    ## 873 SVMPI_44 -1.62929

### Figure 2A

``` r
(Fig2A <- ggplot(Tajima, aes(x=class, y=tajimasD)) + 
    geom_violin(aes(fill=class),trim=F) + 
    ylim(NA,(1.5*max(Tajima$tajimasD))) +
    #scale_fill_manual(values=c("#00BFC4","#F8766D")) +
    scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
    stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
    geom_point(aes(x=class,y=tajimasD),data=sig,fill="white",colour="black",size=3,pch=23) + #,position=position_jitter(width=0.1,seed=1)) +
    geom_hline(yintercept=CI,linetype="longdash",color="darkgray")  +
    geom_label(data=sig,aes(x=2.25,label=ID,y=tajimasD),label.size = 0,size=4,inherit.aes = F) + #, position=position_jitter(width=0,height=0.4,seed=1)) +
    draw_label(paste("Nontoxins",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.45*max(Tajima$tajimasD)),size=10) +
    draw_label(paste("Toxins",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.45*max(Tajima$tajimasD)),size=10) +
    draw_label(paste("R2 =",r2,"\n p =",p,"\n b =",b),x=1.5,y=(1.25*max(Tajima$tajimasD)),size=10) +
    guides(fill="none"))+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

### Does Tajima’s D (All, Positive, and Negative values) correlate with average expression levels?

``` r
x<-Tajima %>% group_by(diff_exprs) %>% summarize(mean=mean(tajimasD),sd=sd(tajimasD),n=length(tajimasD))
kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
diff\_exprs
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
-0.1858890
</td>
<td style="text-align:right;">
0.7786183
</td>
<td style="text-align:right;">
4522
</td>
</tr>
<tr>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
-0.1313354
</td>
<td style="text-align:right;">
0.7598011
</td>
<td style="text-align:right;">
54
</td>
</tr>
</tbody>
</table>

``` r
(m<-summary(lm(Tajima$tajimasD~Tajima$diff_exprs)))
```

    ## 
    ## Call:
    ## lm(formula = Tajima$tajimasD ~ Tajima$diff_exprs)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.79754 -0.57450 -0.00903  0.50781  2.30121 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)           -0.18589    0.01158 -16.059   <2e-16 ***
    ## Tajima$diff_exprsTRUE  0.05455    0.10656   0.512    0.609    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7784 on 4574 degrees of freedom
    ## Multiple R-squared:  5.73e-05,   Adjusted R-squared:  -0.0001613 
    ## F-statistic: 0.2621 on 1 and 4574 DF,  p-value: 0.6087

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="TD ~ Dif_Exp")
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)

library(Hmisc)
```

    ## Loading required package: lattice

    ## Loading required package: survival

    ## Loading required package: Formula

    ## 
    ## Attaching package: 'Hmisc'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     src, summarize

    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, units

``` r
ggplot(Tajima, aes(x=diff_exprs, y=tajimasD,trim=F)) + 
  geom_violin(aes(fill=diff_exprs)) +
  xlab("Differentially Expressed") +
  ylim(NA,(1.5*max(Tajima$tajimasD))) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  draw_label(paste("FALSE",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.45*max(Tajima$tajimasD)),size=10) +
  draw_label(paste("TRUE",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.45*max(Tajima$tajimasD)),size=10) +
  draw_label(paste("R2 =",r2,"\n p =",p),x=1.5,y=(1.25*max(Tajima$tajimasD)),size=10) +
  guides(fill="none")+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
detach('package:Hmisc',unload = T)


all_lm<-lm(Tajima$tajimasD~log(Tajima$average_tpm))
models_tab<-lm_table(all_lm,smr=F,add=T,tab = models_tab,id="TD all ~ expression",aic = T)
AIC(all_lm)
```

    ## [1] 10697.68

``` r
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = Tajima$tajimasD ~ log(Tajima$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.79954 -0.57623 -0.01283  0.50661  2.29979 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             -0.181273   0.026237  -6.909 5.55e-12 ***
    ## log(Tajima$average_tpm) -0.001426   0.008466  -0.168    0.866    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7784 on 4574 degrees of freedom
    ## Multiple R-squared:  6.203e-06,  Adjusted R-squared:  -0.0002124 
    ## F-statistic: 0.02837 on 1 and 4574 DF,  p-value: 0.8662

``` r
lt0<-droplevels.data.frame(subset(Tajima,Tajima$tajimasD<0))
lt0_lm<-lm(lt0$tajimasD~log(lt0$average_tpm))
AIC(lt0_lm)
```

    ## [1] 3676.821

``` r
models_tab<-lm_table(lt0_lm,smr=F,add=T,tab = models_tab,id="TD lower than 0 ~ expression",aic = T)
(m_lt0<-summary(lt0_lm))
```

    ## 
    ## Call:
    ## lm(formula = lt0$tajimasD ~ log(lt0$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.31153 -0.39705  0.07924  0.41477  0.69887 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          -0.657081   0.020617 -31.871   <2e-16 ***
    ## log(lt0$average_tpm) -0.008084   0.006716  -1.204    0.229    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4645 on 2814 degrees of freedom
    ## Multiple R-squared:  0.0005145,  Adjusted R-squared:  0.0001594 
    ## F-statistic: 1.449 on 1 and 2814 DF,  p-value: 0.2288

``` r
gt0<-droplevels.data.frame(subset(Tajima,Tajima$tajimasD>0))
gt0_lm<-lm(gt0$tajimasD~log(gt0$average_tpm))
AIC(gt0_lm)
```

    ## [1] 2280.441

``` r
models_tab<-lm_table(gt0_lm,smr=F,add=T,tab = models_tab,id="TD greater than 0 ~ expression",aic = T)
(m_gt0<-summary(gt0_lm))
```

    ## 
    ## Call:
    ## lm(formula = gt0$tajimasD ~ log(gt0$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.61189 -0.38452 -0.06687  0.28689  1.51050 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          0.601502   0.023986  25.077   <2e-16 ***
    ## log(gt0$average_tpm) 0.001481   0.007625   0.194    0.846    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.462 on 1758 degrees of freedom
    ## Multiple R-squared:  2.145e-05,  Adjusted R-squared:  -0.0005474 
    ## F-statistic: 0.0377 on 1 and 1758 DF,  p-value: 0.8461

``` r
ggplot(Tajima,aes(x=log(average_tpm),y=tajimasD)) + xlab("Average TPM") + ylab("Tajima's D") +
  geom_point() + geom_hline(yintercept=0) + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
  geom_abline(intercept = coef(lt0_lm)[1],slope=coef(lt0_lm)[2],lty=3) +
  geom_abline(intercept = coef(gt0_lm)[1],slope=coef(gt0_lm)[2],lty=2) +
  draw_label(paste("n =",nrow(Tajima),"; R2=",round(m$r.squared,2),"; p =",round(m$coefficients[2,4],2)),x=7.5,y=0,size=10) +
  draw_label(paste("n =",nrow(lt0),"; R2=",round(m_lt0$r.squared,2),"; p =",round(m_lt0$coefficients[2,4],2)),x=7.5,y=-1,size=10) +
  draw_label(paste("n =",nrow(gt0),"; R2=",round(m_gt0$r.squared,2),"; p =",round(m_gt0$coefficients[2,4],2)),x=7.5,y=1,size=10)+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

### Does Tajima’s D (All, Positive, and Negative values) correlate with average expression levels **in toxins**?

``` r
x<-Tajima_tox %>% group_by(diff_exprs) %>% summarize(mean=mean(tajimasD),sd=sd(tajimasD),n=length(tajimasD))
kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
diff\_exprs
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
-0.0389956
</td>
<td style="text-align:right;">
0.8432139
</td>
<td style="text-align:right;">
59
</td>
</tr>
<tr>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
-0.0062333
</td>
<td style="text-align:right;">
0.7284969
</td>
<td style="text-align:right;">
9
</td>
</tr>
</tbody>
</table>

``` r
(m<-summary(lm(Tajima_tox$tajimasD~Tajima_tox$diff_exprs)))
```

    ## 
    ## Call:
    ## lm(formula = Tajima_tox$tajimasD ~ Tajima_tox$diff_exprs)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.5903 -0.5601 -0.1476  0.6231  1.9119 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)               -0.03900    0.10808  -0.361    0.719
    ## Tajima_tox$diff_exprsTRUE  0.03276    0.29707   0.110    0.913
    ## 
    ## Residual standard error: 0.8302 on 66 degrees of freedom
    ## Multiple R-squared:  0.0001842,  Adjusted R-squared:  -0.01496 
    ## F-statistic: 0.01216 on 1 and 66 DF,  p-value: 0.9125

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="Toxins TD ~ Dif_Exp")
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
ggplot(Tajima_tox, aes(x=diff_exprs, y=tajimasD)) + 
  geom_violin(aes(fill=diff_exprs),trim=F) +
  xlab("Differentially Expressed") +
  ylim(NA,(1.5*max(Tajima_tox$tajimasD))) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  draw_label(paste("FALSE",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.45*max(Tajima_tox$tajimasD)),size=10) +
  draw_label(paste("TRUE",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.45*max(Tajima_tox$tajimasD)),size=10) +
  draw_label(paste("R2=",r2,"; p =",p),x=1.5,y=(1.25*max(Tajima_tox$tajimasD)),size=10) +
  guides(fill="none")+
  theme_classic()
```

    ## Warning: Removed 7 rows containing missing values (geom_violin).

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
all_lm<-lm(Tajima_tox$tajimasD~log(Tajima_tox$average_tpm))
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = Tajima_tox$tajimasD ~ log(Tajima_tox$average_tpm))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6281 -0.5553 -0.1395  0.6137  1.7939 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                  -1.0617     0.6247  -1.699   0.0939 .
    ## log(Tajima_tox$average_tpm)   0.1295     0.0778   1.665   0.1007  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8133 on 66 degrees of freedom
    ## Multiple R-squared:  0.0403, Adjusted R-squared:  0.02576 
    ## F-statistic: 2.772 on 1 and 66 DF,  p-value: 0.1007

``` r
models_tab<-lm_table(all_lm,smr=F,add=T,tab = models_tab,id="TOXINS TD all ~ expression",aic = T)
lt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD<0))
lt0_lm<-lm(lt0$tajimasD~log(lt0$average_tpm))
models_tab<-lm_table(lt0_lm,smr=F,add=T,tab = models_tab,id="TOXINS TD lower than 0 ~ expression",aic = T)
(m_lt0<-summary(lt0_lm))
```

    ## 
    ## Call:
    ## lm(formula = lt0$tajimasD ~ log(lt0$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.02890 -0.31715  0.08377  0.34806  0.74393 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)          -1.35321    0.43161  -3.135  0.00341 **
    ## log(lt0$average_tpm)  0.09196    0.05491   1.675  0.10268   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4357 on 36 degrees of freedom
    ## Multiple R-squared:  0.07227,    Adjusted R-squared:  0.0465 
    ## F-statistic: 2.804 on 1 and 36 DF,  p-value: 0.1027

``` r
gt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD>0))
gt0_lm<-lm(gt0$tajimasD~log(gt0$average_tpm))
models_tab<-lm_table(gt0_lm,smr=F,add=T,tab = models_tab,id="TOXINS TD greater than 0 ~ expression",aic = T)
(m_gt0<-summary(gt0_lm))
```

    ## 
    ## Call:
    ## lm(formula = gt0$tajimasD ~ log(gt0$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.73258 -0.39156  0.03817  0.35610  1.15172 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)           0.87139    0.60192   1.448    0.159
    ## log(gt0$average_tpm) -0.01706    0.07306  -0.234    0.817
    ## 
    ## Residual standard error: 0.483 on 28 degrees of freedom
    ## Multiple R-squared:  0.001944,   Adjusted R-squared:  -0.0337 
    ## F-statistic: 0.05453 on 1 and 28 DF,  p-value: 0.8171

``` r
print("all values")
```

    ## [1] "all values"

``` r
AIC(all_lm)
```

    ## [1] 168.8451

``` r
print("lower than 0 values")
```

    ## [1] "lower than 0 values"

``` r
AIC(lt0_lm)
```

    ## [1] 48.64051

``` r
print("greater than 0 values")
```

    ## [1] "greater than 0 values"

``` r
AIC(gt0_lm)
```

    ## [1] 45.39829

### Does Tajima’s D (All, Positive, and Negative values) correlate with average expression levels and toxin families **in toxins**?

``` r
all_lm1<-lm(Tajima_tox$tajimasD~log(Tajima_tox$average_tpm)+Tajima_tox$toxin_family)
models_tab<-lm_table(all_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS TD all ~ expression + toxin family",aic = T)
(m<-summary(all_lm1))
```

    ## 
    ## Call:
    ## lm(formula = Tajima_tox$tajimasD ~ log(Tajima_tox$average_tpm) + 
    ##     Tajima_tox$toxin_family)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.43671 -0.42270 -0.01501  0.43482  1.48185 
    ## 
    ## Coefficients:
    ##                                        Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                            -1.45056    1.09678  -1.323    0.192  
    ## log(Tajima_tox$average_tpm)             0.15637    0.10796   1.448    0.154  
    ## Tajima_tox$toxin_familyCRISP           -0.50254    0.93271  -0.539    0.592  
    ## Tajima_tox$toxin_familyCTL              0.06895    0.56158   0.123    0.903  
    ## Tajima_tox$toxin_familyHYAL             0.37992    1.01796   0.373    0.711  
    ## Tajima_tox$toxin_familyKUN             -0.78116    1.02546  -0.762    0.450  
    ## Tajima_tox$toxin_familyLAAO            -0.94994    0.93250  -1.019    0.313  
    ## Tajima_tox$toxin_familyNUC              0.18100    0.96570   0.187    0.852  
    ## Tajima_tox$toxin_familyPDE              0.67469    0.99843   0.676    0.502  
    ## Tajima_tox$toxin_familyPLA2            -0.19248    0.57000  -0.338    0.737  
    ## Tajima_tox$toxin_familyPLB              0.49295    0.78330   0.629    0.532  
    ## Tajima_tox$toxin_familySVMPI           -0.34391    0.65954  -0.521    0.604  
    ## Tajima_tox$toxin_familySVMPII          -0.05988    0.57122  -0.105    0.917  
    ## Tajima_tox$toxin_familySVMPIII          0.44673    0.53678   0.832    0.409  
    ## Tajima_tox$toxin_familySVSP             0.43859    0.52975   0.828    0.412  
    ## Tajima_tox$toxin_familyuncharacterised  0.39918    0.96570   0.413    0.681  
    ## Tajima_tox$toxin_familyVEGF             1.94648    0.93323   2.086    0.042 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8074 on 51 degrees of freedom
    ## Multiple R-squared:  0.2692, Adjusted R-squared:  0.03997 
    ## F-statistic: 1.174 on 16 and 51 DF,  p-value: 0.3194

``` r
lt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD<0))
lt0_lm1<-lm(lt0$tajimasD~log(lt0$average_tpm)+lt0$toxin_family)
models_tab<-lm_table(lt0_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS TD lower than 0 ~ expression+ toxin family",aic = T)
(m_lt0<-summary(lt0_lm1))
```

    ## 
    ## Call:
    ## lm(formula = lt0$tajimasD ~ log(lt0$average_tpm) + lt0$toxin_family)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.82059 -0.05389  0.00000  0.12640  0.71452 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -1.42903    0.68889  -2.074   0.0489 *
    ## log(lt0$average_tpm)     0.13828    0.07158   1.932   0.0653 .
    ## lt0$toxin_familyCRISP   -0.36251    0.49237  -0.736   0.4687  
    ## lt0$toxin_familyCTL     -0.18746    0.32913  -0.570   0.5743  
    ## lt0$toxin_familyHYAL     0.45623    0.54777   0.833   0.4131  
    ## lt0$toxin_familyKUN     -0.70792    0.55322  -1.280   0.2129  
    ## lt0$toxin_familyLAAO    -0.80169    0.49422  -1.622   0.1178  
    ## lt0$toxin_familyNUC      0.28360    0.51075   0.555   0.5839  
    ## lt0$toxin_familyPLA2    -0.46281    0.33951  -1.363   0.1855  
    ## lt0$toxin_familyPLB      0.16767    0.51817   0.324   0.7491  
    ## lt0$toxin_familySVMPI   -1.33234    0.49398  -2.697   0.0126 *
    ## lt0$toxin_familySVMPII  -0.21755    0.34918  -0.623   0.5391  
    ## lt0$toxin_familySVMPIII -0.33616    0.41197  -0.816   0.4225  
    ## lt0$toxin_familySVSP    -0.34252    0.32669  -1.048   0.3049  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4019 on 24 degrees of freedom
    ## Multiple R-squared:  0.4737, Adjusted R-squared:  0.1886 
    ## F-statistic: 1.662 on 13 and 24 DF,  p-value: 0.1362

``` r
gt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD>0))
gt0_lm1<-lm(gt0$tajimasD~log(gt0$average_tpm)+gt0$toxin_family)
models_tab<-lm_table(gt0_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS TD greater than 0 ~ expression+ toxin family",aic = T)
(m_gt0<-summary(gt0_lm1))
```

    ## 
    ## Call:
    ## lm(formula = gt0$tajimasD ~ log(gt0$average_tpm) + gt0$toxin_family)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.4532 -0.2144  0.0000  0.1855  0.5669 
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)                      0.18578    0.85971   0.216  0.83134   
    ## log(gt0$average_tpm)             0.02086    0.07852   0.266  0.79351   
    ## gt0$toxin_familyCTL              0.57700    0.46846   1.232  0.23392   
    ## gt0$toxin_familyPDE             -0.16412    0.58226  -0.282  0.78125   
    ## gt0$toxin_familyPLA2             0.24438    0.44708   0.547  0.59137   
    ## gt0$toxin_familyPLB              0.22293    0.53891   0.414  0.68400   
    ## gt0$toxin_familySVMPI           -0.14718    0.42009  -0.350  0.73014   
    ## gt0$toxin_familySVMPII           0.09155    0.42964   0.213  0.83366   
    ## gt0$toxin_familySVMPIII          0.21998    0.37547   0.586  0.56523   
    ## gt0$toxin_familySVSP             0.77167    0.39881   1.935  0.06887 . 
    ## gt0$toxin_familyuncharacterised -0.30721    0.54285  -0.566  0.57844   
    ## gt0$toxin_familyVEGF             1.50340    0.49145   3.059  0.00676 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3405 on 18 degrees of freedom
    ## Multiple R-squared:  0.681,  Adjusted R-squared:  0.4861 
    ## F-statistic: 3.494 on 11 and 18 DF,  p-value: 0.00919

``` r
print("all values")
```

    ## [1] "all values"

``` r
AIC(all_lm1)
```

    ## [1] 180.3137

``` r
print("lower than 0 values")
```

    ## [1] "lower than 0 values"

``` r
AIC(lt0_lm1)
```

    ## [1] 51.10008

``` r
print("greater than 0 values")
```

    ## [1] "greater than 0 values"

``` r
AIC(gt0_lm1)
```

    ## [1] 31.17838

### Figure 3A

``` r
(testplot<-ggplot(Tajima_tox,aes(x=log(average_tpm),y=tajimasD, color = toxin_family)) + xlab("Average TPM") + ylab("Tajima's D") +
    geom_point() + #geom_hline(yintercept=0) + 
    geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
    geom_abline(intercept = coef(lt0_lm)[1],slope=coef(lt0_lm)[2],lty=3) +
    geom_abline(intercept = coef(gt0_lm)[1],slope=coef(gt0_lm)[2],lty=2) +
    draw_label(paste("n =",nrow(Tajima_tox),"; R2=",round(m$r.squared,2),"; p =",round(m$coefficients[2,4],2)),x=7.5,y=0,size=10) +
    draw_label(paste("n =",nrow(lt0),"; R2=",round(m_lt0$r.squared,2),"; p =",round(m_lt0$coefficients[2,4],2)),x=7.5,y=-1,size=10) +
    draw_label(paste("n =",nrow(gt0),"; R2=",round(m_gt0$r.squared,2),"; p =",round(m_gt0$coefficients[2,4],2)),x=7.5,y=1,size=10))+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
(Fig3A<-ggplot(Tajima_tox,aes(x=log(average_tpm),y=tajimasD)) + xlab("Average TPM") + ylab("Tajima's D") +
    geom_point() + #geom_hline(yintercept=0) + 
    geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
    geom_abline(intercept = coef(lt0_lm)[1],slope=coef(lt0_lm)[2],lty=3) +
    geom_abline(intercept = coef(gt0_lm)[1],slope=coef(gt0_lm)[2],lty=2) +
    draw_label(paste("n =",nrow(Tajima_tox),"; R2=",round(m$r.squared,2),"; p =",round(m$coefficients[2,4],2)),x=7.5,y=0,size=10) +
    draw_label(paste("n =",nrow(lt0),"; R2=",round(m_lt0$r.squared,2),"; p =",round(m_lt0$coefficients[2,4],2)),x=7.5,y=-1,size=10) +
    draw_label(paste("n =",nrow(gt0),"; R2=",round(m_gt0$r.squared,2),"; p =",round(m_gt0$coefficients[2,4],2)),x=7.5,y=1,size=10))+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

## Global Fst

``` r
Global<-Data[,c("ID","class","average_tpm","stdev_tpm","bcv_tpm","diff_exprs","SPvNP",'Total_variants',"toxin_family")]

ind<- colnames(Global) == "SPvNP"
colnames(Global)[ind]<- "Global"

Global<-Global[which(Global$Total_variants > 0),]
Global<-Global[complete.cases(Global), ]
Global_tox<-droplevels.data.frame(subset(Global,class=="Toxin"))
Global_nontox<-droplevels.data.frame(subset(Global,class=="Nontoxin"))

x<-Global %>% group_by(class) %>% summarize(mean=mean(Global),sd=sd(Global),n=length(Global))

kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
class
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Nontoxin
</td>
<td style="text-align:right;">
0.1894019
</td>
<td style="text-align:right;">
0.2087274
</td>
<td style="text-align:right;">
4508
</td>
</tr>
<tr>
<td style="text-align:left;">
Toxin
</td>
<td style="text-align:right;">
0.2028064
</td>
<td style="text-align:right;">
0.2552510
</td>
<td style="text-align:right;">
68
</td>
</tr>
</tbody>
</table>

### Is there a significant difference in Global Fst between toxins and nontoxins?

``` r
(m<-summary(lm(Global$Global~Global$class)))
```

    ## 
    ## Call:
    ## lm(formula = Global$Global ~ Global$class)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.67906 -0.12697 -0.00232  0.11213  0.81060 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        0.18940    0.00312  60.705   <2e-16 ***
    ## Global$classToxin  0.01341    0.02559   0.524      0.6    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2095 on 4574 degrees of freedom
    ## Multiple R-squared:  5.996e-05,  Adjusted R-squared:  -0.0001587 
    ## F-statistic: 0.2743 on 1 and 4574 DF,  p-value: 0.6005

``` r
models_tab<-lm_table(m,smr=T,add=T,tab=models_tab,id="Global Fst ~ class")

results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(Global_nontox$Global,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","Global")
  data<-Global_tox[,c("class","Global")]
  data<-rbind(data,samp)
  z<-summary(lm(data$Global~data$class))
  results[i,1]<-round(z$r.squared,2)
  results[i,2]<-round(z$coefficients[2,4],2)
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="Fst")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

``` r
###Save the results for plotting and subset toxins which fall outside the 95th percentile of nontoxin estimates.
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
CI<-quantile(Global_nontox$Global,c(0.95))
Global_tox$sig<-NA
for (i in 1:nrow(Global_tox)){
  if (Global_tox$Global[i] >= CI){Global_tox$sig[i]<-TRUE}
  else {Global_tox$sig[i]<-FALSE}
}
sig<-droplevels.data.frame(subset(Global_tox,Global_tox$sig=="TRUE"))
sig[which(sig$sig == T), c("ID","Global")]
```

    ##             ID    Global
    ## 819    VEGF_30 0.8145293
    ## 829    SVSP_52 0.6532508
    ## 833    SVSP_72 0.6013612
    ## 838    SVSP_64 0.6588273
    ## 841 SVMPIII_39 0.5940595
    ## 852 SVMPIII_55 0.5527096
    ## 886 SVMPIII_22 0.6174602

### Figure 2B

``` r
(Fig2B <- ggplot(Global, aes(x=class, y=Global)) + 
    geom_violin(aes(fill=class),trim=F) + 
    ylim(NA,(1.25*max(Global$Global))) +
    #scale_fill_manual(values=c("#00BFC4","#F8766D")) +
    scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
    stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
    geom_point(aes(x=class,y=Global),data=sig,fill="white",colour="black",size=3,pch=23) + #,position=position_jitter(width=0.1,seed=1)) +
    geom_hline(yintercept=CI,linetype="longdash",color="darkgray")  +
    geom_label(data=sig,aes(x=2.25,label=ID,y=Global),label.size = 0,size=4,inherit.aes = F) + #, position=position_jitter(width=0,height=0.4,seed=1)) +
    draw_label(paste("Nontoxins",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.25*max(Global$Global)),size=10) +
    draw_label(paste("Toxins",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.25*max(Global$Global)),size=10) +
    draw_label(paste("R2 =",r2,"\n p =",p,"\n b =",b),x=1.5,y=(max(Global$Global)),size=10) +
    guides(fill="none"))+
  ylab("Fst")+
  xlab("Class")+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

### Does Global Fst correlate with average expression levels?

``` r
x<-Global %>% group_by(diff_exprs) %>% summarize(mean=mean(Global),sd=sd(Global),n=length(Global))
kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
diff\_exprs
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
0.1906812
</td>
<td style="text-align:right;">
0.2092294
</td>
<td style="text-align:right;">
4522
</td>
</tr>
<tr>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
0.0991572
</td>
<td style="text-align:right;">
0.2115429
</td>
<td style="text-align:right;">
54
</td>
</tr>
</tbody>
</table>

``` r
(m<-summary(lm(Global$Global~Global$diff_exprs)))
```

    ## 
    ## Call:
    ## lm(formula = Global$Global ~ Global$diff_exprs)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.68034 -0.12674 -0.00307  0.11211  0.80932 
    ## 
    ## Coefficients:
    ##                        Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            0.190681   0.003112  61.277  < 2e-16 ***
    ## Global$diff_exprsTRUE -0.091524   0.028646  -3.195  0.00141 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2093 on 4574 degrees of freedom
    ## Multiple R-squared:  0.002227,   Adjusted R-squared:  0.002009 
    ## F-statistic: 10.21 on 1 and 4574 DF,  p-value: 0.001408

``` r
models_tab<-lm_table(m,smr = T,add=T,tab = models_tab,id = "Global Fst ~ Dif Expr")
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
ggplot(Global, aes(x=diff_exprs, y=Global)) + 
  geom_violin(aes(fill=diff_exprs),trim=F,scale="width") +
  xlab("Differentially Expressed") +
  ylim(NA,(1.5*max(Global$Global))) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  draw_label(paste("FALSE",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.45*max(Global$Global)),size=10) +
  draw_label(paste("TRUE",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.45*max(Global$Global)),size=10) +
  draw_label(paste("R2=",r2,"; p =",p),x=1.5,y=(1.25*max(Global$Global)),size=10) +
  guides(fill="none")+
  theme_classic()+
  ylab("Fst")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
all_lm<-lm(Global$Global~log(Global$average_tpm))
models_tab<-lm_table(all_lm,add=T,tab = models_tab,id = "Global Fst ~ expression",aic = T)
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = Global$Global ~ log(Global$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.65764 -0.12597 -0.00668  0.10978  0.84055 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             0.139485   0.007013  19.890  < 2e-16 ***
    ## log(Global$average_tpm) 0.017994   0.002263   7.952 2.29e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2081 on 4574 degrees of freedom
    ## Multiple R-squared:  0.01364,    Adjusted R-squared:  0.01342 
    ## F-statistic: 63.23 on 1 and 4574 DF,  p-value: 2.294e-15

``` r
ggplot(Global,aes(x=log(average_tpm),y=Global)) + xlab("Average TPM") + ylab("Global Fst") +
  geom_point() + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
  draw_label(paste("n =",nrow(Global),"; R2=",round(m$r.squared,3),"; p =",round(m$coefficients[2,4],3)),x=7.5,y=-0.5,size=10)+
  theme_classic()+
  ylab("Fst")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

### Does Global Fst correlate with average expression levels **in toxins**?

``` r
x<-Global_tox %>% group_by(diff_exprs) %>% summarize(mean=mean(Global),sd=sd(Global),n=length(Global))
kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
diff\_exprs
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
0.2227856
</td>
<td style="text-align:right;">
0.2614199
</td>
<td style="text-align:right;">
59
</td>
</tr>
<tr>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
0.0718318
</td>
<td style="text-align:right;">
0.1671732
</td>
<td style="text-align:right;">
9
</td>
</tr>
</tbody>
</table>

``` r
(m<-summary(lm(Global_tox$Global~Global_tox$diff_exprs)))
```

    ## 
    ## Call:
    ## lm(formula = Global_tox$Global ~ Global_tox$diff_exprs)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.49242 -0.18558  0.00877  0.16183  0.59174 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                0.22279    0.03279   6.794 3.76e-09 ***
    ## Global_tox$diff_exprsTRUE -0.15095    0.09014  -1.675   0.0987 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2519 on 66 degrees of freedom
    ## Multiple R-squared:  0.04076,    Adjusted R-squared:  0.02623 
    ## F-statistic: 2.805 on 1 and 66 DF,  p-value: 0.09872

``` r
models_tab<-lm_table(m,smr = T,add=T,tab = models_tab,id = "TOXINS Global Fst ~ Dif Expr")
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
ggplot(Global_tox, aes(x=diff_exprs, y=Global)) + 
  geom_violin(aes(fill=diff_exprs),trim=F) +
  xlab("Differentially Expressed") +
  ylim(NA,(1.5*max(Global_tox$Global))) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  draw_label(paste("FALSE",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.45*max(Global_tox$Global)),size=10) +
  draw_label(paste("TRUE",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.45*max(Global_tox$Global)),size=10) +
  draw_label(paste("R2=",r2,"; p =",p),x=1.5,y=(1.25*max(Global_tox$Global)),size=10) +
  guides(fill="none")+
  theme_classic()+
  ylab("Fst")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
all_lm<-lm(Global_tox$Global~log(Global_tox$average_tpm))
models_tab<-lm_table(all_lm,add=T,tab = models_tab,id = "TOXINS Global Fst ~ expression",aic = T)
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = Global_tox$Global ~ log(Global_tox$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.56424 -0.19727  0.00727  0.17124  0.57182 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                 -0.15786    0.19235  -0.821    0.415  
    ## log(Global_tox$average_tpm)  0.04549    0.02396   1.899    0.062 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2504 on 66 degrees of freedom
    ## Multiple R-squared:  0.0518, Adjusted R-squared:  0.03744 
    ## F-statistic: 3.606 on 1 and 66 DF,  p-value: 0.06195

``` r
print("AIC average TPM")
```

    ## [1] "AIC average TPM"

``` r
AIC(all_lm)
```

    ## [1] 8.64199

### Does Global Fst correlate with average expression levels and toxin family **in toxins**?

``` r
all_lm1<-lm(Global_tox$Global~log(Global_tox$average_tpm)+Global_tox$toxin_family)
models_tab<-lm_table(all_lm1,tf=T,add=T,tab = models_tab,id = "TOXINS global Fst ~ expression + toxin family",aic = T)
(m<-summary(all_lm1))
```

    ## 
    ## Call:
    ## lm(formula = Global_tox$Global ~ log(Global_tox$average_tpm) + 
    ##     Global_tox$toxin_family)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5276 -0.1066  0.0000  0.1588  0.4003 
    ## 
    ## Coefficients:
    ##                                        Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                            -0.31168    0.33441  -0.932   0.3557  
    ## log(Global_tox$average_tpm)             0.05909    0.03292   1.795   0.0786 .
    ## Global_tox$toxin_familyCRISP           -0.07922    0.28439  -0.279   0.7817  
    ## Global_tox$toxin_familyCTL             -0.02929    0.17123  -0.171   0.8649  
    ## Global_tox$toxin_familyHYAL             0.27695    0.31038   0.892   0.3764  
    ## Global_tox$toxin_familyKUN              0.25012    0.31267   0.800   0.4274  
    ## Global_tox$toxin_familyLAAO            -0.21785    0.28432  -0.766   0.4471  
    ## Global_tox$toxin_familyNUC              0.05482    0.29444   0.186   0.8531  
    ## Global_tox$toxin_familyPDE              0.29257    0.30442   0.961   0.3411  
    ## Global_tox$toxin_familyPLA2            -0.11511    0.17379  -0.662   0.5107  
    ## Global_tox$toxin_familyPLB             -0.05240    0.23883  -0.219   0.8272  
    ## Global_tox$toxin_familySVMPI           -0.03213    0.20110  -0.160   0.8737  
    ## Global_tox$toxin_familySVMPII          -0.02002    0.17417  -0.115   0.9089  
    ## Global_tox$toxin_familySVMPIII          0.13876    0.16367   0.848   0.4005  
    ## Global_tox$toxin_familySVSP             0.11699    0.16152   0.724   0.4722  
    ## Global_tox$toxin_familyuncharacterised  0.19373    0.29445   0.658   0.5135  
    ## Global_tox$toxin_familyVEGF             0.60590    0.28454   2.129   0.0381 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2462 on 51 degrees of freedom
    ## Multiple R-squared:  0.292,  Adjusted R-squared:  0.06988 
    ## F-statistic: 1.315 on 16 and 51 DF,  p-value: 0.225

``` r
print("AIC toxin families")
```

    ## [1] "AIC toxin families"

``` r
AIC(all_lm1)
```

    ## [1] 18.77811

### Figure 3B

``` r
(Fig3B<-ggplot(Global_tox,aes(x=log(average_tpm),y=Global)) + xlab("Average TPM") + ylab("Global Fst") +
    geom_point() + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
    draw_label(paste("n =",nrow(Global_tox),"; R2=",round(m$r.squared,2),"; p =",round(m$coefficients[2,4],2)),x=7.5,y=-0.15,size=10))+
  theme_classic()+
  ylab("Fst")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

## Tajima’s D Synonymous Mutations

get the data

``` r
Tajima<-Data[,c("ID","class","average_tpm","stdev_tpm","bcv_tpm","diff_exprs","tajimasD_synon","Total_variants","toxin_family")]
colnames(Tajima)[7]<-"tajimasD"
Tajima<-Tajima[which(Tajima$Total_variants != 0),]
Tajima<-Tajima[complete.cases(Tajima), ]
Tajima_tox<-droplevels.data.frame(subset(Tajima,Tajima$class=="Toxin"))
Tajima_nontox<-droplevels.data.frame(subset(Tajima,Tajima$class=="Nontoxin"))

x<-Tajima %>% group_by(class) %>% summarize(mean=mean(tajimasD),sd=sd(tajimasD),
                                            n=length(tajimasD))

kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
class
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Nontoxin
</td>
<td style="text-align:right;">
-0.1301358
</td>
<td style="text-align:right;">
0.7888200
</td>
<td style="text-align:right;">
4508
</td>
</tr>
<tr>
<td style="text-align:left;">
Toxin
</td>
<td style="text-align:right;">
-0.0077297
</td>
<td style="text-align:right;">
0.8446643
</td>
<td style="text-align:right;">
68
</td>
</tr>
</tbody>
</table>

### Are estimates of Tajima’s D for synonymous mutations significantly different than 0 (suggesting selection) for toxins or nontoxins?

``` r
t.test(Tajima_nontox$tajimasD)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  Tajima_nontox$tajimasD
    ## t = -11.077, df = 4507, p-value < 2.2e-16
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1531688 -0.1071028
    ## sample estimates:
    ##  mean of x 
    ## -0.1301358

``` r
models_tab <-t_table(model = t.test(Tajima_nontox$tajimasD), add=T,tab = models_tab,id="NonToxins Synonymous TD t.test")
t.test(Tajima_tox$tajimasD)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  Tajima_tox$tajimasD
    ## t = -0.075463, df = 67, p-value = 0.9401
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2121820  0.1967226
    ## sample estimates:
    ##    mean of x 
    ## -0.007729678

``` r
models_tab <-t_table(model = t.test(Tajima_tox$tajimasD), add=T,tab = models_tab,id="Toxins Synonymous TD t.test")
```

### Is there a significant difference in Tajima’s D for synonymous mutations between toxins and nontoxins?

``` r
(m<-summary(lm(Tajima$tajimasD~Tajima$class)))
```

    ## 
    ## Call:
    ## lm(formula = Tajima$tajimasD ~ Tajima$class)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8135 -0.5986 -0.0191  0.5256  2.2532 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       -0.13014    0.01176 -11.065   <2e-16 ***
    ## Tajima$classToxin  0.12241    0.09648   1.269    0.205    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7897 on 4574 degrees of freedom
    ## Multiple R-squared:  0.0003518,  Adjusted R-squared:  0.0001332 
    ## F-statistic:  1.61 on 1 and 4574 DF,  p-value: 0.2046

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="Synonymous TD ~ class")
results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(Tajima_nontox$tajimasD,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","tajimasD")
  data<-Tajima_tox[,c("class","tajimasD")]
  data<-rbind(data,samp)
  z<-summary(lm(data$tajimasD~data$class))
  results[i,1]<-round(z$r.squared,2)
  results[i,2]<-round(z$coefficients[2,4],2)
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="Tajima's D")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

``` r
#Save the results for plotting and subset toxins which fall outside the 95th percentile of nontoxin estimates.

mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
CI<-quantile(Tajima_nontox$tajimasD,c(0.025,0.975))
Tajima_tox$sig<-NA
for (i in 1:nrow(Tajima_tox)){
  if (Tajima_tox$tajimasD[i] <= CI[1] | Tajima_tox$tajimasD[i] >= CI[2]){Tajima_tox$sig[i]<-TRUE}
  else {Tajima_tox$sig[i]<-FALSE}
}
sig<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$sig=="TRUE"))
print(sig[which(sig$sig == T),c("ID","tajimasD")])
```

    ##             ID tajimasD
    ## 817     BPP_31 -1.74687
    ## 819    VEGF_30  1.62381
    ## 829    SVSP_52  1.62381
    ## 833    SVSP_72  1.48617
    ## 851    SVSP_46  1.48617
    ## 883 SVMPIII_23  1.48819

``` r
ggplot(Tajima, aes(x=class, y=tajimasD)) + 
  geom_violin(aes(fill=class),trim=F) + 
  ylim(NA,(1.5*max(Tajima$tajimasD))) +
  #scale_fill_manual(values=c("#00BFC4","#F8766D")) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  geom_point(aes(x=class,y=tajimasD),data=sig,fill="white",colour="black",size=3,pch=23) + #,position=position_jitter(width=0.1,seed=1)) +
  geom_hline(yintercept=CI,linetype="longdash",color="darkgray")  +
  geom_label(data=sig,aes(x=2.25,label=ID,y=tajimasD),label.size = 0,size=4,inherit.aes = F) + #, position=position_jitter(width=0,height=0.4,seed=1)) +
  draw_label(paste("Nontoxins",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.45*max(Tajima$tajimasD)),size=10) +
  draw_label(paste("Toxins",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.45*max(Tajima$tajimasD)),size=10) +
  draw_label(paste("R2 =",r2,"\n p =",p,"\n b =",b),x=1.5,y=(1.25*max(Tajima$tajimasD)),size=10) +
  guides(fill="none")+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->
\#\#\# Does synonymous Tajima’s D (All, Positive, and Negative values)
correlate with average expression levels and toxin families **in
toxins**?

``` r
all_lm<-lm(Tajima_tox$tajimasD~log(Tajima_tox$average_tpm))
models_tab<-lm_table(all_lm,smr=F,add=T,tab = models_tab,id="TOXINS Synonymous TD all ~ expression",aic = T)
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = Tajima_tox$tajimasD ~ log(Tajima_tox$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.79373 -0.47569 -0.08573  0.49002  1.63095 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                 -0.50218    0.65075  -0.772    0.443
    ## log(Tajima_tox$average_tpm)  0.06236    0.08105   0.769    0.444
    ## 
    ## Residual standard error: 0.8472 on 66 degrees of freedom
    ## Multiple R-squared:  0.008891,   Adjusted R-squared:  -0.006126 
    ## F-statistic: 0.5921 on 1 and 66 DF,  p-value: 0.4444

``` r
lt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD<0))
lt0_lm<-lm(lt0$tajimasD~log(lt0$average_tpm))
models_tab<-lm_table(lt0_lm,smr=F,add=T,tab = models_tab,id="TOXINS Synonymous TD lower than 0 ~ expression",aic = T)
(m_lt0<-summary(lt0_lm))
```

    ## 
    ## Call:
    ## lm(formula = lt0$tajimasD ~ log(lt0$average_tpm))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.1237 -0.4250  0.1695  0.3825  0.6956 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)          -1.35761    0.49725  -2.730   0.0105 *
    ## log(lt0$average_tpm)  0.08342    0.06215   1.342   0.1896  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4718 on 30 degrees of freedom
    ## Multiple R-squared:  0.05666,    Adjusted R-squared:  0.02521 
    ## F-statistic: 1.802 on 1 and 30 DF,  p-value: 0.1896

``` r
gt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD>0))
gt0_lm<-lm(gt0$tajimasD~log(gt0$average_tpm))
models_tab<-lm_table(gt0_lm,smr=F,add=T,tab = models_tab,id="TOXINS Synonymous TD greater than 0 ~ expression",aic = T)
(m_gt0<-summary(gt0_lm))
```

    ## 
    ## Call:
    ## lm(formula = gt0$tajimasD ~ log(gt0$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.72687 -0.45962 -0.01307  0.44653  0.76313 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)          -0.17177    0.69775  -0.246    0.808
    ## log(gt0$average_tpm)  0.13007    0.08872   1.466    0.156
    ## 
    ## Residual standard error: 0.5056 on 24 degrees of freedom
    ## Multiple R-squared:  0.0822, Adjusted R-squared:  0.04396 
    ## F-statistic: 2.149 on 1 and 24 DF,  p-value: 0.1556

``` r
print("all values")
```

    ## [1] "all values"

``` r
AIC(all_lm)
```

    ## [1] 174.4019

``` r
print("lower than 0 values")
```

    ## [1] "lower than 0 values"

``` r
AIC(lt0_lm)
```

    ## [1] 46.66958

``` r
print("greater than 0 values")
```

    ## [1] "greater than 0 values"

``` r
AIC(gt0_lm)
```

    ## [1] 42.23432

``` r
all_lm1<-lm(Tajima_tox$tajimasD~log(Tajima_tox$average_tpm)+Tajima_tox$toxin_family)
(m1<-summary(all_lm1))
```

    ## 
    ## Call:
    ## lm(formula = Tajima_tox$tajimasD ~ log(Tajima_tox$average_tpm) + 
    ##     Tajima_tox$toxin_family)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7663 -0.4407  0.0000  0.3062  1.5046 
    ## 
    ## Coefficients:
    ##                                        Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)                            -1.81322    1.07970  -1.679   0.0992 . 
    ## log(Tajima_tox$average_tpm)             0.09599    0.10628   0.903   0.3707   
    ## Tajima_tox$toxin_familyCRISP           -0.18472    0.91819  -0.201   0.8414   
    ## Tajima_tox$toxin_familyCTL              0.95759    0.55284   1.732   0.0893 . 
    ## Tajima_tox$toxin_familyHYAL             1.55865    1.00211   1.555   0.1260   
    ## Tajima_tox$toxin_familyKUN             -0.07452    1.00950  -0.074   0.9414   
    ## Tajima_tox$toxin_familyLAAO             0.53369    0.91799   0.581   0.5636   
    ## Tajima_tox$toxin_familyNUC              0.59203    0.95066   0.623   0.5362   
    ## Tajima_tox$toxin_familyPDE              0.80491    0.98289   0.819   0.4166   
    ## Tajima_tox$toxin_familyPLA2             0.92252    0.56113   1.644   0.1063   
    ## Tajima_tox$toxin_familyPLB              1.50838    0.77110   1.956   0.0559 . 
    ## Tajima_tox$toxin_familySVMPI            0.10708    0.64928   0.165   0.8697   
    ## Tajima_tox$toxin_familySVMPII           1.04790    0.56233   1.863   0.0682 . 
    ## Tajima_tox$toxin_familySVMPIII          1.44652    0.52843   2.737   0.0085 **
    ## Tajima_tox$toxin_familySVSP             1.31405    0.52151   2.520   0.0149 * 
    ## Tajima_tox$toxin_familyuncharacterised  1.17618    0.95067   1.237   0.2217   
    ## Tajima_tox$toxin_familyVEGF             2.59173    0.91870   2.821   0.0068 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7948 on 51 degrees of freedom
    ## Multiple R-squared:  0.326,  Adjusted R-squared:  0.1146 
    ## F-statistic: 1.542 on 16 and 51 DF,  p-value: 0.1214

``` r
models_tab<-lm_table(all_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS Synonymous TD all ~ expression + toxin family",aic = T)
lt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD<0))
lt0_lm1<-lm(lt0$tajimasD~log(lt0$average_tpm)+lt0$toxin_family)
models_tab<-lm_table(lt0_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS Synonymous TD lower than 0 ~ expression + toxin family",aic = T)
(m_lt01<-summary(lt0_lm1))
```

    ## 
    ## Call:
    ## lm(formula = lt0$tajimasD ~ log(lt0$average_tpm) + lt0$toxin_family)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.80760 -0.12220  0.01127  0.18698  0.61934 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -2.17746    0.84882  -2.565   0.0195 *
    ## log(lt0$average_tpm)     0.13560    0.08729   1.553   0.1377  
    ## lt0$toxin_familyCRISP   -0.17428    0.55231  -0.316   0.7560  
    ## lt0$toxin_familyCTL      0.46680    0.37335   1.250   0.2272  
    ## lt0$toxin_familyKUN      0.08219    0.65099   0.126   0.9009  
    ## lt0$toxin_familyLAAO     0.52614    0.55208   0.953   0.3532  
    ## lt0$toxin_familyNUC      0.68443    0.58820   1.164   0.2598  
    ## lt0$toxin_familyPDE      0.93603    0.62291   1.503   0.1503  
    ## lt0$toxin_familyPLA2     0.61539    0.40751   1.510   0.1484  
    ## lt0$toxin_familyPLB      0.91358    0.59973   1.523   0.1451  
    ## lt0$toxin_familySVMPI   -0.20863    0.44345  -0.470   0.6437  
    ## lt0$toxin_familySVMPII   0.66790    0.39139   1.707   0.1051  
    ## lt0$toxin_familySVMPIII  0.52041    0.39879   1.305   0.2083  
    ## lt0$toxin_familySVSP     0.38365    0.36160   1.061   0.3027  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4779 on 18 degrees of freedom
    ## Multiple R-squared:  0.4193, Adjusted R-squared:  -0.0001529 
    ## F-statistic: 0.9996 on 13 and 18 DF,  p-value: 0.4891

``` r
gt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD>0))
gt0_lm1<-lm(gt0$tajimasD~log(gt0$average_tpm)+gt0$toxin_family)
models_tab<-lm_table(gt0_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS Synonymous TD greater than 0 ~ expression + toxin family",aic = T)
(m_gt01<-summary(gt0_lm1))
```

    ## 
    ## Call:
    ## lm(formula = gt0$tajimasD ~ log(gt0$average_tpm) + gt0$toxin_family)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.8359 -0.3749  0.0000  0.3786  0.6779 
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                      0.32125    1.16819   0.275    0.787
    ## log(gt0$average_tpm)             0.05630    0.16971   0.332    0.744
    ## gt0$toxin_familyHYAL            -0.36109    0.66448  -0.543    0.594
    ## gt0$toxin_familyPLA2             0.03059    0.64132   0.048    0.963
    ## gt0$toxin_familyPLB              0.35208    0.63425   0.555    0.586
    ## gt0$toxin_familySVMPII          -0.04737    0.61461  -0.077    0.940
    ## gt0$toxin_familySVMPIII          0.06650    0.50494   0.132    0.897
    ## gt0$toxin_familySVSP             0.17779    0.42107   0.422    0.678
    ## gt0$toxin_familyuncharacterised -0.68587    0.63279  -1.084    0.294
    ## gt0$toxin_familyVEGF             0.80682    0.73181   1.102    0.287
    ## 
    ## Residual standard error: 0.5469 on 16 degrees of freedom
    ## Multiple R-squared:  0.2839, Adjusted R-squared:  -0.1189 
    ## F-statistic: 0.7048 on 9 and 16 DF,  p-value: 0.6969

``` r
print("all values")
```

    ## [1] "all values"

``` r
AIC(all_lm1)
```

    ## [1] 178.1801

``` r
print("lower than 0 values")
```

    ## [1] "lower than 0 values"

``` r
AIC(lt0_lm1)
```

    ## [1] 55.14524

``` r
print("greater than 0 values")
```

    ## [1] "greater than 0 values"

``` r
AIC(gt0_lm1)
```

    ## [1] 51.782

``` r
(testplot<-ggplot(Tajima_tox,aes(x=log(average_tpm),y=tajimasD, color = toxin_family)) + xlab("Average TPM") + ylab("Tajima's D") +
    geom_point() + #geom_hline(yintercept=0) + 
    geom_abline(intercept = coef(all_lm1)[1],slope=coef(all_lm1)[2],lty=2) +
    geom_abline(intercept = coef(lt0_lm1)[1],slope=coef(lt0_lm1)[2],lty=3) +
    geom_abline(intercept = coef(gt0_lm1)[1],slope=coef(gt0_lm1)[2],lty=2) +
    draw_label(paste("n =",nrow(Tajima_tox),"; R2=",round(m1$r.squared,2),"; p =",round(m1$coefficients[2,4],2)),x=7.5,y=0,size=10) +
    draw_label(paste("n =",nrow(lt0),"; R2=",round(m_lt01$r.squared,2),"; p =",round(m_lt01$coefficients[2,4],2)),x=7.5,y=-1,size=10) +
    draw_label(paste("n =",nrow(gt0),"; R2=",round(m_gt01$r.squared,2),"; p =",round(m_gt01$coefficients[2,4],2)),x=7.5,y=1,size=10))#+
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

``` r
 # theme_classic()
```

``` r
(Fig3A<-ggplot(Tajima_tox,aes(x=log(average_tpm),y=tajimasD)) + xlab("Average TPM") + ylab("Tajima's D") +
    geom_point() + #geom_hline(yintercept=0) + 
    geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
    geom_abline(intercept = coef(lt0_lm)[1],slope=coef(lt0_lm)[2],lty=3) +
    geom_abline(intercept = coef(gt0_lm)[1],slope=coef(gt0_lm)[2],lty=2) +
    draw_label(paste("n =",nrow(Tajima_tox),"; R2=",round(m$r.squared,2),"; p =",round(m$coefficients[2,4],2)),x=7.5,y=0,size=10) +
    draw_label(paste("n =",nrow(lt0),"; R2=",round(m_lt0$r.squared,2),"; p =",round(m_lt0$coefficients[2,4],2)),x=7.5,y=-1,size=10) +
    draw_label(paste("n =",nrow(gt0),"; R2=",round(m_gt0$r.squared,2),"; p =",round(m_gt0$coefficients[2,4],2)),x=7.5,y=1,size=10))#+
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

``` r
 # theme_classic()
```

## Tajima’s D Nonsynonymous Mutations

get data

``` r
Tajima<-Data[,c("ID","class","average_tpm","stdev_tpm","bcv_tpm","diff_exprs","tajimasD_nonsyn",'Total_variants',"toxin_family")]
colnames(Tajima)[7]<-"tajimasD"
Tajima<-Tajima[which(Tajima$Total_variants != 0),]
Tajima<-Tajima[complete.cases(Tajima), ]
Tajima_tox<-droplevels.data.frame(subset(Tajima,Tajima$class=="Toxin"))
Tajima_nontox<-droplevels.data.frame(subset(Tajima,Tajima$class=="Nontoxin"))

x<-Tajima %>% group_by(class) %>% summarize(mean=mean(tajimasD),sd=sd(tajimasD),
                                            n=length(tajimasD))

kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
class
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Nontoxin
</td>
<td style="text-align:right;">
-0.1883360
</td>
<td style="text-align:right;">
0.7085367
</td>
<td style="text-align:right;">
4508
</td>
</tr>
<tr>
<td style="text-align:left;">
Toxin
</td>
<td style="text-align:right;">
-0.0051722
</td>
<td style="text-align:right;">
0.7740596
</td>
<td style="text-align:right;">
68
</td>
</tr>
</tbody>
</table>

### Are estimates of Tajima’s D for nonsynonymous mutations significantly different than 0 (suggesting selection) for toxins or nontoxins?

``` r
t.test(Tajima_nontox$tajimasD)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  Tajima_nontox$tajimasD
    ## t = -17.847, df = 4507, p-value < 2.2e-16
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.2090248 -0.1676473
    ## sample estimates:
    ## mean of x 
    ## -0.188336

``` r
models_tab <-t_table(model = t.test(Tajima_nontox$tajimasD), add=T,tab = models_tab,id="NonToxins Nonynonymous TD t.test")
t.test(Tajima_tox$tajimasD)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  Tajima_tox$tajimasD
    ## t = -0.055101, df = 67, p-value = 0.9562
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1925346  0.1821901
    ## sample estimates:
    ##    mean of x 
    ## -0.005172244

``` r
models_tab <-t_table(model = t.test(Tajima_tox$tajimasD), add=T,tab = models_tab,id="Toxins Nonynonymous TD t.test")
```

### Is there a significant difference in Tajima’s D for nonsynonymous mutations between toxins and nontoxins?

``` r
(m<-summary(lm(Tajima$tajimasD~Tajima$class)))
```

    ## 
    ## Call:
    ## lm(formula = Tajima$tajimasD ~ Tajima$class)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7553 -0.5404  0.1883  0.1883  2.3037 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       -0.18834    0.01057 -17.822   <2e-16 ***
    ## Tajima$classToxin  0.18316    0.08669   2.113   0.0347 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7095 on 4574 degrees of freedom
    ## Multiple R-squared:  0.000975,   Adjusted R-squared:  0.0007566 
    ## F-statistic: 4.464 on 1 and 4574 DF,  p-value: 0.03467

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="Nonsynonymous TD (non normal**) ~ class ")
results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(Tajima_nontox$tajimasD,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","tajimasD")
  data<-Tajima_tox[,c("class","tajimasD")]
  data<-rbind(data,samp)
  z<-summary(lm(data$tajimasD~data$class))
  results[i,1]<-round(z$r.squared,2)
  results[i,2]<-round(z$coefficients[2,4],2)
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="Tajima's D")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

``` r
#Save the results for plotting and subset toxins which fall outside the 95th percentile of nontoxin estimates.

mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
CI<-quantile(Tajima_nontox$tajimasD,c(0.025,0.975))
Tajima_tox$sig<-NA
for (i in 1:nrow(Tajima_tox)){
  if (Tajima_tox$tajimasD[i] <= CI[1] | Tajima_tox$tajimasD[i] >= CI[2]){Tajima_tox$sig[i]<-TRUE}
  else {Tajima_tox$sig[i]<-FALSE}
}
sig<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$sig=="TRUE"))
print(sig[which(sig$sig == T),c("ID","tajimasD")])
```

    ##            ID tajimasD
    ## 819   VEGF_30  1.38110
    ## 820    CTL_78 -1.45138
    ## 821 SVMPII_12  1.48617
    ## 834   SVSP_66 -1.62929
    ## 849   PLA2_43 -1.45138
    ## 870   SVSP_79  2.02297
    ## 873  SVMPI_44 -1.45138

``` r
ggplot(Tajima, aes(x=class, y=tajimasD)) + 
  geom_violin(aes(fill=class),trim=F) + 
  ylim(NA,(1.5*max(Tajima$tajimasD))) +
  #scale_fill_manual(values=c("#00BFC4","#F8766D")) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  geom_point(aes(x=class,y=tajimasD),data=sig,fill="white",colour="black",size=3,pch=23) + #,position=position_jitter(width=0.1,seed=1)) +
  geom_hline(yintercept=CI,linetype="longdash",color="darkgray")  +
  geom_label(data=sig,aes(x=2.25,label=ID,y=tajimasD),label.size = 0,size=4,inherit.aes = F) + #, position=position_jitter(width=0,height=0.4,seed=1)) +
  draw_label(paste("Nontoxins",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.45*max(Tajima$tajimasD)),size=10) +
  draw_label(paste("Toxins",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.45*max(Tajima$tajimasD)),size=10) +
  draw_label(paste("R2 =",r2,"\n p =",p,"\n b =",b),x=1.5,y=(1.25*max(Tajima$tajimasD)),size=10) +
  guides(fill="none")+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->
\#\#\# Non Parametric

``` r
kruskal.test(Tajima$tajimasD~Tajima$class)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Tajima$tajimasD by Tajima$class
    ## Kruskal-Wallis chi-squared = 1.84, df = 1, p-value = 0.1749

``` r
models_tab<-chi2_table(model=kruskal.test(Tajima$tajimasD~Tajima$class),add = T,tab = models_tab,id= "Nonsynonymous TD ~ class, Kruskal test")
```

``` r
results<-matrix(nrow=1000,ncol=2)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(Tajima_nontox$tajimasD,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","tajimasD")
  data<-Tajima_tox[,c("class","tajimasD")]
  data<-rbind(data,samp)
  z<-kruskal.test(data$tajimasD~data$class)
  results[i,1]<-round(z$statistic,2)
  results[i,2]<-round(z$p.value,2)
}
table(results[,2] < 0.05)
```

    ## 
    ## FALSE 
    ##  1000

### Does nonsynonymous Tajima’s D (All, Positive, and Negative values) correlate with average expression levels and toxin families **in toxins**?

``` r
all_lm<-lm(Tajima_tox$tajimasD~log(Tajima_tox$average_tpm))
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = Tajima_tox$tajimasD ~ log(Tajima_tox$average_tpm))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.58339 -0.43957 -0.04252  0.54664  2.13552 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                 -1.03763    0.58504  -1.774   0.0807 .
    ## log(Tajima_tox$average_tpm)  0.13022    0.07286   1.787   0.0785 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7617 on 66 degrees of freedom
    ## Multiple R-squared:  0.04616,    Adjusted R-squared:  0.03171 
    ## F-statistic: 3.194 on 1 and 66 DF,  p-value: 0.0785

``` r
models_tab<-lm_table(all_lm,smr=F,add=T,tab = models_tab,id="TOXINS Nonynonymous TD all ~ expression",aic = T)
lt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD<0))
lt0_lm<-lm(lt0$tajimasD~log(lt0$average_tpm))
(m_lt0<-summary(lt0_lm))
```

    ## 
    ## Call:
    ## lm(formula = lt0$tajimasD ~ log(lt0$average_tpm))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0195 -0.2544  0.1440  0.3726  0.5462 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)          -1.28242    0.45793  -2.800  0.00846 **
    ## log(lt0$average_tpm)  0.08832    0.05814   1.519  0.13826   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4613 on 33 degrees of freedom
    ## Multiple R-squared:  0.06536,    Adjusted R-squared:  0.03704 
    ## F-statistic: 2.308 on 1 and 33 DF,  p-value: 0.1383

``` r
models_tab<-lm_table(lt0_lm,smr=F,add=T,tab = models_tab,id="TOXINS Nonynonymous TD lower than 0 ~ expression",aic = T)
gt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD>0))
gt0_lm<-lm(gt0$tajimasD~log(gt0$average_tpm))
models_tab<-lm_table(gt0_lm,smr=F,add=T,tab = models_tab,id="TOXINS Nonynonymous TD greater than 0 ~ expression",aic = T)
(m_gt0<-summary(gt0_lm))
```

    ## 
    ## Call:
    ## lm(formula = gt0$tajimasD ~ log(gt0$average_tpm))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5158 -0.2511 -0.0090  0.1421  1.1885 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)           1.10861    0.52950   2.094    0.047 *
    ## log(gt0$average_tpm) -0.03859    0.06342  -0.608    0.549  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4017 on 24 degrees of freedom
    ## Multiple R-squared:  0.01519,    Adjusted R-squared:  -0.02585 
    ## F-statistic: 0.3701 on 1 and 24 DF,  p-value: 0.5487

``` r
print("all values")
```

    ## [1] "all values"

``` r
AIC(all_lm)
```

    ## [1] 159.9241

``` r
print("lower than 0 values")
```

    ## [1] "lower than 0 values"

``` r
AIC(lt0_lm)
```

    ## [1] 49.10445

``` r
print("greater than 0 values")
```

    ## [1] "greater than 0 values"

``` r
AIC(gt0_lm)
```

    ## [1] 30.27637

``` r
all_lm1<-lm(Tajima_tox$tajimasD~log(Tajima_tox$average_tpm)+Tajima_tox$toxin_family)
models_tab<-lm_table(all_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS Nonsynonymous TD all ~ expression + toxin family",aic = T)
(m1<-summary(all_lm1))
```

    ## 
    ## Call:
    ## lm(formula = Tajima_tox$tajimasD ~ log(Tajima_tox$average_tpm) + 
    ##     Tajima_tox$toxin_family)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.77480 -0.26831 -0.00999  0.36543  1.94286 
    ## 
    ## Coefficients:
    ##                                        Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                            -0.61514    1.06172  -0.579    0.565  
    ## log(Tajima_tox$average_tpm)             0.12778    0.10451   1.223    0.227  
    ## Tajima_tox$toxin_familyCRISP           -0.80539    0.90290  -0.892    0.377  
    ## Tajima_tox$toxin_familyCTL             -0.52507    0.54363  -0.966    0.339  
    ## Tajima_tox$toxin_familyHYAL            -0.73671    0.98542  -0.748    0.458  
    ## Tajima_tox$toxin_familyKUN             -1.15753    0.99268  -1.166    0.249  
    ## Tajima_tox$toxin_familyLAAO            -1.65823    0.90269  -1.837    0.072 .
    ## Tajima_tox$toxin_familyNUC              0.06622    0.93483   0.071    0.944  
    ## Tajima_tox$toxin_familyPDE              0.53102    0.96651   0.549    0.585  
    ## Tajima_tox$toxin_familyPLA2            -0.65087    0.55178  -1.180    0.244  
    ## Tajima_tox$toxin_familyPLB             -0.49078    0.75826  -0.647    0.520  
    ## Tajima_tox$toxin_familySVMPI           -0.52739    0.63846  -0.826    0.413  
    ## Tajima_tox$toxin_familySVMPII          -0.58291    0.55296  -1.054    0.297  
    ## Tajima_tox$toxin_familySVMPIII         -0.34197    0.51963  -0.658    0.513  
    ## Tajima_tox$toxin_familySVSP            -0.21250    0.51282  -0.414    0.680  
    ## Tajima_tox$toxin_familyuncharacterised -0.26179    0.93483  -0.280    0.781  
    ## Tajima_tox$toxin_familyVEGF             0.87102    0.90339   0.964    0.340  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7816 on 51 degrees of freedom
    ## Multiple R-squared:  0.224,  Adjusted R-squared:  -0.01949 
    ## F-statistic:  0.92 on 16 and 51 DF,  p-value: 0.5524

``` r
lt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD<0))
lt0_lm1<-lm(lt0$tajimasD~log(lt0$average_tpm)+lt0$toxin_family)
models_tab<-lm_table(lt0_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS Nonsynonymous TD lower than 0 ~ expression + toxin family",aic = T)
(m_lt01<-summary(lt0_lm1))
```

    ## 
    ## Call:
    ## lm(formula = lt0$tajimasD ~ log(lt0$average_tpm) + lt0$toxin_family)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.95753 -0.01753  0.03892  0.22370  0.56441 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -1.61717    0.80725  -2.003   0.0571 .
    ## log(lt0$average_tpm)     0.14980    0.07589   1.974   0.0605 .
    ## lt0$toxin_familyCTL      0.05554    0.47858   0.116   0.9086  
    ## lt0$toxin_familyHYAL     0.14624    0.67523   0.217   0.8305  
    ## lt0$toxin_familyKUN     -0.27085    0.68044  -0.398   0.6943  
    ## lt0$toxin_familyLAAO    -0.86284    0.62102  -1.389   0.1780  
    ## lt0$toxin_familyPLA2    -0.36469    0.49534  -0.736   0.4690  
    ## lt0$toxin_familyPLB      0.36289    0.56218   0.646   0.5250  
    ## lt0$toxin_familySVMPI   -1.06056    0.62263  -1.703   0.1020  
    ## lt0$toxin_familySVMPII  -0.06881    0.48030  -0.143   0.8873  
    ## lt0$toxin_familySVMPIII -0.15862    0.49239  -0.322   0.7503  
    ## lt0$toxin_familySVSP    -0.35592    0.48693  -0.731   0.4722  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4384 on 23 degrees of freedom
    ## Multiple R-squared:  0.4115, Adjusted R-squared:   0.13 
    ## F-statistic: 1.462 on 11 and 23 DF,  p-value: 0.2127

``` r
gt0<-droplevels.data.frame(subset(Tajima_tox,Tajima_tox$tajimasD>0))
gt0_lm1<-lm(gt0$tajimasD~log(gt0$average_tpm)+gt0$toxin_family)
models_tab<-lm_table(gt0_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS Nonsynonymous TD greater than 0 ~ expression + toxin family",aic = T)
(m_gt01<-summary(gt0_lm1))
```

    ## 
    ## Call:
    ## lm(formula = gt0$tajimasD ~ log(gt0$average_tpm) + gt0$toxin_family)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5205 -0.1354  0.0000  0.1186  0.9954 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)              1.44746    0.81912   1.767   0.0975 .
    ## log(gt0$average_tpm)    -0.09652    0.08611  -1.121   0.2800  
    ## gt0$toxin_familyCTL      0.16861    0.39478   0.427   0.6754  
    ## gt0$toxin_familyNUC     -0.45702    0.46500  -0.983   0.3413  
    ## gt0$toxin_familyPDE     -0.21146    0.50707  -0.417   0.6826  
    ## gt0$toxin_familyPLA2    -0.15540    0.30879  -0.503   0.6221  
    ## gt0$toxin_familySVMPI    0.19893    0.33207   0.599   0.5581  
    ## gt0$toxin_familySVMPII   0.96384    0.42071   2.291   0.0369 *
    ## gt0$toxin_familySVMPIII  0.08354    0.26736   0.312   0.7590  
    ## gt0$toxin_familySVSP     0.26575    0.27546   0.965   0.3500  
    ## gt0$toxin_familyVEGF     0.78358    0.42071   1.862   0.0822 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3632 on 15 degrees of freedom
    ## Multiple R-squared:  0.4968, Adjusted R-squared:  0.1614 
    ## F-statistic: 1.481 on 10 and 15 DF,  p-value: 0.2381

``` r
print("all values")
```

    ## [1] "all values"

``` r
AIC(all_lm1)
```

    ## [1] 175.8953

``` r
print("lower than 0 values")
```

    ## [1] "lower than 0 values"

``` r
AIC(lt0_lm1)
```

    ## [1] 52.91496

``` r
print("greater than 0 values")
```

    ## [1] "greater than 0 values"

``` r
AIC(gt0_lm1)
```

    ## [1] 30.81612

``` r
(testplot<-ggplot(Tajima_tox,aes(x=log(average_tpm),y=tajimasD, color = toxin_family)) + xlab("Average TPM") + ylab("Tajima's D") +
    geom_point() + #geom_hline(yintercept=0) + 
    geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
    geom_abline(intercept = coef(lt0_lm)[1],slope=coef(lt0_lm)[2],lty=3) +
    geom_abline(intercept = coef(gt0_lm)[1],slope=coef(gt0_lm)[2],lty=2) +
    draw_label(paste("n =",nrow(Tajima_tox),"; R2=",round(m$r.squared,2),"; p =",round(m$coefficients[2,4],2)),x=7.5,y=0,size=10) +
    draw_label(paste("n =",nrow(lt0),"; R2=",round(m_lt0$r.squared,2),"; p =",round(m_lt0$coefficients[2,4],2)),x=7.5,y=-1,size=10) +
    draw_label(paste("n =",nrow(gt0),"; R2=",round(m_gt0$r.squared,2),"; p =",round(m_gt0$coefficients[2,4],2)),x=7.5,y=1,size=10))+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

``` r
(Fig3A<-ggplot(Tajima_tox,aes(x=log(average_tpm),y=tajimasD)) + xlab("Average TPM") + ylab("Tajima's D") +
    geom_point() + #geom_hline(yintercept=0) + 
    geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
    geom_abline(intercept = coef(lt0_lm)[1],slope=coef(lt0_lm)[2],lty=3) +
    geom_abline(intercept = coef(gt0_lm)[1],slope=coef(gt0_lm)[2],lty=2) +
    draw_label(paste("n =",nrow(Tajima_tox),"; R2=",round(m$r.squared,2),"; p =",round(m$coefficients[2,4],2)),x=7.5,y=0,size=10) +
    draw_label(paste("n =",nrow(lt0),"; R2=",round(m_lt0$r.squared,2),"; p =",round(m_lt0$coefficients[2,4],2)),x=7.5,y=-1,size=10) +
    draw_label(paste("n =",nrow(gt0),"; R2=",round(m_gt0$r.squared,2),"; p =",round(m_gt0$coefficients[2,4],2)),x=7.5,y=1,size=10))+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-54-1.png)<!-- -->

## HyPhy Busted

get data

``` r
BUSTED<-Data[,c("ID","class","toxin_family","Total_variants","average_tpm","stdev_tpm","bcv_tpm","diff_exprs","p_value","LRT","UCW3")]

genelist<-scan('list_to_R', character(), quote = '')

for (i in 1:length(genelist)){
  ind<- which(BUSTED$ID == genelist[i])
  BUSTED<-BUSTED[-ind,]
}

#colnames(BUSTED)[ncol(BUSTED)]<- 'LRT'
BUSTED<-BUSTED[which(BUSTED$Total_variants != 0),]

BUSTED<-BUSTED[complete.cases(BUSTED), ]
BUSTED_tox<-droplevels.data.frame(subset(BUSTED,BUSTED$class=="Toxin"))
BUSTED_nontox<-droplevels.data.frame(subset(BUSTED,BUSTED$class=="Nontoxin"))

x<-BUSTED %>% group_by(class) %>% summarize(mean=mean(LRT),sd=sd(LRT),n=length(LRT))

kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
class
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Nontoxin
</td>
<td style="text-align:right;">
0.2869356
</td>
<td style="text-align:right;">
1.037574
</td>
<td style="text-align:right;">
4497
</td>
</tr>
<tr>
<td style="text-align:left;">
Toxin
</td>
<td style="text-align:right;">
0.6380741
</td>
<td style="text-align:right;">
1.458140
</td>
<td style="text-align:right;">
68
</td>
</tr>
</tbody>
</table>

### Is there a significant difference in BUSTED between toxins and nontoxins?

#### Non Parametric test

independent 2-group Mann-Whitney U Test wilcox.test(y\~A) where y is
numeric and A is A binary factor

Kruskal Wallis Test One Way Anova by Ranks kruskal.test(y\~A) \# where
y1 is numeric and A is a factor

``` r
wilcox.test(BUSTED$LRT~BUSTED$class)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  BUSTED$LRT by BUSTED$class
    ## W = 126159, p-value = 0.0005966
    ## alternative hypothesis: true location shift is not equal to 0

``` r
kruskal.test(BUSTED$LRT~BUSTED$class)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  BUSTED$LRT by BUSTED$class
    ## Kruskal-Wallis chi-squared = 11.787, df = 1, p-value = 0.0005965

``` r
models_tab<-chi2_table(model = kruskal.test(BUSTED$LRT~BUSTED$class),add = T,tab = models_tab,id = "BUSTED LRT ~ class, Kruskal test")
```

``` r
results<-matrix(nrow=1000,ncol=3)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(BUSTED_nontox$LRT,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","LRT")
  data<-BUSTED_tox[,c("class","LRT")]
  data<-rbind(data,samp)
  z<-kruskal.test(data$LRT~data$class)
  results[i,1]<-round(z$statistic,2)
  results[i,2]<-round(z$p.value,2)
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="BUSTED Likelihood Ratio")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

#### SLM

``` r
library(car)
```

    ## Loading required package: carData

    ## 
    ## Attaching package: 'car'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode

``` r
library(compositions)
#(m<-summary(lm(LRT~class,data = BUSTED[which(BUSTED$LRT > 0),])))
(m<-summary(lm(LRT~class,data = BUSTED)))
```

    ## 
    ## Call:
    ## lm(formula = LRT ~ class, data = BUSTED)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -5.8850 -0.2869 -0.2869 -0.2869 12.3574 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.28694    0.01558   18.41  < 2e-16 ***
    ## classToxin   0.35114    0.12768    2.75  0.00598 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.045 on 4563 degrees of freedom
    ## Multiple R-squared:  0.001655,   Adjusted R-squared:  0.001436 
    ## F-statistic: 7.564 on 1 and 4563 DF,  p-value: 0.005979

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="BUSTED LRT(non normal **) ~ class")
results<-matrix(nrow=1000,ncol=3)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample(BUSTED_nontox$LRT,x$n[2], replace=T)
  samp<-cbind(rep("Nontoxin",x$n[2]),samp)
  colnames(samp)<-c("class","LRT")
  data<-BUSTED_tox[,c("class","LRT")]
  data<-rbind(data,samp)
  z<-summary(lm(data$LRT~data$class))
  results[i,1]<-round(z$r.squared,2)
  results[i,2]<-round(z$coefficients[2,4],2)
  t<-shapiro.test(z$residuals)
  results[i,3]<-t$p.value
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="BUSTED Likelihood Ratio")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

``` r
b1<-round(pct[2]*0.01,2)
```

``` r
#Save the results for plotting and subset toxins which fall outside the 95th percentile of nontoxin estimates.
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]


p<-round(kruskal.test(BUSTED$LRT~BUSTED$class)$p.value,2)

chi2<-round(kruskal.test(BUSTED$LRT~BUSTED$class)$statistic,2)
CI<-quantile(BUSTED_nontox$LRT,c(0.95))
BUSTED_tox$sig<-NA
for (i in 1:nrow(BUSTED_tox)){
  if (BUSTED_tox$LRT[i] >= CI){BUSTED_tox$sig[i]<-TRUE}
  else {BUSTED_tox$sig[i]<-FALSE}
}
sig<-droplevels.data.frame(subset(BUSTED_tox,BUSTED_tox$sig=="TRUE"))
citx<-sig[which(sig$sig == T),c("ID","class","toxin_family","p_value","LRT","UCW3")]
citx[order(citx$toxin_family),]
```

    ##             ID class toxin_family     p_value      LRT       UCW3
    ## 848     CTL_86 Toxin          CTL 0.099788179 3.223117   34.50596
    ## 861     CTL_80 Toxin          CTL 0.107366619 3.076718 8877.19394
    ## 862   HYAL_112 Toxin         HYAL 0.094410191 3.333918   96.05210
    ## 826    SVMPI_2 Toxin        SVMPI 0.076789867 3.747071   32.05170
    ## 881   SVMPI_40 Toxin        SVMPI 0.091081231 3.405713   30.41514
    ## 846  SVMPII_15 Toxin       SVMPII 0.092293762 3.379263  407.99673
    ## 841 SVMPIII_39 Toxin      SVMPIII 0.005254159 9.111176  341.56629
    ## 844    SVSP_58 Toxin         SVSP 0.142038429 2.517021  583.88977

### Figure 2C

``` r
(Fig2C <- ggplot(BUSTED, aes(x=class, y=LRT)) + 
    geom_violin(aes(fill=class),trim=F,scale = "width") + 
    ylim(NA,(1.25*max(BUSTED$LRT))) +
    #scale_fill_manual(values=c("#00BFC4","#F8766D")) +
    scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
    stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
    geom_point(aes(x=class,y=LRT),data=sig,fill="white",colour="black",size=3,pch=23) + #,position=position_jitter(width=0.1,seed=1)) +
    geom_hline(yintercept=CI,linetype="longdash",color="darkgray")  +
    labs(y = "LRT_Busted", x= "Class")+
    geom_label(data=sig,aes(x=2.25,label=ID,y=LRT),label.size = 0,size=4,inherit.aes = F) +
#    position=position_jitter(width=0,height=0.4,seed=1)) +
    draw_label(paste("Nontoxins",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.25*max(BUSTED$LRT)),size=10) +
    draw_label(paste("Toxins",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.25*max(BUSTED$LRT)),size=10) +
    draw_label(paste("Chi^2 =",chi2,"\n p =",p,"\n b =",b),x=1.5,y=(max(BUSTED$LRT)),size=10) +
    guides(fill="none"))+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

### Do BUSTED Likelihood Ratios correlate with average expression levels?

#### Non parametric

#### Non parametric test

``` r
kruskal.test(BUSTED$LRT~BUSTED$diff_exprs)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  BUSTED$LRT by BUSTED$diff_exprs
    ## Kruskal-Wallis chi-squared = 0.062586, df = 1, p-value = 0.8025

``` r
models_tab<-chi2_table(model = kruskal.test(BUSTED$LRT~BUSTED$diff_exprs),add = T,tab = models_tab,id = "BUSTED LRT ~ Dif Exp, Kruskal test")
```

``` r
results<-matrix(nrow=1000,ncol=3)
set.seed(2019)
for(i in 1:1000) {
  samp<-sample_n(BUSTED_nontox[,c('diff_exprs','LRT')],length(which(BUSTED$diff_exprs == "TRUE")), replace=T)
# samp<-cbind(rep("Nontoxin",x$n[2]),samp)
#  colnames(samp)<-c("class","LRT")
  data<-BUSTED[which(BUSTED$diff_exprs == "TRUE"),c("diff_exprs","LRT")]
  data<-rbind(data,samp)
  z<-kruskal.test(data$LRT~data$diff_exprs)
  results[i,1]<-round(z$statistic,2)
  results[i,2]<-round(z$p.value,2)
}
pie_results<-unlist(as.data.frame(table(results[,2] < 0.05),row.names = c("FALSE","TRUE"))[2])
lbls<-c("False","True")
pct<-as.vector(round(pie_results/sum(pie_results)*100,2))
lbls<-paste(lbls,pct)
lbls<-paste(lbls,"%",sep="")
pie(pie_results,labels=lbls,col=c("red","green"),main="BUSTED Likelihood Ratio")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-62-1.png)<!-- -->

``` r
b<-round(pct[2]*0.01,2)
```

#### SLM

``` r
x<-BUSTED %>% group_by(diff_exprs) %>% summarize(mean=mean(LRT),sd=sd(LRT),n=length(LRT))
kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
diff\_exprs
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
0.2944282
</td>
<td style="text-align:right;">
1.0509219
</td>
<td style="text-align:right;">
4511
</td>
</tr>
<tr>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
0.1031940
</td>
<td style="text-align:right;">
0.3876709
</td>
<td style="text-align:right;">
54
</td>
</tr>
</tbody>
</table>

``` r
(m<-summary(lm(BUSTED$LRT~BUSTED$diff_exprs)))
```

    ## 
    ## Call:
    ## lm(formula = BUSTED$LRT ~ BUSTED$diff_exprs)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -5.8925 -0.2944 -0.2944 -0.2944 12.3499 
    ## 
    ## Coefficients:
    ##                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            0.29443    0.01557  18.912   <2e-16 ***
    ## BUSTED$diff_exprsTRUE -0.19123    0.14314  -1.336    0.182    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.046 on 4563 degrees of freedom
    ## Multiple R-squared:  0.000391,   Adjusted R-squared:  0.0001719 
    ## F-statistic: 1.785 on 1 and 4563 DF,  p-value: 0.1816

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="BUSTED LRT(non normal **) ~ Dif_Exp")
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
ggplot(BUSTED, aes(x=diff_exprs, y=LRT)) + 
  geom_violin(aes(fill=diff_exprs),trim=F,scale="width") +
  xlab("Differentially Expressed") +
  ylab("LRT_Busted") +
  ylim(NA,(1.5*max(BUSTED$LRT))) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  draw_label(paste("FALSE",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.45*max(BUSTED$LRT)),size=10) +
  draw_label(paste("TRUE",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.45*max(BUSTED$LRT)),size=10) +
  draw_label(paste("R2=",r2,"; p =",p),x=1.5,y=(1.25*max(BUSTED$LRT)),size=10) +
  guides(fill="none")+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->

#### Non parametric

``` r
kruskal.test(BUSTED$LRT,log(BUSTED$average_tpm))
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  BUSTED$LRT and log(BUSTED$average_tpm)
    ## Kruskal-Wallis chi-squared = 4144.1, df = 4126, p-value = 0.4184

``` r
cor(log(BUSTED$average_tpm),BUSTED$LRT,method = "spearman")
```

    ## [1] -0.002665209

#### SLM

``` r
all_lm<-lm(BUSTED$LRT~log(BUSTED$average_tpm))
AIC(all_lm)
```

    ## [1] 13368.06

``` r
models_tab<-lm_table(all_lm,smr=F,add=T,tab = models_tab,id="BUSTED LRT ~ expression",aic = T)
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = BUSTED$LRT ~ log(BUSTED$average_tpm))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -5.8893 -0.2934 -0.2910 -0.2886 12.3518 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             0.284526   0.035284   8.064 9.36e-16 ***
    ## log(BUSTED$average_tpm) 0.002743   0.011382   0.241     0.81    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.046 on 4563 degrees of freedom
    ## Multiple R-squared:  1.273e-05,  Adjusted R-squared:  -0.0002064 
    ## F-statistic: 0.05807 on 1 and 4563 DF,  p-value: 0.8096

``` r
ggplot(BUSTED,aes(x=log(average_tpm),y=LRT)) + xlab("Average TPM") + ylab("BUSTED Likelihood Ratio") +
  geom_point() + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
  draw_label(paste("n =",nrow(BUSTED),"; R2=",round(m$r.squared,3),"; p =",round(m$coefficients[2,4],3)),x=7.5,y=-10,size=10)+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

``` r
hist(BUSTED$LRT[which(BUSTED$LRT != 0)])
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

``` r
all_lm<-lm(BUSTED$LRT[which(BUSTED$LRT != 0)]~log(BUSTED$average_tpm)[which(BUSTED$LRT != 0)])
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = BUSTED$LRT[which(BUSTED$LRT != 0)] ~ log(BUSTED$average_tpm)[which(BUSTED$LRT != 
    ##     0)])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -6.9410 -1.3388 -0.7155  0.5576 11.3028 
    ## 
    ## Coefficients:
    ##                                                  Estimate Std. Error t value
    ## (Intercept)                                      1.350865   0.132841  10.169
    ## log(BUSTED$average_tpm)[which(BUSTED$LRT != 0)] -0.003232   0.042145  -0.077
    ##                                                 Pr(>|t|)    
    ## (Intercept)                                       <2e-16 ***
    ## log(BUSTED$average_tpm)[which(BUSTED$LRT != 0)]    0.939    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.903 on 992 degrees of freedom
    ## Multiple R-squared:  5.93e-06,   Adjusted R-squared:  -0.001002 
    ## F-statistic: 0.005882 on 1 and 992 DF,  p-value: 0.9389

``` r
ggplot(BUSTED[which(BUSTED$LRT != 0),],aes(x=log(average_tpm),y=LRT)) + xlab("Average TPM") + ylab("BUSTED Likelihood Ratio") +
  geom_point() + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
  draw_label(paste("n =",nrow(BUSTED),"; R2=",round(m$r.squared,3),"; p =",round(m$coefficients[2,4],3)),x=7.5,y=-10,size=10)+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-66-2.png)<!-- -->

### Do BUSTED Likelihood Ratios correlate with average expression levels **in toxins**?

get data

``` r
x<-BUSTED_tox %>% group_by(diff_exprs) %>% summarize(mean=mean(LRT),sd=sd(LRT),n=length(LRT))
kable(x) %>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", font_size = 10), full_width = F)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
diff\_exprs
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FALSE
</td>
<td style="text-align:right;">
0.7354095
</td>
<td style="text-align:right;">
1.5437869
</td>
<td style="text-align:right;">
59
</td>
</tr>
<tr>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:right;">
-0.0000137
</td>
<td style="text-align:right;">
0.0000412
</td>
<td style="text-align:right;">
9
</td>
</tr>
</tbody>
</table>

#### Non parametric

``` r
kruskal.test(BUSTED_tox$LRT,BUSTED_tox$diff_exprs)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  BUSTED_tox$LRT and BUSTED_tox$diff_exprs
    ## Kruskal-Wallis chi-squared = 4.7568, df = 1, p-value = 0.02918

``` r
models_tab<-chi2_table(model = kruskal.test(BUSTED_tox$LRT,BUSTED_tox$diff_exprs),add = T,tab = models_tab,id = "TOXIN BUSTED LRT ~ Dif Exp, Kruskal test")
```

#### SLM

``` r
(m<-summary(lm(BUSTED_tox$LRT~BUSTED_tox$diff_exprs)))
```

    ## 
    ## Call:
    ## lm(formula = BUSTED_tox$LRT ~ BUSTED_tox$diff_exprs)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7376 -0.7354 -0.7320  0.0000  8.3758 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                 0.7354     0.1884   3.903 0.000225 ***
    ## BUSTED_tox$diff_exprsTRUE  -0.7354     0.5179  -1.420 0.160302    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.447 on 66 degrees of freedom
    ## Multiple R-squared:  0.02965,    Adjusted R-squared:  0.01495 
    ## F-statistic: 2.017 on 1 and 66 DF,  p-value: 0.1603

``` r
models_tab<-lm_table(m,smr=T,add=T,tab = models_tab,id="TOXIN BUSTED LRT ~ Dif_Exp")
mean_1<-round(x$mean[1],2)
mean_2<-round(x$mean[2],2)
sd_1<-round(x$sd[1],2)
sd_2<-round(x$sd[2],2)
n_1<-x$n[1]
n_2<-x$n[2]
p<-round(m$coefficients[2,4],2)
r2<-round(m$r.squared,2)
ggplot(BUSTED_tox, aes(x=diff_exprs, y=LRT)) + 
  geom_violin(aes(fill=diff_exprs),trim=F) +
  xlab("Differentially Expressed") +
  ylim(NA,(1.2*max(BUSTED_tox$LRT))) +
  scale_fill_manual(values=c("dodgerblue2","firebrick2")) +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  draw_label(paste("FALSE",paste(mean_1,"±",sd_1),paste("n =",n_1),sep='\n'),x=1,y=(1.45*max(BUSTED_tox$LRT)),size=10) +
  draw_label(paste("TRUE",paste(mean_2,"±",sd_2),paste("n =",n_2),sep='\n'),x=2,y=(1.45*max(BUSTED_tox$LRT)),size=10) +
  draw_label(paste("R2=",r2,"; p =",p),x=1.5,y=(1.25*max(BUSTED_tox$LRT)),size=10) +
  guides(fill="none")+
  theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

``` r
all_lm<-lm(BUSTED_tox$LRT~log(BUSTED_tox$average_tpm))
AIC(all_lm)
```

    ## [1] 248.9039

``` r
models_tab<-lm_table(all_lm,smr=F,add=T,tab = models_tab,id="TOXIN BUSTED LRT ~ expression",aic = T)
(m<-summary(all_lm))
```

    ## 
    ## Call:
    ## lm(formula = BUSTED_tox$LRT ~ log(BUSTED_tox$average_tpm))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.8210 -0.6553 -0.5762  0.0048  8.4325 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                 -0.01810    1.12545  -0.016    0.987
    ## log(BUSTED_tox$average_tpm)  0.08276    0.14017   0.590    0.557
    ## 
    ## Residual standard error: 1.465 on 66 degrees of freedom
    ## Multiple R-squared:  0.005254,   Adjusted R-squared:  -0.009818 
    ## F-statistic: 0.3486 on 1 and 66 DF,  p-value: 0.5569

``` r
print("AIC BUSTED ~ Expression")
```

    ## [1] "AIC BUSTED ~ Expression"

``` r
AIC(all_lm)
```

    ## [1] 248.9039

#### SLM

``` r
all_lm1<-lm(BUSTED_tox$LRT~log(BUSTED_tox$average_tpm)+BUSTED_tox$toxin_family)
models_tab<-lm_table(all_lm1,tf=T,smr=F,add=T,tab = models_tab,id="TOXINS BUSTED LRT ~ expression + toxin family",aic = T)
(m1<-summary(all_lm1))
```

    ## 
    ## Call:
    ## lm(formula = BUSTED_tox$LRT ~ log(BUSTED_tox$average_tpm) + BUSTED_tox$toxin_family)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2872 -0.6787 -0.2418  0.0000  8.2078 
    ## 
    ## Coefficients:
    ##                                        Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                            -0.83762    2.07168  -0.404   0.6877  
    ## log(BUSTED_tox$average_tpm)             0.11875    0.20393   0.582   0.5629  
    ## BUSTED_tox$toxin_familyCRISP            1.31803    1.76178   0.748   0.4578  
    ## BUSTED_tox$toxin_familyCTL              0.65122    1.06076   0.614   0.5420  
    ## BUSTED_tox$toxin_familyHYAL             3.52915    1.92280   1.835   0.0723 .
    ## BUSTED_tox$toxin_familyKUN              0.66665    1.93698   0.344   0.7321  
    ## BUSTED_tox$toxin_familyLAAO            -0.27023    1.76139  -0.153   0.8787  
    ## BUSTED_tox$toxin_familyNUC              0.02263    1.82408   0.012   0.9902  
    ## BUSTED_tox$toxin_familyPDE              0.13870    1.88591   0.074   0.9417  
    ## BUSTED_tox$toxin_familyPLA2             0.42073    1.07666   0.391   0.6976  
    ## BUSTED_tox$toxin_familyPLB              0.44024    1.47955   0.298   0.7673  
    ## BUSTED_tox$toxin_familySVMPI            2.15257    1.24580   1.728   0.0901 .
    ## BUSTED_tox$toxin_familySVMPII           0.38783    1.07897   0.359   0.7207  
    ## BUSTED_tox$toxin_familySVMPIII          0.74111    1.01392   0.731   0.4682  
    ## BUSTED_tox$toxin_familySVSP             0.16349    1.00064   0.163   0.8709  
    ## BUSTED_tox$toxin_familyuncharacterised  0.02265    1.82409   0.012   0.9901  
    ## BUSTED_tox$toxin_familyVEGF            -0.20810    1.76275  -0.118   0.9065  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.525 on 51 degrees of freedom
    ## Multiple R-squared:  0.1674, Adjusted R-squared:  -0.09386 
    ## F-statistic: 0.6407 on 16 and 51 DF,  p-value: 0.8352

``` r
print("AIC BUSTED ~ Expression")
```

    ## [1] "AIC BUSTED ~ Expression"

``` r
AIC(all_lm)
```

    ## [1] 248.9039

### Figure 3C

``` r
(testplot<-ggplot(BUSTED_tox,aes(x=log(average_tpm),y=LRT, color = toxin_family)) + xlab("Average TPM") + ylab("LRT Busted model") +
    geom_point() + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
    draw_label(paste("n =",nrow(BUSTED_tox)," ; R2=",round(m$r.squared,2),"; p =",round(m$coefficients[2,4],2)),x=7.5,y=3,size=10))+
    theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

``` r
(Fig3C<-ggplot(BUSTED_tox,aes(x=log(average_tpm),y=LRT)) + xlab("Average TPM") + ylab("LRT Busted model") +
    geom_point() + geom_abline(intercept = coef(all_lm)[1],slope=coef(all_lm)[2],lty=2) +
    draw_label(paste("n =",nrow(BUSTED_tox),"; R2=",round(m$r.squared,2),"; p =",round(m$coefficients[2,4],2)),x=7.5,y=3,size=10))+
    theme_classic()
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

## COWPLOTS

``` r
Fig2A<-Fig2A+theme_classic()
Fig2B<-Fig2B+theme_classic()
Fig2C<-Fig2C+theme_classic()
Fig3A<-Fig3A+theme_classic()
Fig3B<-Fig3B+theme_classic()
Fig3C<-Fig3C+theme_classic()

plot_grid(Fig2A,Fig2B,Fig2C,Fig3A,Fig3B,Fig3C,align="hv",nrow=2,ncol=3,labels="AUTO")
```

![](Cgodm_selection_BMC_files/figure-gfm/unnamed-chunk-73-1.png)<!-- -->

``` r
output<-plot_grid(Fig2A,Fig2B,Fig2C,Fig3A,Fig3B,Fig3C,align="hv",nrow=2,ncol=3,labels="AUTO")

ggsave("Cgodm_selection_plots.svg",output,width = 30, height = 20, units = "cm")
#pdf("Cgodm_selection_plots.pdf",width = 16, height=10)
#plot_grid(Fig2A,Fig2B,Fig2C,Fig3A,Fig3B,Fig3C,align="hv",nrow=2,ncol=3,labels="AUTO")
#dev.off() 

colnames(models_tab)[2]<-"Statistic"
write.csv(models_tab,"models_summary.csv")
```
