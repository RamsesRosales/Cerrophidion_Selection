lm_table<-function(model,tf=F,smr=F,add=F,tab,id="NO",aic=F){
  d<-"NA"
  if(aic){
    d<- as.numeric(AIC(model))
  }
  if(!smr){
    model<-summary(model)
  }
  if(!tf){
    a = strsplit(as.character(model$call),"=")[2]
    b = model$r.squared
    c = model$coefficients[2,4]
  }
  if(tf){
    a = as.character(strsplit(as.character(model$call),"=")[2])
    b = as.numeric(model$r.squared)
    c = as.numeric(pf(model$fstatistic[1],model$fstatistic[2],model$fstatistic[3],lower.tail = F))
  }
  output<- data.frame("formula"=a,"R.squared"=b,"p-value"=c,"id" = id,aic=d)
  if(!add){
    return(output)
  }
  if(add){
    tab[(nrow(tab)+1),]<- output[1,]
    return(tab)
  }
}

chi2_table<-function(model,add=F,tab,id="NO"){
  a = as.character(model$data.name)
  b = as.numeric(model$statistic)
  c = as.numeric(model$p.value)
  output<- data.frame("formula"= a,"Chisquared"= b,"p-value"= c,"id" = id,aic = "NA")
  if(!add){
    return(output)
  }
  if(add){
    tab[(nrow(tab)+1),]<- output[1,]
    return(tab)
  }
}


t_table<-function(model,add=F,tab,id="NO"){
  a = as.character(model$data.name)
  b = as.numeric(model$statistic)
  c = as.numeric(model$p.value)
  d = paste0("sd = ",model$stderr)
  output<- data.frame("formula"= a,"t.statistic"= b,"p-value"= c,"id" = id,sd = d)
  if(!add){
    return(output)
  }
  if(add){
    tab[(nrow(tab)+1),]<- output[1,]
    return(tab)
  }
}
f_table<-function(model,add=F,tab,id="NO"){
  a = as.character(model$data.name)
  b = as.numeric(model$statistic)
  c = as.numeric(model$p.value)
  output<- data.frame("formula"= a,"t.statistic"= b,"p-value"= c,"id" = id,aic = "NA")
  if(!add){
    return(output)
  }
  if(add){
    tab[(nrow(tab)+1),]<- output[1,]
    return(tab)
  }
}

#kk_table<-function(model,add=F,tab,id="NO"){
#  a = as.character(model$data.name)
#  b = as.numeric(model$statistic)
#  c = as.numeric(model$p.value)
#  output<- data.frame("formula"= a,"Chisquared"= b,"p-value"= c,"id" = id,aic = "NA")
#  if(!add){
#    return(output)
#  }
#  if(add){
#    tab[(nrow(tab)+1),]<- output[1,]
#    return(tab)
#  }
#}


