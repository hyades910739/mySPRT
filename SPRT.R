##mySPRT
## version 0.87, 2017.4.28. 04:15 A.M.
## onlt Bernoulli test now
## author: Hyades Lai

library(magrittr)


########## test data 1########## 
set.seed(689)
test2 = replicate(20,rbinom(1,1,0.8))
sprt2 = mySPRT(test2,0.5,0.8,0.05,0.2)
print(sprt2)
plot(sprt2)

########## test data 2########## 
set.seed(689)
test3 = replicate(20,rbinom(1,1,0.5))
sprt3 = mySPRT(test3,0.35,0.55,0.05,0.1)
print(sprt3)
plot(sprt3)

###############################

mySPRT = function(data,H0,H1,alpha,beta,distribution="bernoulli"){
  
  n = length(data);
  likelihood = lh.f(data,H0,H1,distribution);
  boundary = boundary.cal(alpha,beta,log=T);
  
  
  ##fake SPRT start:
  anyReject = which(likelihood>boundary$A)
  anyAccept = which(likelihood<boundary$B)
  if(length(anyReject)==0 & length(anyAccept)==0) decision = "Continue Trials!"    
  if(length(anyReject)!=0 & length(anyAccept)==0) decision = "Reject H0!"
  if(length(anyReject)==0 & length(anyAccept)!=0) decision = "Do not reject H0!"
  if(length(anyReject)!=0 & length(anyAccept)!=0){
    if( min(anyReject) < min(anyAccept) ){
      decision = "Reject H0!"
    }else{
      decision = "Do not reject H0!"
    }
  } 
  
  ##prepare for result:
  res = list(
    decision = decision,
    llr =  exp(likelihood),
    log.llr = likelihood,
    boundary.A = boundary$A,
    boundary.B = boundary$B,
    H0 = H0,
    H1 = H1,
    alpha = alpha,
    beta = beta,
    n=n,
    distribution = distribution
  )
  class(res) = "mySPRT"
  
  return(res)
}


# return vector of liklihood ration
lh.f = function(data,H0,H1,distribution="bernoulli"){
  n = length(data)
  if(distribution =="bernoulli"){
    lhr = sapply(1:n,function(x){
      f1 = data[1:x] %>% dbinom(.,1,H1,log=T) %>% sum()
      f2 = data[1:x] %>% dbinom(.,1,H0,log=T) %>% sum()
      (f1-f2)
    })
    return(lhr) 
  }  
}

# return list of boundary.A boundary.B
boundary.cal = function(alpha,beta,log=T){
  A = (1-beta)/alpha;
  B = beta/(1-alpha);
  if(log){
    A = log(A);
    B = log(B);
  }
  boundary = list(A=A,B=B);
  return(boundary)
}

#S3 methods:
plot.mySPRT = function(sprt){

  A = sprt$boundary.A %>% exp()
  B = sprt$boundary.B %>% exp()
  A.text = paste("boundary.A: ",A,sep="")
  B.text = paste("boundary.B: ",round(B,4),sep="")
  
  plot(1:sprt$n,sprt$llr,pch=19,lwd=3,xlab = "Trials",ylab="Likelihood-Ratio",
       ylim = c(min(B,sprt$llr),max(A,sprt$llr)))
  abline(A,0,lwd=1,lty=2,col="red")
  abline(B,0,lwd=1,lty=2,col="blue")
  par(xpd=TRUE)
  legend(x=(sprt$n/4),y=max(A,sprt$llr)*1.2,c(A.text,B.text),lty=2,col=c("red","blue"))
  par(xpd=F)
}

print.mySPRT = function(sprt){
 cat("#SPRT test:\n",
     "decision:",sprt$decision,"\n\n",
     "---\n",
     "distribution = ",sprt$distribution,"\n",
     "n = ",sprt$n,"\n",
     "boundary.A = ",sprt$boundary.A,"\n",
     "boundary.B = ",sprt$boundary.B,"\n",
     "alpha = ",sprt$alpha,"\n",
     "beta = ",sprt$beta,"\n") 
  
}
