##mySPRT
## version 2.87, 2017.5.5. 10:35 A.M. 
## only Bernoulli test now
## add : argument: previous. Use previous to sent previous llr to the new test.
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
## decision:  continue trials!
newdata = c(0,1,0,1,1,1,0,1,0,1);
sprt4 = mySPRT(newdata,0.35,0.55,0.05,0.1,previous = sprt3)
plot(sprt4)
###############################

##########testSPRT############

testSPRT = function(repeatT,n,actual,h0,h1,alpha,beta){
  test = replicate(repeatT,replicate(n,rbinom(1,1,actual)))
  
  apply(test,2,FUN=function(x){
    mySPRT(data=x,H0=h0,H1=h1,alpha=alpha,beta=beta)$decision
  }) %>% table() %>% '/'(repeatT)
}

testSPRT2 = function(repeatT,n,actual,h0,h1,alpha,beta){
  test = replicate(repeatT,replicate(n,rbinom(1,1,actual)))
  header = paste("sample size \nalpha=",alpha,",beta=",beta,",repat=",repeatT,sep="")
  
  apply(test,2,FUN=function(x){
    mySPRT(data=x,H0=h0,H1=h1,alpha=alpha,beta=beta)$n
  })  %T>% boxplot(.,main=header) %>% summary()
}

testSPRT2(1000,1000,0.5,0.5,0.8,0.05,0.1)
testSPRT2(1000,1000,0.5,0.5,0.8,0.1,0.1)
testSPRT2(1000,1000,0.5,0.5,0.6,0.05,0.1)


set.seed(689)
##test for type I error /w repeat 100, 1000 times:
testSPRT(100,1000,0.5,0.5,0.8,0.05,0.1)
testSPRT(1000,1000,0.5,0.5,0.8,0.05,0.1)

##test for type II error /w repeat 100, 1000 times:
testSPRT(100,1000,0.5,0.8,0.5,0.05,0.1)
testSPRT(1000,1000,0.5,0.8,0.5,0.05,0.1)


############for test end ##############

#main function:
mySPRT = function(data,H0,H1,alpha,beta,distribution="bernoulli",previous=NULL){
  boundary = boundary.cal(alpha,beta,log=T);  
####for llr
  if( class(previous)=="mySPRT"){
    if(previous$decision != "Continue Trials!") return(previous);
    if(previous$distribution != distribution) retrun("fause! Different distribution been applied.");
    prelogllr = previous$log.llr[length(previous$log.llr)];
    llh2 =  lh.f2.pre(data,H0,H1,boundary,distribution="bernoulli",pre.log.llr =prelogllr );
    llh2 = c(previous$log.llr,llh2);
  }else{
    llh2 =  lh.f2(data,H0,H1,boundary,distribution="bernoulli");     
  }
  n2 = length(llh2);  
####decisions:
  if(llh2[length(llh2)]>boundary$A) decision = "Reject H0!";
  if(llh2[length(llh2)]<boundary$B) decision = "Do not reject H0!";
  if(llh2[length(llh2)]<boundary$A & llh2[n2]>boundary$B ) decision = "Continue Trials!";

####prepare for result:
  res = list(
    decision = decision,
    llr = exp(llh2),
    log.llr = llh2,
    boundary.A = boundary$A,
    boundary.B = boundary$B,
    H0 = H0,
    H1 = H1,
    alpha = alpha,
    beta = beta,
    n=n2,
    distribution = distribution
  )
  class(res) = "mySPRT";
  return(res);
}


#return vector of liklihood ration
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

# calculate likelihood ratio:
lh.f2 = function(data,H0,H1,boundary,distribution="bernoulli"){
  n = length(data);
  f1 = c();
  f2 = c();
  lh = c();

  if(distribution =="bernoulli"){
    for(i in 1:n){
      f1[i] = data[1:i] %>% dbinom(. ,1,H1,log=T) %>% sum();
      f2[i] = data[1:i] %>% dbinom(. ,1,H0,log=T) %>% sum();
      lh[i] = f1[i]-f2[i];
      if(lh[i]>boundary$A){
        return (lh)
      }else if(lh[i]<boundary$B){
        return (lh)
      }
    } 
    return(lh)
  }  
}


# calculate likelihood ratio if previous exist:
lh.f2.pre = function(data,H0,H1,boundary,distribution="bernoulli",pre.log.llr){
  n = length(data);
  f1 = c();
  f2 = c();
  lh = c();
  
  if(distribution =="bernoulli"){
    for(i in 1:n){
      f1[i] = data[1:i] %>% dbinom(. ,1,H1,log=T) %>% sum();
      f2[i] = data[1:i] %>% dbinom(. ,1,H0,log=T) %>% sum();
      lh[i] = f1[i]-f2[i]+pre.log.llr;
      if(lh[i]>boundary$A){
        return (lh)
      }else if(lh[i]<boundary$B){
        return (lh)
      } 
    }  
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


############ S3 methods: ############
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
