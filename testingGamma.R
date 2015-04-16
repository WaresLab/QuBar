
# set parameters to plot, we should have the number of haplotypes
# know the real number of samples
# and have haplotype diversity from the original field population

numhaps<-c(1,2,3,4,5,6,7,8)
numSamp<-c(3,6,9,10,12,16,18,19)
Hapdiv=0.7 #minvalue 0.001 if no information; this is PRIOR information

# set plot space depending number of datasets
par(mfrow=c(2,4))

# loop through all the haplotypes and generate plots
for(i in 1:length(numhaps)){
x=1
cdf=0
indprob=0
array<-NULL
while (cdf<0.99) {
  cdfprev<-cdf
  #  cdf<-pgamma(x,1,Hapdiv) 
  #  if use 1 as shape parameter keeping shape parameter constant doesn't account for increased variance (?) as numhaps go up, 
  #  e.g. error may be higher as you observe more...once it is working run it by somebody mathier.
  cdf<-pgamma(x,numhaps[i],Hapdiv) #might be that numhaps is actually the shape parameter!!!! or: something else...non-gamma.
  indprob<-cdf-cdfprev
    
  happrob<-numhaps[i]+(x-1)
  array<-c(array,happrob)
  array<-c(array,cdf)
  array<-c(array,indprob)
  #  print(happrob)
  #  print(cdf)
  x=x+1
   
  }
probs<-t(matrix(array,nrow=3))
probs
title<-paste(c("n=",numSamp[i]))
plot(probs[,1],probs[,3],col='red',main=title)
}
