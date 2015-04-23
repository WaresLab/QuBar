
# set parameters to plot, we should have the number of haplotypes
# know the real number of samples
# and have haplotype diversity from the original field population

# Population 1 (theta=2, no growth)
numhaps<-c(1,2,3,4,6,6,9) # I will have this vector for real
numsamp<-c(2,4,8,16,32,64,128) #niter in MS
hapdiv<-c(0.5504204, 0.3377538, 0.8915536, 0.9258679) #minvalue 0.001 if no information; this is PRIOR information

# set plot space depending number of datasets
par(mfrow=c(2,4),oma=c(0,0,2,0))

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
  cdf<-pgamma(x,numhaps[i],hapdiv[4]) #might be that numhaps is actually the shape parameter!!!! or: something else...non-gamma.
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
title<-paste(c("n=",numsamp[i],sep=""))
plot(probs[,1],probs[,3],col='red',main=title)
mtext("Pop4:Theta 10-no growth",outer=T)
}
