# load libraries
library(lattice)

# set parameters to plot, we should have the number of haplotypes
# know the real number of samples
# and have haplotype diversity from the original field population

# Number of haplotypes for each simulation, these number were originated in the script "runningSimulations"
numhapsPop1<-c(1,2,3,4,6,6,9) 
numhapsPop2<-c(1,1,3,5,4,8,10) 
numhapsPop3<-c(2,4,5,8,15,17,26) 
numhapsPop4<-c(2,3,5,9,16,21,33) 

# Number of samples, same for all populations
numsamp<-c(2,4,8,16,32,64,128) 

# Haplotype diversity for each of the four populations
hapdiv<-c(0.5504204, 0.3377538, 0.8915536, 0.9258679) #minvalue 0.001 if no information; this is PRIOR information

# We can make this more efficient, but for now, let's keep a loop for each population:

# set dataframe to fill with the predicted and observed values,
# for each theta and each sampling size
c(rep("Pop1",7),rep("Pop2",7),rep("Pop3",7),rep("Pop4",7))->col1
c(rep("two",14),rep("ten",14))->col2
c(rep("growth",7),rep("nogrowth",7),rep("growth",7),rep("nogrowth",7))->col3
cbind(col1,col2,col3)->mycols
c("Population","theta","growth")->namesCols
as.data.frame(mycols)->toFill
names(toFill)<-namesCols

# set plot space depending number of datasets, for each theta we have 7 plots
par(mfrow=c(2,4),oma=c(0,0,2,0)) #oma is to set space for the main tittle

# Loop population 1
for(i in 1:length(numhapsPop1)){
    x=1
    cdf=0
    indprob=0
    array<-NULL
  while (cdf<0.99) {
    cdfprev<-cdf
    #  cdf<-pgamma(x,1,Hapdiv) 
    #  if use 1 as shape parameter keeping shape parameter constant doesn't account for increased variance (?) as numhaps go up, 
    #  e.g. error may be higher as you observe more...once it is working run it by somebody mathier.
    cdf<-pgamma(x,numhapsPop1[i],hapdiv[1]) #might be that numhaps is actually the shape parameter!!!! or: something else...non-gamma.
    indprob<-cdf-cdfprev
      
    happrob<-numhapsPop1[i]+(x-1)
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
  max<-max(probs[,3])
  maxPred<-probs[which(probs[,3]==max),1]
  toFill$Pred.value[i]<-maxPred
  toFill$Obs.n[i]<-numsamp[i]
  plot(probs[,1],probs[,3],pch=19,col='red',main=title)
  mtext("Pop1:Theta 2-growth",outer=T)
  dev.copy(pdf, "Pop1_theta2_growth.pdf", width=14, height=8)
  dev.off()
}

# Loop population 2
par(mfrow=c(2,4),oma=c(0,0,2,0))
for(i in 1:length(numhapsPop2)){
  x=1
  cdf=0
  indprob=0
  array<-NULL
  while (cdf<0.99) {
    cdfprev<-cdf
    #  cdf<-pgamma(x,1,Hapdiv) 
    #  if use 1 as shape parameter keeping shape parameter constant doesn't account for increased variance (?) as numhaps go up, 
    #  e.g. error may be higher as you observe more...once it is working run it by somebody mathier.
    cdf<-pgamma(x,numhapsPop2[i],hapdiv[2]) #might be that numhaps is actually the shape parameter!!!! or: something else...non-gamma.
    indprob<-cdf-cdfprev
    
    happrob<-numhapsPop2[i]+(x-1)
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
  max<-max(probs[,3])
  maxPred<-probs[which(probs[,3]==max),1]
  toFill$Pred.value[i+7]<-maxPred
  toFill$Obs.n[i+7]<-numsamp[i]
  plot(probs[,1],probs[,3],col='red',pch=19,main=title)
  mtext("Pop2:Theta 2-no growth",outer=T)
  dev.copy(pdf, "Pop2_theta2_nogrowth.pdf", width=14, height=8)
  dev.off()
}

# Loop population 3
par(mfrow=c(2,4),oma=c(0,0,2,0))
for(i in 1:length(numhapsPop3)){
  x=1
  cdf=0
  indprob=0
  array<-NULL
  while (cdf<0.99) {
    cdfprev<-cdf
    #  cdf<-pgamma(x,1,Hapdiv) 
    #  if use 1 as shape parameter keeping shape parameter constant doesn't account for increased variance (?) as numhaps go up, 
    #  e.g. error may be higher as you observe more...once it is working run it by somebody mathier.
    cdf<-pgamma(x,numhapsPop3[i],hapdiv[3]) #might be that numhaps is actually the shape parameter!!!! or: something else...non-gamma.
    indprob<-cdf-cdfprev
    
    happrob<-numhapsPop3[i]+(x-1)
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
  max<-max(probs[,3])
  maxPred<-probs[which(probs[,3]==max),1]
  toFill$Pred.value[i+14]<-maxPred
  toFill$Obs.n[i+14]<-numsamp[i]
  plot(probs[,1],probs[,3],col='red',pch=19,main=title)
  mtext("Pop3:Theta 10-growth",outer=T)
  dev.copy(pdf, "Pop3_theta10_growth.pdf", width=14, height=8)
  dev.off()
}

# Loop population 4
par(mfrow=c(2,4),oma=c(0,0,2,0))
for(i in 1:length(numhapsPop4)){
  x=1
  cdf=0
  indprob=0
  array<-NULL
  while (cdf<0.99) {
    cdfprev<-cdf
    #  cdf<-pgamma(x,1,Hapdiv) 
    #  if use 1 as shape parameter keeping shape parameter constant doesn't account for increased variance (?) as numhaps go up, 
    #  e.g. error may be higher as you observe more...once it is working run it by somebody mathier.
    cdf<-pgamma(x,numhapsPop4[i],hapdiv[4]) #might be that numhaps is actually the shape parameter!!!! or: something else...non-gamma.
    indprob<-cdf-cdfprev
    
    happrob<-numhapsPop4[i]+(x-1)
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
  max<-max(probs[,3])
  maxPred<-probs[which(probs[,3]==max),1]
  toFill$Pred.value[i+21]<-maxPred
  toFill$Obs.n[i+21]<-numsamp[i]
    plot(probs[,1],probs[,3],col='red',pch=19,main=title)
  mtext("Pop4:Theta 10-no growth",outer=T)
  dev.copy(pdf, "Pop4_theta10_nogrowth.pdf", width=14, height=8)
  dev.off()
}

# now we have the dataframe "toFill" with all the values we need to plot
xyplot(Pred.value~Obs.n|theta*growth,data=toFill)
xyplot(Pred.value~Obs.n|theta*growth,data=toFill,panel=function(x,y,...){
  panel.xyplot(x,y,...)
  panel.lines(x,x)})