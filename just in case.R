factor(miniwak$segSites)->miniwak$segSites
keep<-subset(miniwak,miniwak$n==2|miniwak$n==4|miniwak$n==8|miniwak$n==16|miniwak$n==32)
row.names(keep)<-NULL

# we need to see which n has the maximum probablity for each segregating site according Wakeley's formula
sumdata<-ddply(keep, c("Theta","segSites"), summarise, maxprob = max(prob))

# add flag to the bests n's to keep track later
sumdata$best.n<-"best"
sumdata$Theta<-NULL; sumdata$segSites<-NULL

# so we merge with Wakeley's data to flag the best n
merge(miniwak,sumdata,by.x="prob",by.y="maxprob",all.x=T)->newdata

# keep only the flag, that is the best.n data
goodwak<-newdata[complete.cases(newdata),]
row.names(goodwak)<-NULL
goodwak$best.n<-NULL
names(goodwak)<-c("prob","Theta","best.n","segSites")
# goodwak is the summary of Wakeley for the higher probilities
