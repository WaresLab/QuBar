#theta=2

bigresult<-rep(NA,30)
for (j in 2:30){
  result<-rep(NA,30)
  for (i in 2:30) {
    PSK<-((-1)^i)*(choose((j-1),(i-1)))*((i-1)/(2+i-1))*(2/(2+i-1))^1
    PSK->result[i]
    }
  sum(result,na.rm=T)->ourpsk
  bigresult[j]<-ourpsk
}

# theta 10

theta=10
maxn=50
obsvdk=10

bigresult<-rep(NA,maxn)
for (j in 2:maxn){
  result<-rep(NA,maxn)
  for (i in 2:maxn) {
    PSK<-((-1)^i)*(choose((j-1),(i-1)))*((i-1)/(theta+i-1))*(theta/(theta+i-1))^obsvdk
    PSK->result[i]
  }
  sum(result,na.rm=T)->ourpsk
  bigresult[j]<-ourpsk
}



# comparing with excel
maxn=100
theta=2
k=4
i=60

PSK<-((-1)^i)*(choose((maxn-1),(i-1)))*((i-1)/(theta+i-1))*((theta/(theta+i-1))^k)