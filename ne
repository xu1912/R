se=0.9
sp=0.8
spi=1-sp
dbinom(5,5,se) + dbinom(4,5,se)*pbinom(1,95,spi,lower.tail=F)+ dbinom(3,5,se)*dbinom(2,95,spi) + dbinom(2,5,se)*dbinom(3,95,spi) + 
dbinom(1,5,se)*dbinom(4,95,spi) + dbinom(0,5,se)*dbinom(5,95,spi)

se=0.9
sp=0.999
spi=1-sp
N=10000
K=N*0.05+1
P=se^K
rp=c()
for(i in 1:K){
	P=dbinom(i-1,K,se)*pbinom(K-i+1,N-K,spi,lower.tail=F)+P
	rp[i]=dbinom(i-1,K,se)*pbinom(K-i+1,N-K,spi,lower.tail=F)
}
P
