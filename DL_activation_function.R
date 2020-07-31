library(ggplot2)
x=seq(-10,10,by=0.05)
alpha=1.67326324
lambda=1.05070098
y=x
for(i in 1:length(x)){
if(x[i]>0){
	y[i]=lambda*x[i]
}else{
	y[i]=lambda*(alpha*exp(x[i])-alpha)
}
}
df1=data.frame(x,y,gp="SELU")

for(i in 1:length(x)){
if(x[i]>0){
	y[i]=x[i]
}else{
	y[i]=0
}
}
df2=data.frame(x,y,gp="ReLU")

df=rbind(df1,df2)
ggplot(data=df,aes(x=x,y=y, group=gp))+
geom_line(aes(color=gp), size=1.25)+
scale_color_brewer(name="activation function",palette="Dark2")+
theme_classic()+
theme(text = element_text(size=20), legend.position = c(0.22, 0.9))
