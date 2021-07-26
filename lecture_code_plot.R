##Mathematical annotation
x<-seq(-10,10,length=400)
y1<-dnorm(x)
y2<-dnorm(x,m=3)
plot(x,y2,xlim=c(-3,8),type="l",
xlab=expression(Z==frac(mu[1]-mu[2],sigma/sqrt(n))),
ylab="Density")
lines(x,y1)

polygon(c(1.96,1.96,x[240:400],10),
c(0,dnorm(1.96,m=3),y2[240:400],0),
col="grey80",lty=0)

polygon(c(-1.96,-1.96,x[161:1],-10),
c(0,dnorm(-1.96,m=0),y1[161:1],0),
col="grey30",lty=0)

polygon(c(1.96,1.96,x[240:400],10),
c(0,dnorm(1.96,m=0),y1[240:400],0),
col="grey30")

legend(4.2,.4,fill=c("grey80","grey30"),
legend=expression(P(abs(Z)>1.96,H[1])==0.85,
P(abs(Z)>1.96,H[0])==0.05),bty="n")
text(0,.2,quote(H[0]:~~mu[1]==mu[2]))
text(3,.2,quote(H[1]:~~mu[1]==mu[2]+delta))

##ggplot
library(ggplot2)

ggplot(mpg, aes(displ, hwy)) + 
  geom_point(size=3, color="blue") +
  geom_smooth(method = lm, se=F, color="red")

ggplot(mpg, aes(displ, hwy, color = class)) + 
  geom_point(size=3) +
  geom_smooth(method = lm, se=F) +
  scale_color_brewer(palette="Dark2")

ggplot(mpg, aes(displ, hwy, color = class)) + 
  geom_point() +
  geom_smooth(method = lm, se=F) +
  labs(x = "Engine displacement (litres)", y = "Highway miles per gallon", colour="Type of car") +
  theme_classic()+
  theme(axis.text = element_text(size=14), legend.position=c(0.8,0.8), legend.background=element_blank(),
		legend.text= element_text(size=14), axis.title = element_text(size=20))

ggplot(mpg, aes(displ, hwy)) + 
  geom_point() +
  geom_smooth(method = lm)

ggplot(mpg, aes(displ, hwy, colour = class)) + 
  geom_boxplot()

ggplot(mpg, aes(displ, hwy, fill = class)) + 
  geom_boxplot()

ggplot(mpg, aes(displ, hwy, colour = class, fill = class)) + 
  geom_boxplot() +
  scale_fill_brewer(palette="Dark2")

##Multiplot
source("X://R Course (sponsored by BERD) 2020//Part II//Topic8//multiplot.r")

p1=ggplot(mpg, aes(displ, hwy, colour = class)) + 
  geom_point() +
  geom_smooth(method = lm)

p2=ggplot(mpg, aes(displ, hwy, colour = class)) + 
  geom_boxplot()

multiplot(p1, p2, cols=1)

##Graph quality
windows(height=4, width=6)
x=rnorm(10)
y=5*x+rnorm(10)

pdf("example.pdf")
plot(y~x)
dev.off()

png("example.png", width=4, height=4, units='in', res=300)
plot(y~x)
dev.off()

png("example.png", width=4, height=4, units='in', res=100)
plot(y~x)
dev.off()
