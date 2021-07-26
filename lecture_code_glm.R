d=read.table("LogisticRegression.txt",header=T,sep="\t")
str(d)
d$smoking=factor(d$smoking, levels=c("Nonsmoker", "Smoker"))
d$gender=as.factor(d$gender)
d$impaired=as.factor(d$impaired)
table(d$smoking,d$impaired,d$gender)
fit.1=glm(impaired~smoking, family="binomial",data=d)
fit.1$coefficients
fit.1$coefficients[1]
fit.1$coefficients[2]
exp(fit.1$coefficients)
confint(fit.1)
inv.logit=function(x){
	return(exp(x)/(1+exp(x)))
}
inv.logit(fit.1$coefficients)

library(oddsratio)
or_glm(data=d, model=fit.1)

summary(fit.1)
summary(fit.1)$coefficients

fit.2=glm(impaired~gender+smoking, family="binomial",data=d)
summary(fit.2)
odds.ratio(fit.2, gender)
or_glm(data=d, model=fit.2)

fit.3=glm(impaired~smoking, family="binomial",data=d, subset= gender=="Male")

##Improvement test
glm.0 <- glm(impaired~1, family="binomial",data=d)
anova(fit.1, glm.0, test="F")

library(lmtest)
lrtest(fit.1, glm.0)

##goodness-of-fit test
library(LogisticDx)
g1=gof(fit.1)
g1
unclass(g1)
dev_ns=ob_ns*log(ob_ns/exp_ns)*2

##Poisson
n=40
##generate data
X=rnorm(n, mean=0.5, sd=1)
lam=exp(0.2+0.5*X)
Y=rpois(n,lam)
plot(Y~X)

##fit model
fit.p=glm(Y~X, family="poisson")
summary(fit.p)

##Deviance test
with(fit.p, cbind(res.deviance = deviance, df = df.residual,
  p = pchisq(deviance, df.residual, lower.tail=FALSE)))

##improvement test
fit.p0=update(fit.p, .~.-X)
anova(fit.p0, fit.p, test="Chisq")

##NB model
Y=rnbinom(n, mu=lam, size=2)
library(MASS)
fit.nb=glm.nb(Y~X)
fit.nb
summary(fit.nb)

##odTest
library(pscl)
odTest(fit.nb)
