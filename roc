library(ggplot2)

pth=c(12,48,29,24,7.4,16.3,37.3,2.3,46,5,28,16,10.3,8.8,40.8,37.9,14.0,12)
ca=c(1.01,1.09,1.1,1.1,1.04,1.15,1.13,1.07,1.07,1.06,1.17,1.06,1.05,1.05,1.07,1.05,1.37,1.1)
cor.test(pth,ca,method="spearman")

d=data.frame(pth,ca)
ggplot(d, aes(x=pth, y=ca)) +
  geom_point() +
  geom_smooth(method = lm) +
  labs(x = "PTH(2H)", y = "Ionized Ca(6-8H)") +
  theme_classic()+
  theme(text = element_text(size=20))


d$hypocalcemia=0
d$hypocalcemia[d$ca<1.05]=1

library(pROC)
fit=roc(d, hypocalcemia, pth)
plot(fit,print.auc=TRUE)

coords(fit, "best")

library(OptimalCutpoints)
optimal.cutpoint.Youden<-optimal.cutpoints(X = "pth", status = "hypocalcemia", tag.healthy = "0", data=d, direction=">",
methods = "Youden",  pop.prev = NULL, control = control.cutpoints(), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
plot(optimal.cutpoint.Youden)

d$hypocalcemia_pred=0
d$hypocalcemia_pred[d$pth<=12]=1

library(caret)
confusionMatrix(as.factor(d$hypocalcemia_pred), as.factor(d$hypocalcemia), positive="1")
