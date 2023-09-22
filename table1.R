library(table1)
library(ggplot2)

d=read.table("demographic_data.txt", header=T, sep="\t", na.strings = "")
##Recode data
d$PageTime=d$Time.since.last.known.normal..mins. / 60
d$Stroke="Y"
d$Stroke[is.na(d$Type.of.Stroke)]="N"

d$PageHour=ceiling(d$PageTime)
d$PageHour[d$PageHour==0]=1
d$PageHour[d$PageHour>8]="Unknown"
d$PageHour=factor(d$PageHour, levels=c("1","2","3","4","5","6","7","8","Unknown"), labels=c("1","2","3","4","5","6","7","8","Unknown"))

d$AIS.Training=factor(d$AIS.Training..1.yes..0.no., levels=c(0,1), labels=c("N", "Y"))
d$AIS.Training_pvalue=factor(d$AIS.Training..1.yes..0.no., levels=c(0,1,2), labels=c("N", "Y", "P-value"))
d$mRS=factor(d$mRS.baseline, levels=c(0,1,2,3,4,5), labels=c("0", "1", "2", "3", "4", "5"))

##change the display name of the variabe in table
label(d$PageHour)="Page Hour"
label(d$AIS.Training)="AIS Training"
label(d$AIS.Training_pvalue)="AIS Training"

##rndr defines the testing method for each variable by d$Group_pvalue. Notice d is the data frame. In other words, the dataframe name must be d.
##d$Group_pvalue is also a fixed name.
rndr <- function(x, name, ...) {
    if (length(x) == 0) {
        y <- d[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
		if(shapiro.test(y)$p.value<0.05){
			if(length(levels(droplevels(d$Group_pvalue)))>2){
            		p <- kruskal.test(y ~ droplevels(d$Group_pvalue))$p.value
			}else{
				p <- wilcox.test(y ~ droplevels(d$Group_pvalue))$p.value
			}
			p=paste(scales::pvalue(p),"**",sep="")
		}else{
			if(length(levels(droplevels(d$Group_pvalue)))>2){
				p = summary(aov(y ~ droplevels(d$Group_pvalue)))[[1]][1,5]
			}else{
				p = t.test(y ~ droplevels(d$Group_pvalue))$p.value
			}
			p=scales::pvalue(p)
		}
        } else {
		yd=table(y, droplevels(d$Group_pvalue))
		if(sum(yd<5)/length(yd)>0.2 | sum(yd==0)>0){
			p <- fisher.test(table(y, droplevels(d$Group_pvalue)))$p.value
			p=paste(scales::pvalue(p),"*",sep="")
		}else{
            	p <- chisq.test(table(y, droplevels(d$Group_pvalue)))$p.value
			p=scales::pvalue(p)
		}
        }
        s[2] <- sub("<", "&lt;", p)
        s
    } else {
        render.default(x=x, name=name, ...)
    }
}

##default table1 function
rndr.strat <- function(label, n, ...) {
    ifelse(n==0, label, render.strat.default(label, n, ...))
}

##define what to show for continuous varialbe. here is Mean SD and Median Q1 Q3
my.render.cont <- function(x) {
    with(stats.apply.rounding(stats.default(x), digits=2), c("",
        "Mean (SD)"=sprintf("%s (&plusmn;%s)", MEAN, SD), "Median [Q1, Q3]"=sprintf("%s [%s, %s]", MEDIAN, q25,q75)))
}

##Define the group/column variable in Table 1. 
d$Group=factor(d$Treatment, levels=c("Tr1", "Tr2", "Tr3"), labels=c("Tr1", "Tr2", "Tr3"))

##Define the group/column variable in Table 1, when a p-value column is needed.
d$Group_pvalue=factor(d$Treatment, levels=c("Tr1", "Tr2", "Tr3", "P-value"), labels=c("Tr1", "Tr2", "Tr3", "P-value"))

##This is the Table 1 working function. No p-value/comparion among column variable.
table1(~Age + Sex + DM + HTN + CAD + Afib + BMI +
		Dyslipidemia + Current.Smoking, data=d, render.continuous=my.render.cont)

##This is the Table 1 working function, with p-value/comparion among column variable.
table1(~Age + Sex + DM + HTN + CAD + Afib + BMI +
		Dyslipidemia + Current.Smoking | Group_pvalue, data=d, droplevels=F,
		render.continuous=my.render.cont, render.strat=rndr.strat, overall=T, render=rndr)

