library(table1)
library(ggplot2)

d=read.table("demographic_data.txt", header=T, sep="\t", na.strings = "")
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

label(d$PageHour)="Page Hour"
label(d$AIS.Training)="AIS Training"
label(d$AIS.Training_pvalue)="AIS Training"

my.render.cont <- function(x) {
    with(stats.apply.rounding(stats.default(x), digits=2), c("",
        "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD), "Median [Q1, Q3]"=sprintf("%s [%s , %s]", MEDIAN, q25,q75)))
}

table1(~Age + Sex + PageHour + SystemsCategory + Neuro.Experience + Days.since.admit + AIS.Training + DM + HTN + CAD + Afib + FocalvsNon.focalsymptoms +
		Dyslipidemia + Current.Smoking + Prior.stroke + NIHSS + mRS +Type.of.Stroke + Stroke, data=d, render.continuous=my.render.cont)
 

rndr <- function(x, name, ...) {
    if (length(x) == 0) {
        y <- d[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ d$AIS.Training_pvalue)$p.value
	    p = summary(aov(y ~ droplevels(d$TYPE_pvalue)))[[1]][1,5]
        } else {
            p <- chisq.test(table(y, droplevels(d$AIS.Training_pvalue)))$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
        s
    } else {
        render.default(x=x, name=name, ...)
    }
}
rndr_fisher <- function(x, name, ...) {
    if (length(x) == 0) {
        y <- d[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ d$AIS.Training_pvalue)$p.value
        } else {
            p <- fisher.test(table(y, droplevels(d$AIS.Training_pvalue)))$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
        s
    } else {
        render.default(x=x, name=name, ...)
    }
}

rndr_fisher_simul <- function(x, name, ...) {
    if (length(x) == 0) {
        y <- d[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ d$AIS.Training_pvalue)$p.value
        } else {
            p <- fisher.test(table(y, droplevels(d$AIS.Training_pvalue)), simulate.p.value=TRUE)$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
        s
    } else {
        render.default(x=x, name=name, ...)
    }
}

rndr.strat <- function(label, n, ...) {
    ifelse(n==0, label, render.strat.default(label, n, ...))
}


table1(~Time.to.CT.MR.from.stroke.alert..min. + Time.to.Neurology.eval | AIS.Training_pvalue, data=d, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F)
table1(~IV.tPA + Time.to.IV.Tpa.from.last.known.normal + Time.to.IV.tPA.from.stroke.alert + Final.Diagnosis | AIS.Training_pvalue, data=d, droplevels=F, render=rndr_fisher, render.strat=rndr.strat, overall=F)

table1(~PageHour  + DM + HTN + CAD + Afib +
		Dyslipidemia + Current.Smoking + Prior.stroke + NIHSS + mRS +Type.of.Stroke + Stroke | AIS.Training_pvalue, 
		data=d, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F)
