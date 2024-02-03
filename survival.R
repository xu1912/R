library(survival)
library(xlsx)
library(stringr)
library(sqldf)

d_surv=read.xlsx("C:/Users/cxu2/Documents/OvarianCancer/TCGA_OS_PFS_OV.xlsx", sheetIndex=1, header=T)
d_pho=read.table("F:/CancerCenter/PDC/ov_tcga_pub2011/data_mrna_mirna_merged_zscores.txt", header=T, sep="\t")
inx=which(d_pho$Hugo_Symbol=="NADK2")
d_stat1=t(d_pho[inx,-c(1:2)])
colnames(d_stat1)=d_pho[inx,1]
rownames(d_stat1)=str_replace_all(rownames(d_stat1),"\\.01","")

d_surv$bcr_patient_barcode=str_replace_all(d_surv$bcr_patient_barcode,"-","\\.")

rownames(d_stat1) %in% d_surv$bcr_patient_barcode
ds=cbind(d_surv[match(rownames(d_stat1), d_surv$bcr_patient_barcode),], d_stat1)

ds$race_group="Other"
ds$race_group[ds$race=="WHITE"]="White"
ds$race_group[ds$race=="BLACK OR AFRICAN AMERICAN"]="Black"
ds$race_group[ds$race=="[Not Available]"]=NA

ds$stage="S1/2"
ds$stage[ds$clinical_stage=="[Not Available]"]=NA
ds$stage[ds$clinical_stage=="Stage IIIA"|ds$clinical_stage=="Stage IIIB"|ds$clinical_stage=="Stage IIIC"]="S3"
ds$stage[ds$clinical_stage=="Stage IV"]="S4"

fit=coxph(Surv(OS.time, OS)~age_at_initial_pathologic_diagnosis, data=ds)
summary(fit)
temp <- cox.zph(fit) 
print(temp)
fit=survdiff(Surv(OS.time, OS)~age_at_initial_pathologic_diagnosis, data=ds)
1 - pchisq(fit$chisq, length(fit$n) - 1)

##Optimal cutoff
cutoff=sort(ds$NADK2)
pval=c()
for(i in 2:(length(cutoff)-1)){

ds$NADK2_g="high"
ds$NADK2_g[ds$NADK2<cutoff[i]]="low"
ds$NADK2_g=factor(ds$NADK2_g, levels=c("low","high"))
fit=survdiff(Surv(OS.time, OS)~NADK2_g, data=ds)
pval[i-1]=1 - pchisq(fit$chisq, length(fit$n) - 1)

}
fit <- survfit(Surv(OS.time, OS) ~ NADK2_g, data = ds)
ggsurvplot(fit, data = ds, legend.title = "NADK2", legend.labs = c("Low", "High"), font.legend=c(14),
		pval = TRUE, pval.method = TRUE, conf.int = FALSE, xlab="Time (Day)",
		pval.coord=c(0.1,0.1), pval.method.coord=c(0.1,0.2))
