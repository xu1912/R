library(xlsx)
library(edgeR)
library(network)
d=read.xlsx("TMT_data.xlsx", sheetIndex=2, header=T)
d_blank=d[which(d$Accession=="Blank"),]
d_n=d[-which(d$Accession=="Blank"),]
rownames(d_n)=d_n$Accession

data_raw=t(d_n[,4:1614])

colnames(data_raw)=d_n$Accession

boxplot(log2(data_raw), col = as.color(d_n$pool), 
        notch = TRUE, main = 'Raw data',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity')
legend(-2, 16.6, c("Pool1","Pool2","Pool3","Pool4","Pool5","Pool6"), pch = 16, col = c(1,2,3,4,5,6), horiz=T)

plotMDS(log2(data_raw),  col = as.color(d_n$pool), main = "PC clusters group by pools")
legend(5.8, -0.35, c("Pool1","Pool2","Pool3","Pool4","Pool5","Pool6"), pch = 16, col = c(1,2,3,4,5,6), bty = "n")

d_blank[is.na(d_blank)] <- 0

for(i in 1:nrow(d_n)){
	
	tv=d_n[i,4:1614]-d_blank[which(d_blank$pool==d_n[i,3]), 4:1614]
	tv[tv<0]=0
	d_n[i,4:1614]=tv

}


data_bk=t(d_n[,4:1614])

colnames(data_bk)=d_n$Accession

format(round(colSums(data_bk, na.rm=T), digits = 0), big.mark = ",")


data_mix=data_bk[,which(d_n$race=="Mix")]
data_mix_c=data_mix[complete.cases(data_mix),]

data_neg=data_bk[,which(d_n$race=="Negative")]
data_neg_c=data_neg[complete.cases(data_neg),]


nm_mix=colSums(data_mix_c)
nm_neg=colSums(data_neg_c)

nm_neg=nm_neg*mean(nm_mix)/mean(nm_neg)
nm_v=(nm_neg+nm_mix)/2
norm_facts=mean(nm_v)/nm_v
names(norm_facts)=c("pool1", "pool2", "pool3", "pool4", "pool5", "pool6")

for(i in 1:nrow(d_n)){
	d_n[i,4:1614]=d_n[i,4:1614]*norm_facts[which(names(norm_facts)==d_n[i,3])]
}


data_nm=t(d_n[,4:1614])

colnames(data_nm)=d_n$Accession

boxplot(log2(data_nm), col = as.color(d_n$pool), 
        notch = TRUE, main = 'Normalized data',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity')
legend(-2, 17, c("Pool1","Pool2","Pool3","Pool4","Pool5","Pool6"), pch = 16, col = c(1,2,3,4,5,6), horiz=T)

plotDensities(log2(data_nm), col = rep(c('red', 'green', 'blue'), 6), main = 'Raw data')

data_nm_c=data_nm[complete.cases(data_nm),]
sl_tmm <- calcNormFactors(data_nm_c)
data_nm_tmm <- sweep(data_nm, 2, sl_tmm, FUN = "/")

boxplot(log2(data_nm_tmm), col = as.color(d_n$pool), 
        notch = TRUE, main = 'Normalized data',
        xlab = 'TMT Samples', ylab = 'log2 of Intensity')
legend(-2, 17, c("Pool1","Pool2","Pool3","Pool4","Pool5","Pool6"), pch = 16, col = c(1,2,3,4,5,6), horiz=T)



plotDensities(log2(data_nm_tmm), col = rep(c('red', 'green', 'blue'), 6), main = 'Raw data')

plotMDS(log2(data_nm_tmm),  col = as.color(d_n$pool), main = "PC clusters group by pools")
legend(6.5, -0.6, c("Pool1","Pool2","Pool3","Pool4","Pool5","Pool6"), pch = 16, col = c(1,2,3,4,5,6), bty = "n")

data_nm_tmm=data_nm_tmm[ , !colnames(data_nm_tmm) %in% c("Mix_p1", "Mix_p2", "Mix_p3", "Mix_p4", "Mix_p5", "Mix_p6", "Negative_p1", "Negative_p2","Negative_p3","Negative_p4","Negative_p5","Negative_p6")]

str(data_nm_tmm) 
kw_pvalue=c()
kw_e_pvalue=c()
kw_ee_pvalue=c()
aov_pvalue=c()
aov_e_pvalue=c()
aov_ee_pvalue=c()
missR=c()
for(i in 1:nrow(data_nm_tmm)){
	
	protein=data_nm_tmm[i,]
	race=d_n$race[match(names(protein), d_n$Accession)]
	dt=data.frame(protein, race)
	missR[i]=sum(is.na(dt$protein))/46
	kw_pvalue[i]=kruskal.test(protein~race, data=dt)$p.value
	dt$log2Protein=log2(dt$protein)
	fit.l=aov(log2Protein~race, data=dt[is.finite(dt$log2Protein),])
	aov_pvalue[i]=summary(fit.l)[[1]]$"Pr(>F)"[1]
	dt_e=dt[names(protein)!= "AA9",]
	kw_e_pvalue[i]=kruskal.test(protein~race, data=dt_e)$p.value
	fit.le=aov(log2Protein~race, data=dt_e[is.finite(dt_e$log2Protein),])
	aov_e_pvalue[i]=summary(fit.le)[[1]]$"Pr(>F)"[1]		
	dt_ee=dt[names(protein)!= "AA9" & names(protein)!= "A1",]
	kw_ee_pvalue[i]=kruskal.test(protein~race, data=dt_ee)$p.value
	fit.lee=aov(log2Protein~race, data=dt_ee[is.finite(dt_ee$log2Protein),])
	aov_ee_pvalue[i]=summary(fit.lee)[[1]]$"Pr(>F)"[1]		

}

p.adjust(aov_pvalue[missR<0.1], method="BH")

res=data.frame(kw_pvalue, kw_e_pvalue, kw_ee_pvalue, aov_pvalue, aov_e_pvalue, aov_ee_pvalue, missR)
row.names(res)=row.names(data_nm_tmm)

write.xlsx(res, "protein_difference_in_race_groups.xlsx", row.names=T, sheetName="p-values2", col.names=T, append=T)
write.xlsx(data_nm_tmm, "protein_after_normalization.xlsx", row.names=T, sheetName="protein level", col.names=T)


res_sig=res[which(res$aov_pvalue<0.05),]

row.names(res_sig)
library(DescTools)

missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm[which(row.names(data_nm_tmm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	dt=data.frame(protein, race)
	missR[i]=sum(is.na(dt$protein))/46
	dt$log2Protein=log2(dt$protein)
	dti=dt[is.finite(dt$log2Protein),]
	fit.l=aov(log2Protein~race, data=dti)
	aov_pvalue[i]=summary(fit.l)[[1]]$"Pr(>F)"[1]
	tukey.test <- TukeyHSD(fit.l, ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, "post_hoc_all_samples_Tukey_pairwise.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])
	dti$race=as.factor(dti$race)
	dunnett.res=DunnettTest(log2Protein~race, data=dti, control="White")$White
	colnames(dunnett.res)[2]="95% CI lower bound"
	colnames(dunnett.res)[3]="95% CI upper bound"
	dunnett.res=data.frame(dunnett.res, check.names=F)
	write.xlsx(dunnett.res, "post_hoc_all_samples_Dunnett_White.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])

}

write.xlsx(res, "protein_difference_in_race_groups.xlsx", row.names=T, sheetName="p-values2", col.names=T, append=T)
write.xlsx(data_nm_tmm, "protein_after_normalization.xlsx", row.names=T, sheetName="protein level", col.names=T)

res_sig=res[which(res$aov_e_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm[which(row.names(data_nm_tmm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	dt=data.frame(protein, race)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9",]
	dti=dt_e[is.finite(dt_e$log2Protein),]
	fit.l=aov(log2Protein~race, data=dti)
	aov_pvalue[i]=summary(fit.l)[[1]]$"Pr(>F)"[1]
	tukey.test <- TukeyHSD(fit.l, ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, "post_hoc_exclude_AA9_Tukey_pairwise.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])
	dti$race=as.factor(dti$race)
	dunnett.res=DunnettTest(log2Protein~race, data=dti, control="White")$White
	colnames(dunnett.res)[2]="95% CI lower bound"
	colnames(dunnett.res)[3]="95% CI upper bound"
	dunnett.res=data.frame(dunnett.res, check.names=F)
	write.xlsx(dunnett.res, "post_hoc_exclude_AA9_Dunnett_White.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])

}

res_sig=res[which(res$aov_ee_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm[which(row.names(data_nm_tmm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	dt=data.frame(protein, race)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9" & names(protein)!= "A1",]
	dti=dt_e[is.finite(dt_e$log2Protein),]
	fit.l=aov(log2Protein~race, data=dti)
	aov_pvalue[i]=summary(fit.l)[[1]]$"Pr(>F)"[1]
	tukey.test <- TukeyHSD(fit.l, ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, "post_hoc_exclude_AA9_AA1_Tukey_pairwise.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])
	dti$race=as.factor(dti$race)
	dunnett.res=DunnettTest(log2Protein~race, data=dti, control="White")$White
	colnames(dunnett.res)[2]="95% CI lower bound"
	colnames(dunnett.res)[3]="95% CI upper bound"
	dunnett.res=data.frame(dunnett.res, check.names=F)
	write.xlsx(dunnett.res, "post_hoc_exclude_AA9_AA1_Dunnett_White.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])

}



p.adjust(aov_e_pvalue, method="BH")

data_nm_tmm=data_nm_tmm[ , !colnames(data_nm_tmm) %in% c("Mix_p1", "Mix_p2", "Mix_p3", "Mix_p4", "Mix_p5", "Mix_p6", "Negative_p1", "Negative_p2","Negative_p3","Negative_p4","Negative_p5","Negative_p6")]
str(data_nm_tmm)

data_nm_tmm_com=data_nm_tmm[complete.cases(data_nm_tmm),]

res.pca <- prcomp(t(data_nm_tmm_com), scale = TRUE)
fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

fviz_pca_ind(res.pca, repel = TRUE)

pcmx=res.pca$x[,1:2]
 
aov_pvalue=c()
aov_e_pvalue=c()
aov_ee_pvalue=c()
missR=c()
for(i in 1:nrow(data_nm_tmm)){
	
	protein=data_nm_tmm[i,]
	race=d_n$race[match(names(protein), d_n$Accession)]
	pc12=pcmx[match(names(protein),rownames(pcmx)), ]
	dt=data.frame(protein, race, pc12)
	missR[i]=sum(is.na(dt$protein))/46
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9",]
	dt_ee=dt[names(protein)!= "AA9" & names(protein)!= "A1",]
	fit.l=aov(log2Protein~race+PC1, data=dt[is.finite(dt$log2Protein),])
	if(nrow(summary(fit.l)[[1]])<3){
		aov_pvalue[i]=NA
		aov_e_pvalue[i]=NA
	}else{
		aov_pvalue[i]=summary(fit.l)[[1]]$"Pr(>F)"[1]
		fit.le=aov(log2Protein~race+PC1, data=dt_e[is.finite(dt_e$log2Protein),])
		aov_e_pvalue[i]=summary(fit.le)[[1]]$"Pr(>F)"[1]
		fit.lee=aov(log2Protein~race+PC1, data=dt_ee[is.finite(dt_ee$log2Protein),])
		aov_ee_pvalue[i]=summary(fit.lee)[[1]]$"Pr(>F)"[1]
	}	
}

p.adjust(aov_pvalue[missR<0.1], method="BH")

res=data.frame(aov_pvalue, aov_e_pvalue, aov_ee_pvalue, missR)
row.names(res)=row.names(data_nm_tmm)
write.xlsx(res, "protein_difference_in_race_groups.xlsx", row.names=T, sheetName="pc adjusted p-values", col.names=T, append=T)

which(res$aov_pvalue<0.05)
which(res$aov_e_pvalue<0.05)
which(res$aov_ee_pvalue<0.05)

res_sig=res[which(res$aov_ee_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm[which(row.names(data_nm_tmm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	pc12=pcmx[match(names(protein),rownames(pcmx)), ]
	dt=data.frame(protein, race, pc12)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9" & names(protein)!= "A1",]
	dti=dt_e[is.finite(dt_e$log2Protein),]
	dti$race=as.factor(dti$race)
	if (sum("White" %in% dti$race)>0){
		relevel(dti$race, ref="White")
		dunnett_tag=1
	}else{
		dunnett_tag=0
	}
	fit.l=aov(log2Protein~race+PC1, data=dti)
	tukey.test <- TukeyHSD(fit.l, "race", ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, paste("post_hoc_PC_adjusted_exclude_AA9_AA1/", row.names(res_sig)[i],".xlsx"), row.names=T, col.names=T, append=T, sheetName="Tukey_pairwise")

	if(dunnett_tag==1){
		dunnett.res=contrast(lsmeans(fit.l, list(~race)), method="trt.vs.ctrl", adjust="dunnett", ref="White")[[1]]
		dunnett.res=data.frame(dunnett.res, check.names=F)
		write.xlsx(dunnett.res[,-5], paste("post_hoc_PC_adjusted_exclude_AA9_AA1/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}else{
		dunnett.res="Missing data in White population"
		write.xlsx(dunnett.res, paste("post_hoc_PC_adjusted_exclude_AA9_AA1/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}

}


res_sig=res[which(res$aov_e_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm[which(row.names(data_nm_tmm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	pc12=pcmx[match(names(protein),rownames(pcmx)), ]
	dt=data.frame(protein, race, pc12)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9",]
	dti=dt_e[is.finite(dt_e$log2Protein),]
	dti$race=as.factor(dti$race)
	if (sum("White" %in% dti$race)>0){
		relevel(dti$race, ref="White")
		dunnett_tag=1
	}else{
		dunnett_tag=0
	}
	fit.l=aov(log2Protein~race+PC1, data=dti)
	tukey.test <- TukeyHSD(fit.l, "race", ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, paste("post_hoc_PC_adjusted_exclude_AA9/", row.names(res_sig)[i],".xlsx"), row.names=T, col.names=T, append=T, sheetName="Tukey_pairwise")

	if(dunnett_tag==1){
		dunnett.res=contrast(lsmeans(fit.l, list(~race)), method="trt.vs.ctrl", adjust="dunnett", ref="White")[[1]]
		dunnett.res=data.frame(dunnett.res, check.names=F)
		write.xlsx(dunnett.res[,-5], paste("post_hoc_PC_adjusted_exclude_AA9/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}else{
		dunnett.res="Missing data in White population"
		write.xlsx(dunnett.res, paste("post_hoc_PC_adjusted_exclude_AA9/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}

}

res_sig=res[which(res$aov_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm[which(row.names(data_nm_tmm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	pc12=pcmx[match(names(protein),rownames(pcmx)), ]
	dt=data.frame(protein, race, pc12)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt
	dti=dt_e[is.finite(dt_e$log2Protein),]
	dti$race=as.factor(dti$race)
	if (sum("White" %in% dti$race)>0){
		relevel(dti$race, ref="White")
		dunnett_tag=1
	}else{
		dunnett_tag=0
	}
	fit.l=aov(log2Protein~race+PC1, data=dti)
	tukey.test <- TukeyHSD(fit.l, "race", ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, paste("post_hoc_PC_adjusted_all_samples/", row.names(res_sig)[i],".xlsx"), row.names=T, col.names=T, append=T, sheetName="Tukey_pairwise")

	if(dunnett_tag==1){
		dunnett.res=contrast(lsmeans(fit.l, list(~race)), method="trt.vs.ctrl", adjust="dunnett", ref="White")[[1]]
		dunnett.res=data.frame(dunnett.res, check.names=F)
		write.xlsx(dunnett.res[,-5], paste("post_hoc_PC_adjusted_all_samples/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}else{
		dunnett.res="Missing data in White population"
		write.xlsx(dunnett.res, paste("post_hoc_PC_adjusted_all_samples/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}

}



##Remove AA12, W12, and I12
data_nm_tmm_rm=data_nm_tmm[ , !colnames(data_nm_tmm) %in% c("AA12", "W12", "I12", "Mix_p1", "Mix_p2", "Mix_p3", "Mix_p4", "Mix_p5", "Mix_p6", "Negative_p1", "Negative_p2","Negative_p3","Negative_p4","Negative_p5","Negative_p6")]
str(data_nm_tmm_rm) 
kw_pvalue=c()
kw_e_pvalue=c()
kw_ee_pvalue=c()
aov_pvalue=c()
aov_e_pvalue=c()
aov_ee_pvalue=c()
missR=c()
for(i in 1:nrow(data_nm_tmm_rm)){
	
	protein=data_nm_tmm_rm[i,]
	race=d_n$race[match(names(protein), d_n$Accession)]
	dt=data.frame(protein, race)
	missR[i]=sum(is.na(dt$protein))/43
	kw_pvalue[i]=kruskal.test(protein~race, data=dt)$p.value
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9",]
	dt_ee=dt[names(protein)!= "AA9" & names(protein)!= "A1",]
	fit.l=aov(log2Protein~race, data=dt[is.finite(dt$log2Protein),])
	if(nrow(summary(fit.l)[[1]])==1){
		aov_pvalue[i]=NA
		aov_e_pvalue[i]=NA
	}else{
		aov_pvalue[i]=summary(fit.l)[[1]]$"Pr(>F)"[1]
		fit.le=aov(log2Protein~race, data=dt_e[is.finite(dt_e$log2Protein),])
		aov_e_pvalue[i]=summary(fit.le)[[1]]$"Pr(>F)"[1]
		fit.lee=aov(log2Protein~race, data=dt_ee[is.finite(dt_ee$log2Protein),])
		aov_ee_pvalue[i]=summary(fit.lee)[[1]]$"Pr(>F)"[1]
	}
	
	kw_e_pvalue[i]=kruskal.test(protein~race, data=dt_e)$p.value
	kw_ee_pvalue[i]=kruskal.test(protein~race, data=dt_ee)$p.value
		
}

p.adjust(aov_pvalue[missR<0.1], method="BH")

p.adjust(aov_e_pvalue, method="BH")
res=data.frame(kw_pvalue, kw_e_pvalue, kw_ee_pvalue, aov_pvalue, aov_e_pvalue, aov_ee_pvalue, missR)
row.names(res)=row.names(data_nm_tmm_rm)
write.xlsx(res, "protein_difference_in_race_groups_rm_serous.xlsx", row.names=T, sheetName="p-values", col.names=T, append=T)


res_sig=res[which(res$aov_ee_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm_rm[which(row.names(data_nm_tmm_rm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	dt=data.frame(protein, race)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9" & names(protein)!= "A1",]
	dti=dt_e[is.finite(dt_e$log2Protein),]
	fit.l=aov(log2Protein~race, data=dti)
	tukey.test <- TukeyHSD(fit.l, ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, "Remove_serous/post_hoc_exclude_AA9_AA1_Tukey_pairwise.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])
	dti$race=as.factor(dti$race)
	dunnett.res=DunnettTest(log2Protein~race, data=dti, control="White")$White
	colnames(dunnett.res)[2]="95% CI lower bound"
	colnames(dunnett.res)[3]="95% CI upper bound"
	dunnett.res=data.frame(dunnett.res, check.names=F)
	write.xlsx(dunnett.res, "Remove_serous/post_hoc_exclude_AA9_AA1_Dunnett_White.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])

}

res_sig=res[which(res$aov_e_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm_rm[which(row.names(data_nm_tmm_rm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	dt=data.frame(protein, race)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9",]
	dti=dt_e[is.finite(dt_e$log2Protein),]
	fit.l=aov(log2Protein~race, data=dti)
	tukey.test <- TukeyHSD(fit.l, ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, "Remove_serous/post_hoc_exclude_AA9_Tukey_pairwise.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])
	dti$race=as.factor(dti$race)
	dunnett.res=DunnettTest(log2Protein~race, data=dti, control="White")$White
	colnames(dunnett.res)[2]="95% CI lower bound"
	colnames(dunnett.res)[3]="95% CI upper bound"
	dunnett.res=data.frame(dunnett.res, check.names=F)
	write.xlsx(dunnett.res, "Remove_serous/post_hoc_exclude_AA9_Dunnett_White.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])

}

res_sig=res[which(res$aov_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm_rm[which(row.names(data_nm_tmm_rm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	dt=data.frame(protein, race)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt
	dti=dt_e[is.finite(dt_e$log2Protein),]
	fit.l=aov(log2Protein~race, data=dti)
	tukey.test <- TukeyHSD(fit.l, ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, "Remove_serous/post_hoc_all_samples_Tukey_pairwise.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])
	dti$race=as.factor(dti$race)
	dunnett.res=DunnettTest(log2Protein~race, data=dti, control="White")$White
	colnames(dunnett.res)[2]="95% CI lower bound"
	colnames(dunnett.res)[3]="95% CI upper bound"
	dunnett.res=data.frame(dunnett.res, check.names=F)
	write.xlsx(dunnett.res, "Remove_serous/post_hoc_all_samples_Dunnett_White.xlsx", row.names=T, col.names=T, append=T, sheetName=row.names(res_sig)[i])

}

library(factoextra)
data_nm_tmm_rm_com=data_nm_tmm_rm[complete.cases(data_nm_tmm_rm),]
res.pca <- prcomp(t(data_nm_tmm_rm_com), scale = TRUE)
fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

fviz_pca_ind(res.pca, repel = TRUE)

pcmx=res.pca$x[,1:2]

str(data_nm_tmm_rm) 
aov_pvalue=c()
aov_e_pvalue=c()
aov_ee_pvalue=c()
missR=c()
for(i in 1:nrow(data_nm_tmm_rm)){
	
	protein=data_nm_tmm_rm[i,]
	race=d_n$race[match(names(protein), d_n$Accession)]
	pc12=pcmx[match(names(protein),rownames(pcmx)), ]
	dt=data.frame(protein, race, pc12)
	missR[i]=sum(is.na(dt$protein))/43
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9",]
	dt_ee=dt[names(protein)!= "AA9" & names(protein)!= "A1",]
	fit.l=aov(log2Protein~race+PC1, data=dt[is.finite(dt$log2Protein),])
	if(nrow(summary(fit.l)[[1]])<3){
		aov_pvalue[i]=NA
		aov_e_pvalue[i]=NA
	}else{
		aov_pvalue[i]=summary(fit.l)[[1]]$"Pr(>F)"[1]
		fit.le=aov(log2Protein~race+PC1, data=dt_e[is.finite(dt_e$log2Protein),])
		aov_e_pvalue[i]=summary(fit.le)[[1]]$"Pr(>F)"[1]
		fit.lee=aov(log2Protein~race+PC1, data=dt_ee[is.finite(dt_ee$log2Protein),])
		aov_ee_pvalue[i]=summary(fit.lee)[[1]]$"Pr(>F)"[1]
	}	
}

p.adjust(aov_pvalue[missR<0.1], method="BH")

res=data.frame(aov_pvalue, aov_e_pvalue, aov_ee_pvalue, missR)
row.names(res)=row.names(data_nm_tmm_rm)
write.xlsx(res, "protein_difference_in_race_groups_rm_serous.xlsx", row.names=T, sheetName="pc adjusted p-values", col.names=T, append=T)

res_sig=res[which(res$aov_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm_rm[which(row.names(data_nm_tmm_rm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	pc12=pcmx[match(names(protein),rownames(pcmx)), ]
	dt=data.frame(protein, race, pc12)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt
	dti=dt_e[is.finite(dt_e$log2Protein),]
	dti$race=as.factor(dti$race)
	if (sum("White" %in% dti$race)>0){
		relevel(dti$race, ref="White")
		dunnett_tag=1
	}else{
		dunnett_tag=0
	}
	fit.l=aov(log2Protein~race+PC1, data=dti)
	tukey.test <- TukeyHSD(fit.l, "race", ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, paste("Remove_serous/post_hoc_PC_adjusted_all_samples/", row.names(res_sig)[i],".xlsx"), row.names=T, col.names=T, append=T, sheetName="Tukey_pairwise")

	if(dunnett_tag==1){
		dunnett.res=contrast(lsmeans(fit.l, list(~race)), method="trt.vs.ctrl", adjust="dunnett", ref="White")[[1]]
		dunnett.res=data.frame(dunnett.res, check.names=F)
		write.xlsx(dunnett.res[,-5], paste("Remove_serous/post_hoc_PC_adjusted_all_samples/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}else{
		dunnett.res="Missing data in White population"
		write.xlsx(dunnett.res, paste("Remove_serous/post_hoc_PC_adjusted_all_samples/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}

}

res_sig=res[which(res$aov_ee_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm_rm[which(row.names(data_nm_tmm_rm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	pc12=pcmx[match(names(protein),rownames(pcmx)), ]
	dt=data.frame(protein, race, pc12)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9" & names(protein)!= "A1",]
	dti=dt_e[is.finite(dt_e$log2Protein),]
	dti$race=as.factor(dti$race)
	if (sum("White" %in% dti$race)>0){
		relevel(dti$race, ref="White")
		dunnett_tag=1
	}else{
		dunnett_tag=0
	}
	fit.l=aov(log2Protein~race+PC1, data=dti)
	tukey.test <- TukeyHSD(fit.l, "race", ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, paste("Remove_serous/post_hoc_PC_adjusted_exclude_AA9_AA1/", row.names(res_sig)[i],".xlsx"), row.names=T, col.names=T, append=T, sheetName="Tukey_pairwise")

	if(dunnett_tag==1){
		dunnett.res=contrast(lsmeans(fit.l, list(~race)), method="trt.vs.ctrl", adjust="dunnett", ref="White")[[1]]
		dunnett.res=data.frame(dunnett.res, check.names=F)
		write.xlsx(dunnett.res[,-5], paste("Remove_serous/post_hoc_PC_adjusted_exclude_AA9_AA1/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}else{
		dunnett.res="Missing data in White population"
		write.xlsx(dunnett.res, paste("Remove_serous/post_hoc_PC_adjusted_exclude_AA9_AA1/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}

}


res_sig=res[which(res$aov_e_pvalue<0.05),]

row.names(res_sig)
missR=c()
for(i in 1:nrow(res_sig)){

	protein=data_nm_tmm_rm[which(row.names(data_nm_tmm_rm)==row.names(res_sig)[i]),]
	race=d_n$race[match(names(protein), d_n$Accession)]
	pc12=pcmx[match(names(protein),rownames(pcmx)), ]
	dt=data.frame(protein, race, pc12)
	dt$log2Protein=log2(dt$protein)
	dt_e=dt[names(protein)!= "AA9",]
	dti=dt_e[is.finite(dt_e$log2Protein),]
	dti$race=as.factor(dti$race)
	if (sum("White" %in% dti$race)>0){
		relevel(dti$race, ref="White")
		dunnett_tag=1
	}else{
		dunnett_tag=0
	}
	fit.l=aov(log2Protein~race+PC1, data=dti)
	tukey.test <- TukeyHSD(fit.l, "race", ordered = T)
	tukey.res=tukey.test$race
	colnames(tukey.res)[2]="95% CI lower bound"
	colnames(tukey.res)[3]="95% CI upper bound"
	tukey.res=data.frame(tukey.res, check.names=F)
	write.xlsx(tukey.res, paste("Remove_serous/post_hoc_PC_adjusted_exclude_AA9/", row.names(res_sig)[i],".xlsx"), row.names=T, col.names=T, append=T, sheetName="Tukey_pairwise")

	if(dunnett_tag==1){
		dunnett.res=contrast(lsmeans(fit.l, list(~race)), method="trt.vs.ctrl", adjust="dunnett", ref="White")[[1]]
		dunnett.res=data.frame(dunnett.res, check.names=F)
		write.xlsx(dunnett.res[,-5], paste("Remove_serous/post_hoc_PC_adjusted_exclude_AA9/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}else{
		dunnett.res="Missing data in White population"
		write.xlsx(dunnett.res, paste("Remove_serous/post_hoc_PC_adjusted_exclude_AA9/", row.names(res_sig)[i],".xlsx"), row.names=F, col.names=T, append=T, sheetName="Dunnett_White")
	}

}

