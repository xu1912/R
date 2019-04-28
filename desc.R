sum_conti=function(vb){
	dt=data.frame(d[,colnames(d)==vb], d$tag)
	print(t.test(dt[,1]~dt[,2]))
	print(sd(dt[dt[,2]==0,1],na.rm=T))
	print(sd(dt[dt[,2]==1,1],na.rm=T))
}

sum_conti("Plt")
sum_conti("Total_chol")
sum_conti("LDL")
sum_conti("HDL")
sum_conti("TG")
sum_conti("AG")
sum_conti("HbA1c")
sum_conti("ICH")
