library("pheatmap")

d2=read.xlsx("ap_ans.xlsx", header=T, sheetIndex=2)
d2m=as.matrix(d2[,-1])
d2m[upper.tri(d2m, diag = FALSE)]=t(d2m)[upper.tri(t(d2m), diag = FALSE)]

pheatmap((d2m), labels_col=d2$Name, labels_row=d2$Name, cluster_rows=F, 
		cluster_cols=F, color=colorRampPalette(c("white","blue"))(100),
		fontsize_row = 15,fontsize_col = 15)

pheatmap((d2m), labels_col=d2$Name, labels_row=d2$Name, cluster_rows=T, 
		cluster_cols=T, color=colorRampPalette(c("white","blue"))(100),
		fontsize_row = 15,fontsize_col = 15)
