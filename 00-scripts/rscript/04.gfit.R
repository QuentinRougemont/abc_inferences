#!/usr/bin/Rscript

obs    <- read.table("../results/OBS.ABC.stat.txt", header=T)
target <-obs[-c(1:3,12,13,18:25)]

ppc_sc <- read.table('gfit.outfile', header=F)
ppc_sc <- ppc_sc[-c(1:3,12,13,18:25)]
vect   <-target

stat_names <-colnames(vect)

#plot stat obs et stat sim:
#Change the RIVER NAME in pdf file eahc time a new river is plotted
pdf("PPC_SC.pdf",width=40,height=40)
par(mfrow=c(5,4))
for (i in 1:ncol(ppc_sc)){
	hist(ppc_sc[,i], xlab=stat_names[i], main=stat_names[i])
	abline(v=vect[,i], col='red')
	pval=mean(abs(ppc_sc[,i]) >= abs(vect[,i]))
	write.table(t(pval),file="ppc_sc_Goodness.txt", append=T,quote=F,row.names=F, col.names=F)
	}
dev.off()
