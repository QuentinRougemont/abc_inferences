##!/usr/bin/Rscript
if("abc" %in% rownames(installed.packages()) == FALSE) {install.packages("abc") }
library(abc)

target=as.numeric(read.table("OBS.ABC.stat.txt",skip=2,h=F))
target=target[-c(1:3,12,13,18:25)]
obs=target[1:19]

nlinesFul <- as.numeric(strsplit(system("wc -l sc.homom.homon.ABC.stat.txt ", intern=T), " ")[[1]][1])
M_SC3 <- matrix(scan("sc.homom.homon.ABC.stat.txt"), byrow=T, nrow=nlinesFul)

means1 <- colMeans(M_SC3, na.rm=TRUE)
for (j in 1:ncol(M_SC3)){
     M_SC3[is.na(M_SC3[, j]), j] <- means1[j]
}

M_SC3b <- M_SC3[,-c(1:3,12,13,18:25)]

priorfile <- matrix(scan("sc.homom.homon.priorfile.txt"), byrow=T, nrow=nlinesFul)
colnames(priorfile) = c("N1","N2","Na","Tsplit","Tsc","M1","M2")

abc.postSC <- abc(target  = obs, param = priorfile, sumstat = M_SC3b, tol = 1000/1e6, transf=c(rep("logit",7) ),
         logit.bounds = rbind(range(priorfile[, 1]),range(priorfile[, 2]), range(priorfile[, 3]), range(priorfile[, 4]), range(priorfile[, 5]),
         range(priorfile[, 6]),range(priorfile[, 7])),
         hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

simpost<-500000
newsamp<-sample(1:(dim(abc.postSC$adj)[1]),size=simpost,replace=T,prob=abc.postSC$weights) 
newsamp<-abc.postSC$adj[newsamp,]
write.table(newsamp,"prior.from.post.homo", quote=F, row.names=F, col.names=F)
