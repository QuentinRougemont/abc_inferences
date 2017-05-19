##!/usr/bin/Rscript
if("abc" %in% rownames(installed.packages()) == FALSE) {install.packages("abc") }
library(abc)


nlocus <- 1
target=as.numeric(read.table("OBS.ABC.stat.txt",skip=2,h=F))
target=target[-c(1:3,12,13,18:25)]
obs=target[1:19]

nlinesFul <- as.numeric(strsplit(system("wc -l sc.heterom.heteron.ABC.stat.txt ", intern=T), " ")[[1]][1])
M_SC3 <- matrix(scan("sc.heterom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul)

means1 <- colMeans(M_SC3, na.rm=TRUE)
for (j in 1:ncol(M_SC3)){
     M_SC3[is.na(M_SC3[, j]), j] <- means1[j]
}

M_SC3b <- M_SC3[,-c(1:3,12,13,18:25)]

priorfile <- matrix(scan("sc.heterom.heteron.priorfile.txt"), byrow=T, nrow=nlinesFul)
colnames(priorfile) = c("N1","N2","Na","Tsplit","Tsc","Shape1Ne","Shape2Ne","propNtrN1","proprNtrN2","M1","M2","shape1M1","shape2M1","shape1M2","shape2M2","propNtrlM1","propNtrlM2")

abc.postSC <- abc(target  = obs, param = priorfile, sumstat = M_SC3b, tol = 1000/1e6, transf=c(rep("logit",17) ),
         logit.bounds = rbind(range(priorfile[, 1]),range(priorfile[, 2]), range(priorfile[, 3]), range(priorfile[, 4]), range(priorfile[, 5]),
         range(priorfile[, 6]),range(priorfile[, 7]),range(priorfile[, 8]),range(priorfile[, 9]),range(priorfile[, 10]), range(priorfile[, 11]),range(priorfile[,12]),range(priorfile[,13]),range(priorfile[,14]),range(priorfile[,15]),range(priorfile[,16]),range(priorfile[,17])),
         hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)


simpost <- 500000
simpost2 <- simpost/2
newsamp <- sample(1:(dim(abc.postSC$adj)[1]),size=simpost,replace=T,prob=abc.postSC$weights) 

newsamp <- abc.postSC$adj[newsamp,]
write.table(newsamp,"prior.from.post.het.het_v1", quote=F, row.names=F, col.names=F)

#n1.2 <- matrix(ncol = 2 , nrow = nrow(newsamp))
n1 <- rbeta(nlocus,newsamp[,6], newsamp[,7]) * newsamp[,1]
n2 <- rbeta(nlocus,newsamp[,6], newsamp[,7]) * newsamp[,2]
#1.tmp  <- sample(n1, simpost2)
#n2.tmp  <- sample(n2, simpost2)
#n1.tmp1 <- sample(newsamp[,1], simpost2)
#n2.tmp1 <- sample (newsamp[,2],simpost2)

newsamp[,1] <- n1
newsamp[,2] <- n2

write.table(newsamp,"prior.from.post.het.het_v2", quote=F, row.names=F, col.names=F)
