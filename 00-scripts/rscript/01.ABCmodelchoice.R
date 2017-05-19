source('/home/quentin/script/rscript/cv4abc.R')

target=as.numeric(read.table("OBS.ABC.stat.txt",skip=2,h=F))

nlinesFul=1e6

M_IM1=matrix(scan("im.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul) #les 1e6 simuls sont mergés ensemble au préalable avec du bash et sed.
M_IM2=matrix(scan("im.heterom.homon.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul) 
M_IM3=matrix(scan("im.homom.heteron.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul) 
M_IM4=matrix(scan("im.heterom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 
M_AM1=matrix(scan("am.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul) #les 1e6 sAMuls sont mergés ensemble au préalable avec du bash et sed.
M_AM2=matrix(scan("am.heterom.homon.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul) 
M_AM3=matrix(scan("am.homom.heteron.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul) 
M_AM4=matrix(scan("am.heterom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 
M_SC1=matrix(scan("sc.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul) #les 1e6 sSCuls sont mergés ensemble au préalable avec du bash et sed.
M_SC2=matrix(scan("sc.heterom.homon.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul) 
M_SC3=matrix(scan("sc.homom.heteron.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul) 
M_SC4=matrix(scan("sc.heterom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 
M_SI1=matrix(scan("si.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul) #les 1e6 sSIuls sont mergés ensemble au préalable avec du bash et sed.
M_SI3=matrix(scan("si.homom.heteron.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul) 

colSums(is.na(M_IM1))
means2 <- colMeans(M_IM1, na.rm=TRUE)
for (j in 1:ncol(M_IM1)){
     M_IM1[is.na(M_IM1[, j]), j] <- means2[j]
 }
colSums(is.na(M_IM2))
means2 <- colMeans(M_IM2, na.rm=TRUE)
for (j in 1:ncol(M_IM2)){
     M_IM2[is.na(M_IM2[, j]), j] <- means2[j]
 }
colSums(is.na(M_IM3))
means1 <- colMeans(M_IM3, na.rm=TRUE)
for (j in 1:ncol(M_IM3)){
     M_IM3[is.na(M_IM3[, j]), j] <- means1[j]
 }
colSums(is.na(M_IM4))
means3 <- colMeans(M_IM4, na.rm=TRUE)
for (j in 1:ncol(M_IM4)){
     M_IM4[is.na(M_IM4[, j]), j] <- means3[j]
 }

means2 <- colMeans(M_AM1, na.rm=TRUE)
for (j in 1:ncol(M_AM1)){
     M_AM1[is.na(M_AM1[, j]), j] <- means2[j]
 }

means2 <- colMeans(M_AM2, na.rm=TRUE)
for (j in 1:ncol(M_AM2)){
     M_AM2[is.na(M_AM2[, j]), j] <- means2[j]
 }

means1 <- colMeans(M_AM3, na.rm=TRUE)
for (j in 1:ncol(M_AM3)){
     M_AM3[is.na(M_AM3[, j]), j] <- means1[j]
 }

means3 <- colMeans(M_AM4, na.rm=TRUE)
for (j in 1:ncol(M_AM4)){
     M_AM4[is.na(M_AM4[, j]), j] <- means3[j]
 }

means2 <- colMeans(M_SC1, na.rm=TRUE)
for (j in 1:ncol(M_SC1)){
     M_SC1[is.na(M_SC1[, j]), j] <- means2[j]
 }

means2 <- colMeans(M_SC2, na.rm=TRUE)
for (j in 1:ncol(M_SC2)){
     M_SC2[is.na(M_SC2[, j]), j] <- means2[j]
 }
means1 <- colMeans(M_SC3, na.rm=TRUE)
for (j in 1:ncol(M_SC3)){
     M_SC3[is.na(M_SC3[, j]), j] <- means1[j]
 }
means3 <- colMeans(M_SC4, na.rm=TRUE)
for (j in 1:ncol(M_SC4)){
     M_SC4[is.na(M_SC4[, j]), j] <- means3[j]
 }

means2 <- colMeans(M_SI1, na.rm=TRUE)
for (j in 1:ncol(M_SI1)){
     M_SI1[is.na(M_SI1[, j]), j] <- means2[j]
 }
means1 <- colMeans(M_SI3, na.rm=TRUE)
for (j in 1:ncol(M_SI3)){
     M_SI3[is.na(M_SI3[, j]), j] <- means1[j]
 }

M_IM1b=M_IM1[,-c(1:3,12,13,18:25)]
M_IM2b=M_IM2[,-c(1:3,12,13,18:25)]
M_IM3b=M_IM3[,-c(1:3,12,13,18:25)]
M_IM4b=M_IM4[,-c(1:3,12,13,18:25)]
M_SC1b=M_SC1[,-c(1:3,12,13,18:25)]
M_SC2b=M_SC2[,-c(1:3,12,13,18:25)]
M_SC3b=M_SC3[,-c(1:3,12,13,18:25)]
M_SC4b=M_SC4[,-c(1:3,12,13,18:25)]
M_AM1b=M_AM1[,-c(1:3,12,13,18:25)]
M_AM2b=M_AM2[,-c(1:3,12,13,18:25)]
M_AM3b=M_AM3[,-c(1:3,12,13,18:25)]
M_AM4b=M_AM4[,-c(1:3,12,13,18:25)]
M_SI1b=M_SI1[,-c(1:3,12,13,18:25)]
M_SI3b=M_SI3[,-c(1:3,12,13,18:25)]

target=target[-c(1:3,12,13,18:25)]
obs=matrix(rep(target[1:19],20), byrow=T, nrow=20)

#Comparaison All:
x=as.factor(c(rep("A",nlinesFul),rep("B",nlinesFul),rep("C",nlinesFul),rep("D",nlinesFul),rep("E",nlinesFul),rep("F",nlinesFul),rep("G",nlinesFul),
			rep("H",nlinesFul),rep("I",nlinesFul),rep("J",nlinesFul),rep("K",nlinesFul),rep("L",nlinesFul),rep("M",nlinesFul),rep("N",nlinesFul)))

z=rbind(M_SI1b,M_SI3b,M_IM1b,M_IM2b,M_IM3b,M_IM4b,M_AM1b,M_AM2b,M_AM3b,M_AM4b,M_SC1b,M_SC2b,M_SC3b,M_SC4b)
res5=model_selection_abc_nnet(target=obs, x=x, sumstat=z, tol=3500/(length(x)), noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output="ALL.MODELS.29.06.16")

mean_all=apply(res5,2,mean)
sd_all=apply(res5,2,sd)
write.table(mean_all,"mean_all",quote=F,row.names=c("SI1","SI2","AM1","AM2","AM3","AM4","IM1","IM2","IM3","IM4","SC1","SC2","SC3","SC4"),col.names=F)
