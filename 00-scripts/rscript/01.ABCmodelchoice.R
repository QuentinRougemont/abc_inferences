source('00-scripts/rscript/cv4abc.R')

colon_count = paste(" awk -F ' ' '{print NF }' ", "00.data/im.homom.homon.ABC.stat.txt" , " |head -1 " )
ncol = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])

target=as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))

#load simulations
IM.1=matrix(scan("00.data/im.homom.homon.ABC.stat.txt"    ), byrow=T, ncol=ncol)
IM.2=matrix(scan("00.data/im.heterom.homon.ABC.stat.txt"  ), byrow=T, ncol=ncol)
IM.3=matrix(scan("00.data/im.homom.heteron.ABC.stat.txt"  ), byrow=T, ncol=ncol)
IM.4=matrix(scan("00.data/im.heterom.heteron.ABC.stat.txt"), byrow=T, ncol=ncol)
AM.1=matrix(scan("00.data/am.homom.homon.ABC.stat.txt"    ), byrow=T, ncol=ncol)
AM.2=matrix(scan("00.data/am.heterom.homon.ABC.stat.txt"  ), byrow=T, ncol=ncol)
AM.3=matrix(scan("00.data/am.homom.heteron.ABC.stat.txt"  ), byrow=T, ncol=ncol)
AM.4=matrix(scan("00.data/am.heterom.heteron.ABC.stat.txt"), byrow=T, ncol=ncol)
SC.1=matrix(scan("00.data/sc.homom.homon.ABC.stat.txt"    ), byrow=T, ncol=ncol)
SC.2=matrix(scan("00.data/sc.heterom.homon.ABC.stat.txt"  ), byrow=T, ncol=ncol)
SC.3=matrix(scan("00.data/sc.homom.heteron.ABC.stat.txt"  ), byrow=T, ncol=ncol)
SC.4=matrix(scan("00.data/sc.heterom.heteron.ABC.stat.txt"), byrow=T, ncol=ncol)
SI.1=matrix(scan("00.data/si.homom.homon.ABC.stat.txt"    ), byrow=T, ncol=ncol)
SI.3=matrix(scan("00.data/si.homom.heteron.ABC.stat.txt"  ), byrow=T, ncol=ncol)

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

nlinesFul=min(nrow(IM.1), nrow(IM.2), nrow(IM.3), nrow(IM.4),
              nrow(SI.1), nrow(SI.3),
              nrow(AM.1), nrow(AM.2), nrow(AM.3), nrow(AM.4),
              nrow(SC.1), nrow(SC.2), nrow(SC.3), nrow(SC.4))

M_IM1b=IM.1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_IM2b=IM.2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_IM3b=IM.3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_IM4b=IM.4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_SC1b=SC.1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_SC2b=SC.2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_SC3b=SC.3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_SC4b=SC.4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_AM1b=AM.1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_AM2b=AM.2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_AM3b=AM.3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_AM4b=AM.4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_SI1b=SI.1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
M_SI3b=SI.3[c(1:nlinesFul),-c(1:3,12,13,18:25)]

target=target[-c(1:3,12,13,18:25)]
obs=matrix(rep(target[1:19],20), byrow=T, nrow=20)

#Comparaison All:
x=as.factor(c(rep(1:14,each=nlinesFul))

z=rbind(M_SI1b,M_SI3b,M_IM1b,M_IM2b,M_IM3b,M_IM4b,M_AM1b,M_AM2b,M_AM3b,M_AM4b,M_SC1b,M_SC2b,M_SC3b,M_SC4b)
res5=model_selection_abc_nnet(target=obs, x=x, sumstat=z, tol=3500/(length(x)), noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output="ALL.MODELS.29.06.16")

mean_all=apply(res5,2,mean)
sd_all=apply(res5,2,sd)
write.table(mean_all,"mean_all",quote=F,row.names=c("SI1","SI2","AM1","AM2","AM3","AM4","IM1","IM2","IM3","IM4","SC1","SC2","SC3","SC4"),col.names=F)
