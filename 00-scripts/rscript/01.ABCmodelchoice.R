source('00-scripts/rscript/cv4abc.R')

colon_count = paste(" awk -F ' ' '{print NF }' ", "00.data/im.homom.homon.ABC.stat.txt" , " |head -1 " )
ncol = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])

target=as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))

#load simulations
IM1=matrix(scan("00.data/im.homom.homon.ABC.stat.txt"    ), byrow=T, ncol=ncol)
IM2=matrix(scan("00.data/im.heterom.homon.ABC.stat.txt"  ), byrow=T, ncol=ncol)
IM3=matrix(scan("00.data/im.homom.heteron.ABC.stat.txt"  ), byrow=T, ncol=ncol)
IM4=matrix(scan("00.data/im.heterom.heteron.ABC.stat.txt"), byrow=T, ncol=ncol)
AM1=matrix(scan("00.data/am.homom.homon.ABC.stat.txt"    ), byrow=T, ncol=ncol)
AM2=matrix(scan("00.data/am.heterom.homon.ABC.stat.txt"  ), byrow=T, ncol=ncol)
AM3=matrix(scan("00.data/am.homom.heteron.ABC.stat.txt"  ), byrow=T, ncol=ncol)
AM4=matrix(scan("00.data/am.heterom.heteron.ABC.stat.txt"), byrow=T, ncol=ncol)
SC1=matrix(scan("00.data/sc.homom.homon.ABC.stat.txt"    ), byrow=T, ncol=ncol)
SC2=matrix(scan("00.data/sc.heterom.homon.ABC.stat.txt"  ), byrow=T, ncol=ncol)
SC3=matrix(scan("00.data/sc.homom.heteron.ABC.stat.txt"  ), byrow=T, ncol=ncol)
SC4=matrix(scan("00.data/sc.heterom.heteron.ABC.stat.txt"), byrow=T, ncol=ncol)
SI1=matrix(scan("00.data/si.homom.homon.ABC.stat.txt"    ), byrow=T, ncol=ncol)
SI3=matrix(scan("00.data/si.homom.heteron.ABC.stat.txt"  ), byrow=T, ncol=ncol)

colSums(is.na(IM1))
means2 <- colMeans(IM1, na.rm=TRUE)
for (j in 1:ncol(IM1)){
     IM1[is.na(IM1[, j]), j] <- means2[j]
 }
colSums(is.na(IM2))
means2 <- colMeans(IM2, na.rm=TRUE)
for (j in 1:ncol(IM2)){
     IM2[is.na(IM2[, j]), j] <- means2[j]
 }
colSums(is.na(IM3))
means1 <- colMeans(IM3, na.rm=TRUE)
for (j in 1:ncol(IM3)){
     IM3[is.na(IM3[, j]), j] <- means1[j]
 }
colSums(is.na(IM4))
means3 <- colMeans(IM4, na.rm=TRUE)
for (j in 1:ncol(IM4)){
     IM4[is.na(IM4[, j]), j] <- means3[j]
 }

means2 <- colMeans(AM1, na.rm=TRUE)
for (j in 1:ncol(AM1)){
     AM1[is.na(AM1[, j]), j] <- means2[j]
 }

means2 <- colMeans(AM2, na.rm=TRUE)
for (j in 1:ncol(AM2)){
     AM2[is.na(AM2[, j]), j] <- means2[j]
 }

means1 <- colMeans(AM3, na.rm=TRUE)
for (j in 1:ncol(AM3)){
     AM3[is.na(AM3[, j]), j] <- means1[j]
 }

means3 <- colMeans(AM4, na.rm=TRUE)
for (j in 1:ncol(AM4)){
     AM4[is.na(AM4[, j]), j] <- means3[j]
 }

means2 <- colMeans(SC1, na.rm=TRUE)
for (j in 1:ncol(SC1)){
     SC1[is.na(SC1[, j]), j] <- means2[j]
 }

means2 <- colMeans(SC2, na.rm=TRUE)
for (j in 1:ncol(SC2)){
     SC2[is.na(SC2[, j]), j] <- means2[j]
 }
means1 <- colMeans(SC3, na.rm=TRUE)
for (j in 1:ncol(SC3)){
     SC3[is.na(SC3[, j]), j] <- means1[j]
 }
means3 <- colMeans(SC4, na.rm=TRUE)
for (j in 1:ncol(SC4)){
     SC4[is.na(SC4[, j]), j] <- means3[j]
 }

means2 <- colMeans(SI1, na.rm=TRUE)
for (j in 1:ncol(SI1)){
     SI1[is.na(SI1[, j]), j] <- means2[j]
 }
means1 <- colMeans(SI3, na.rm=TRUE)
for (j in 1:ncol(SI3)){
     SI3[is.na(SI3[, j]), j] <- means1[j]
 }

nlinesFul=min(nrow(IM1), nrow(IM2), nrow(IM3), nrow(IM4),
              nrow(SI1), nrow(SI3),
              nrow(AM1), nrow(AM2), nrow(AM3), nrow(AM4),
              nrow(SC1), nrow(SC2), nrow(SC3), nrow(SC4))

IM1b=IM1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
IM2b=IM2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
IM3b=IM3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
IM4b=IM4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC1b=SC1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC2b=SC2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC3b=SC3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC4b=SC4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM1b=AM1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM2b=AM2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM3b=AM3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM4b=AM4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SI1b=SI1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SI3b=SI3[c(1:nlinesFul),-c(1:3,12,13,18:25)]

target=target[-c(1:3,12,13,18:25)]
obs=matrix(rep(target[1:19],20), byrow=T, nrow=20)

#Comparaison All:
x=as.factor(c(rep(1:14,each=nlinesFul)))

z=rbind(SI1b,SI3b,
        IM1b,IM2b,IM3b,IM4b,
        AM1b,AM2b,AM3b,AM4b,
        SC1b,SC2b,SC3b,SC4b)

res5=model_selection_abc_nnet(target=obs, 
                              x=x, 
                              sumstat=z, 
                              tol=0.00025,
                              noweight=F,
                              rejmethod=F,
                              nb.nnet=50,
                              size.nnet=15,
                              output="ALL.MODELS.29.06.16")

mean_all=apply(res5,2,mean)
sd_all=apply(res5,2,sd)
write.table(mean_all,
            "mean_all",
            quote=F,
            row.names=c("SI1","SI2","AM1","AM2","AM3","AM4","IM1","IM2","IM3","IM4","SC1","SC2","SC3","SC4"),
            col.names=F)
