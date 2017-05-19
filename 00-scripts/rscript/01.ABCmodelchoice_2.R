source('/home/quentin/script/rscript/cv4abc.R')

target=as.numeric(read.table("OBS.ABC.stat.txt",skip=2,h=F))

nlinesFul=1e6

M_PA1 <-matrix(scan("pan.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul) #les 1e6 simuls sont mergés ensemble au préalable avec du bash et sed.
M_IM1 <-matrix(scan("im.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul) #les 1e6 simuls sont mergés ensemble au préalable avec du bash et sed.
M_BO1 <-matrix(scan("bot.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul) #les 1e6 sSIuls sont mergés ensemble au préalable avec du bash et sed.

colSums(is.na(M_PA1))
means2 <- colMeans(M_PA1, na.rm=TRUE)
for (j in 1:ncol(M_PA1)){
     M_PA1[is.na(M_PA1[, j]), j] <- means2[j]
 }

colSums(is.na(M_IM1))
means2 <- colMeans(M_IM1, na.rm=TRUE)
for (j in 1:ncol(M_IM1)){
     M_IM1[is.na(M_IM1[, j]), j] <- means2[j]
 }
for (j in 1:ncol(M_BO1)){
     M_BO1[is.na(M_BO1[, j]), j] <- means2[j]
 }
M_PA1b=M_PA1[,-c(1:3,12,13,18:25)]
M_IM1b=M_IM1[,-c(1:3,12,13,18:25)]
M_BO1b=M_BO1[,-c(1:3,12,13,18:25)]

target=target[-c(1:3,12,13,18:25)]
obs=matrix(rep(target[1:19],10), byrow=T, nrow=10)

#Comparaison All:
x=as.factor(c(rep("A",nlinesFul),rep("B",nlinesFul),rep("C",nlinesFul))),

z=rbind(M_BO1b,M_IM1b,M_PA1b)
res5=model_selection_abc_nnet(target=obs, x=x, sumstat=z, tol=1000/(length(x)), noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output="ALL.MODELS.29.06.16")

mean_all=apply(res5,2,mean)
sd_all=apply(res5,2,sd)
write.table(mean_all,"mean_all",quote=F,row.names=c("BO1","IM1","PA1b"),col.names=F)
