
library(abc)
library(weights)
library(laeken)
library(grid)
library(hexbin)


target <- as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))
obs <- target[-c(1:3,12,13,18:25)]
colon_count = paste(" awk -F ' ' '{print NF }' ", "00.data/im.homom.homon.ABC.stat.txt" , " |head -1 " )
ncol = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])

M_IM1=matrix(scan("00.data/im.homom.homon.ABC.stat.txt"  ), byrow=T, ncol=ncol)
means2 <- colMeans(M_IM1, na.rm=TRUE)
for (j in 1:ncol(M_IM1)){
             M_IM1[is.na(M_IM1[, j]), j] <- means2[j]
 }
nlinesFul=nrow(M_IM1)
M_IM1b=M_IM1[c(1:nlinesFul),-c(1:3,12,13,18:25)]

colon_count_pri = paste(" awk -F ' ' '{print NF }' ", "00.data/im.homom.homon.priorfile.txt" , " |head -1 " )
ncol_prior = as.numeric(strsplit(system ( colon_count_pri , intern=T), " ")[[1]] [1])

priorfile=matrix(scan("00.data/im.homom.homon.priorfile.txt"), byrow=T, ncol=ncol_prior)
colnames(priorfile)=c("N1","N2","Na","Tsplit","M1","M2")
priorfile <- priorfile[1:nlinesFul,]

model <- M_IM1b

logit.bound <- NULL
log.bound <- NULL
for (i in 1:ncol(priorfile))
{
log.bound <-range(priorfile[,i])
logit.bound <- rbind(logit.bound, log.bound)

}

abc.posterior <- abc(target  = obs, param = priorfile, sumstat = model, 
              tol = 1000/1e6, transf=c(rep("logit",ncol(priorfile) ) ),
              logit.bounds = logit.bound,
              hcorr=T, method  = "neuralnet", numnet=50, sizenet=15
               )
#export
write.table(summary(abc.posterior),"im.1.parameter.estim.txt", quote=F,row.names=F,col.names=F)
write.table(abc.posterior$weights,"im.1.weight.homo.homo",quote=F,row.names=F,col.names=F)
write.table(abc.posterior$adj.values,"im.1.adjval.homo.homo",quote=F,row.names=F)
pop <- strsplit(getwd(),"03.abc/")[[1]][2]

#One possible graph:
pdf(file=paste(pop,"im.1.pior.posterior.1000post.pdf", sep="."),width=12,height=8)
par(mfrow=c(2,3))

hist(priorfile[,1],
     breaks=seq(0,max(priorfile[,2]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,0.5),
     main=expression(theta["daughter pop1"]),cex.main=1.5,
     xlab=expression((theta[1])),
     ylab="probability density")
#hist((abc.posterior$unadj.values[,1] ),  breaks=seq(0,max(priorfile[,2]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist((abc.posterior$adj.values[,1]), breaks=seq(0,max(priorfile[,2]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights)
box()

legend(
 legend = c("prior distribution", "adjusted posterior"),
col = c("grey",rgb(0.8,0,0.3,0.2)),
pch = c(15,15),
x = "topleft",
cex = 1,
bty ="n")
box()

#legend(
# legend = c("prior distribution","unadjust posterior (rejection)", "adjusted posterior(neural net)"),
#col = c("grey",rgb(1,0,0,0.5),rgb(1,0.5,0,0.7)),
#pch = c(15,15,15),
#x = "topleft",
#cex = 1,
#bty ="n")
#box()
hist(priorfile[,2],
     breaks=seq(0,max(priorfile[,2]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,0.5),
     main=expression(theta["daugther pop2"]),cex.main=1.5,
     xlab=expression((theta[2])),
     ylab="probability density")
#hist((abc.posterior$unadj.values[,2] ),  breaks=seq(0,max(priorfile[,3]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
#wtd.hist((abc.posterior$adj.values[,2]), breaks=seq(0,max(priorfile[,2]+0.1),0.5), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc.posterior$weights)
wtd.hist((abc.posterior$adj.values[,2]), breaks=seq(0,max(priorfile[,2]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights)
box()

hist(priorfile[,3],
     breaks=seq(0,max(priorfile[,3]+0.1),0.5),
     col="grey",
     freq=F,
     ylim=c(0,0.5),
     main=expression(theta["Ancestral pop"]),cex.main=1.5,
     xlab=expression((theta["Ancestral pop"])),
     ylab="probability density")
#hist(abc.posterior$unadj.values[,3],   breaks=seq(0,max(priorfile[,4]+0.1),0.1), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,3], breaks=seq(0,max(priorfile[,3]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()

hist(priorfile[,4],
	breaks=seq(0,max(priorfile[,4]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,0.5),
	main="divergence time( 4N generations)",cex.main=1.5,
	xlab= expression(tau),
	ylab="probability density")
#hist(abc.posterior$unadj.values[,4],   breaks=seq(0,max(priorfile[,5]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,4], breaks=seq(0,max(priorfile[,4]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
#
hist(priorfile[,5],
	breaks=seq(0,max(priorfile[,5]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,0.25),
	main="Effective migration rate ",cex.main=1.5,
	xlab= "M1<-2",
	ylab="probab density")
#hist(abc.posterior$unadj.values[,5],   breaks=seq(0,max(priorfile[,5]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,5], breaks=seq(0,max(priorfile[,5]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
                  
hist(priorfile[,6],
        breaks=seq(0,max(priorfile[,6]+0.1),1),
        col="lightgrey",
        freq=F,
        ylim=c(0,0.25),
        main="Effective migration rate ",cex.main=1.5,
        xlab= "M2<-1" ,
        ylab="probab density")
#hist(abc.posterior$unadj.values[,6],   breaks=seq(0,max(priorfile[,7]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,6], breaks=seq(0,max(priorfile[,6]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))

box()
dev.off()

