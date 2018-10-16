
#if("abc" %in% rownames(installed.packages()) == FALSE) {install.packages("abc") }
#if("weigths" %in% rownames(installed.packages()) == FALSE) {install.packages("weights") }
#if("laeken" %in% rownames(installed.packages()) == FALSE) {install.packages("laeken") }
#if("grid" %in% rownames(installed.packages()) == FALSE) {install.packages("grid") }
#if("hexbin" %in% rownames(installed.packages()) == FALSE) {install.packages("hexbin") }

library(abc)
library(weights)
library(laeken)
library(grid)
library(hexbin)

target=as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))
target=target[-c(1:3,12,13,18:25)]
obs=target

colon_count = paste(" awk -F ' ' '{print NF }' ", "00.data/sc.homom.homon.ABC.stat.txt" , " |head -1 " )
ncol = as.numeric(strsplit(system ( colon_count , intern=T), " ")[[1]] [1])

#load simulations
SC.1=matrix(scan("00.data/sc.homom.homon.ABC.stat.txt"    ), byrow=T, ncol=ncol)

# Replace missing data
for(i in 1:ncol(SC.1)){
  SC.1[which(SC.1[,i]=="NaN"),i]=mean(SC.1[,i], na.rm=T)
}

nlinesFul=nrow(SC.1)

M_SC1b=SC.1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
model <- M_SC1b
#Priorfile
colon_count_pri = paste(" awk -F ' ' '{print NF }' ", "00.data/sc.homom.homon.priorfile.txt" , " |head -1 " )
ncol_prior = as.numeric(strsplit(system ( colon_count_pri , intern=T), " ")[[1]] [1])

priorfile=matrix(scan("00.data/sc.homom.homon.priorfile.txt"), byrow=T, ncol=ncol_prior)
colnames(priorfile)=c("N1","N2","Na","Tsplit","Tsc","M1","M2")
priorfile <- priorfile[c(1:nlinesFul),]

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

write.table(summary(abc.posterior),"sc.1.parameter.estim.txt", quote=F,row.names=F,col.names=F)
write.table(abc.posterior$weights,"sc.1.weight.homo.homo",quote=F,row.names=F,col.names=F)
write.table(abc.posterior$adj.values,"sc.1.adjval.homo.homo",quote=F,row.names=F)
pop <- strsplit(getwd(),"03.abc/")[[1]][2]

pdf(file=paste(pop,"sc.1.pior.posterior.1000post.pdf",sep=".") ,width=12,height=8)
par(mfrow=c(3,3))

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
	main=paste(expression(tau),"=Tsplit/4Nref ( generations)", sep=" "),cex.main=1.5,
	xlab= expression(tau),
	ylab="probability density")
#hist(abc.posterior$unadj.values[,4],   breaks=seq(0,max(priorfile[,5]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,4], breaks=seq(0,max(priorfile[,4]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
hist(priorfile[,5],
	breaks=seq(0,max(priorfile[,5]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,0.5),
	main="Time of Secondary Contact ( 4N generations)",cex.main=1.5,
	xlab= expression(tau),
	ylab="probability density")
#hist(abc.posterior$unadj.values[,4],   breaks=seq(0,max(priorfile[,5]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,5], breaks=seq(0,max(priorfile[,5]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
#
hist(priorfile[,6],
	breaks=seq(0,max(priorfile[,6]+0.1),1),
	col="lightgrey",
	freq=F,
	ylim=c(0,0.25),
	main="Effective migration rate (4Nm) ",cex.main=1.5,
	xlab= "M1<-2",
	ylab="probab density")
#hist(abc.posterior$unadj.values[,5],   breaks=seq(0,max(priorfile[,5]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,6], breaks=seq(0,max(priorfile[,6]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()
                  
hist(priorfile[,7],
        breaks=seq(0,max(priorfile[,7]+0.1),1),
        col="lightgrey",
        freq=F,
        ylim=c(0,0.25),
        main="Effective migration rate ",cex.main=1.5,
        xlab= "M2<-1" ,
        ylab="probab density")
#hist(abc.posterior$unadj.values[,6],   breaks=seq(0,max(priorfile[,7]+0.1),0.5), col=rgb(1,0,0,0.5), freq=F, add=T, )
wtd.hist(abc.posterior$adj.values[,7], breaks=seq(0,max(priorfile[,7]+0.1),0.5), col=rgb(0.8,0,0.3,0.2), freq=F, add=T, weight=abc.posterior$weights/sum(abc.posterior$weights))
box()

dev.off()
