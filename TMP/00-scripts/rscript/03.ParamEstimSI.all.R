#ParamEstimSC.R
source('/home/quentin/script/rscript/cv4abc.R')
if("abc" %in% rownames(installed.packages()) == FALSE) {install.packages("abc") }
library(abc)
target=as.numeric(read.table("OBS.ABC.stat.txt",skip=2,h=F))
target=target[-c(1:3,12,13,18:25)]
obs=target[1:19]

nlinesFul=as.numeric(strsplit(system("wc -l si.homom.heteron.ABC.stat.txt ", intern=T), " ")[[1]][1])
M_SC3=matrix(scan("si.homom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 

means1 <- colMeans(M_SC3, na.rm=TRUE)
for (j in 1:ncol(M_SC3)){
     M_SC3[is.na(M_SC3[, j]), j] <- means1[j]
 }
M_SC3b=M_SC3[,-c(1:3,12,13,18:25)]

#Hetero Hetero
priorfile=matrix(scan("si.homom.heteron.priorfile.txt"), byrow=T, nrow=nlinesFul) 
colnames(priorfile)=c("N1","N2","Na",      
						"Tsplit",
						"shape1Ne",        
						"shape2Ne",        
						"propNtrlNe1",     
						"propNtrlNe2"    
						)

abc.postSC <- abc(target  = obs, param = priorfile, sumstat = M_SC3b, tol = 1000/1e6, transf=c(rep("logit",8) ), 
	logit.bounds = rbind(range(priorfile[, 1]),range(priorfile[, 2]), range(priorfile[, 3]), range(priorfile[, 4]), range(priorfile[, 5]),range(priorfile[, 6]),range(priorfile[, 7]),range(priorfile[, 8])),
		hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

write.table(abc.postSC$weights,"weight.si.homom.heteron",quote=F,row.names=F,col.names=F)
write.table(abc.postSC$adj.values,"adjval.si.homom.heteron",quote=F,row.names=F)

z<-summary(abc.postSC)
write.table(z,"posterior.si.tol.0.01",quote=F,col.names=T,row.names=F)

pdf(file="si.homom.hetn.tol0.01.pdf",width=12,height=8)
par(mfrow=c(2,2))
prior.tetha1     <- density( (priorfile[,1] ) )
posterior.tetha1 <- density( abc.postSC$adj.values[,1] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta1 <- expression(theta[1])
#theta0 <- expression(theta[0])
plot( prior.tetha1,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta1,
      lwd=2, 
      main="theta1/thetaRef")  #or main= "Effective mutation rate (lf)")
lines(posterior.tetha1, col="blue",lwd=2)

legend(
 legend = c("prior distribution","posterior"),
col = c("black","blue"),
lty = c(1,1), lwd = c(1,1),
x = "topright",
cex = 1,
bty ="n")

prior.tetha2     <- density( priorfile[,2] )
posterior.tetha2 <- density( abc.postSC$adj.values[,2] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta <- expression(theta[2])
plot( prior.tetha2,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta,
      lwd=2, 
      main="theta2/thetaRef") #main="Effective mutation rate (lp)")
lines(posterior.tetha2, col="blue",lwd=2)
text(-1,8,"IM",cex=1)

prior.tethaA     <- density( priorfile[,3] )
posterior.tethaA <- density( abc.postSC$adj.values[,3] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta <- expression(theta[A])
plot( prior.tethaA,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta,
      lwd=2, 
      main="thetaA/thetaRef") # main="Effective mutation rate (ancestral population)")
lines(posterior.tethaA, col="blue",lwd=2)

prior.T     <- density( priorfile[,4] )
posterior.T <- density( abc.postSC$adj.values[,4] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=print_tau,
      lwd=2, 
      main="Divergence time (4Ngenerations)")
lines(posterior.T, col="blue",lwd=2)

dev.off()

###Homo Homom
M_SC3=matrix(scan("si.homom.homon.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 

means1 <- colMeans(M_SC3, na.rm=TRUE)
for (j in 1:ncol(M_SC3)){
     M_SC3[is.na(M_SC3[, j]), j] <- means1[j]
 }
M_SC3b=M_SC3[,-c(1:3,12,13,18:25)]

#Hetero Hetero
priorfile=matrix(scan("si.homom.homon.priorfile.txt"), byrow=T, nrow=nlinesFul) 
colnames(priorfile)=c("N1","N2","Na","Tsplit")

abc.postSC <- abc(target  = obs, param = priorfile, sumstat = M_SC3b, tol = 1000/1e6, transf=c(rep("logit",4) ), 
	logit.bounds = rbind(range(priorfile[, 1]),range(priorfile[, 2]), range(priorfile[, 3]), range(priorfile[, 4])),
		hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

write.table(abc.postSC$weights,"weight.si.m3.homom.homon",quote=F,row.names=F,col.names=F)
write.table(abc.postSC$adj.values,"adjval.si.m3.homom.homon",quote=F,row.names=F)

z<-summary(abc.postSC)
write.table(z,"posterior.si.m3.tol.0.01",quote=F,col.names=T,row.names=F)

pdf(file="si.m3.homom.homon.tol0.01.pdf",width=12,height=8)
par(mfrow=c(2,2))
prior.tetha1     <- density( (priorfile[,1] ) )
posterior.tetha1 <- density( abc.postSC$adj.values[,1] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta1 <- expression(theta[1])
#theta0 <- expression(theta[0])
plot( prior.tetha1,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta1,
      lwd=2, 
      main="theta1/thetaRef")  #or main= "Effective mutation rate (lf)")
lines(posterior.tetha1, col="blue",lwd=2)

legend(
 legend = c("prior distribution","posterior"),
col = c("black","blue"),
lty = c(1,1), lwd = c(1,1),
x = "topright",
cex = 1,
bty ="n")

prior.tetha2     <- density( priorfile[,2] )
posterior.tetha2 <- density( abc.postSC$adj.values[,2] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta <- expression(theta[2])
plot( prior.tetha2,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta,
      lwd=2, 
      main="theta2/thetaRef") #main="Effective mutation rate (lp)")
lines(posterior.tetha2, col="blue",lwd=2)
text(-1,8,"IM",cex=1)

prior.tethaA     <- density( priorfile[,3] )
posterior.tethaA <- density( abc.postSC$adj.values[,3] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta <- expression(theta[A])
plot( prior.tethaA,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta,
      lwd=2, 
      main="thetaA/thetaRef") # main="Effective mutation rate (ancestral population)")
lines(posterior.tethaA, col="blue",lwd=2)

prior.T     <- density( priorfile[,4] )
posterior.T <- density( abc.postSC$adj.values[,4] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=print_tau,
      lwd=2, 
      main="Divergence time (4Ngenerations)")
lines(posterior.T, col="blue",lwd=2)

dev.off()

