#ParamEstimSC.R
source('/home/quentin/script/rscript/cv4abc.R')
if("abc" %in% rownames(installed.packages()) == FALSE) {install.packages("abc") }
library(abc)
target=as.numeric(read.table("OBS.ABC.stat.txt",skip=2,h=F))
target=target[-c(1:3,12,13,18:25)]
obs=target[1:19]

nlinesFul=as.numeric(strsplit(system("wc -l sc.homom.homon.ABC.stat.txt ", intern=T), " ")[[1]][1])
M_SC1=matrix(scan("sc.homom.homon.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 

means1 <- colMeans(M_SC1, na.rm=TRUE)
for (j in 1:ncol(M_SC1)){
     M_SC1[is.na(M_SC1[, j]), j] <- means1[j]
 }
M_SC1b=M_SC1[,-c(1:3,12,13,18:25)]

#Hetero Hetero
priorfile=matrix(scan("sc.homom.homon.priorfile.txt"), byrow=T, nrow=nlinesFul) 
colnames(priorfile)=c("N1","N2","Na","Tsplit","Tsc","M1","M2")

abc.postSC <- abc(target  = obs, param = priorfile, sumstat = M_SC1b, tol = 1000/1e6, transf=c(rep("logit",7) ), 
	logit.bounds = rbind(range(priorfile[, 1]),range(priorfile[, 2]), range(priorfile[, 3]), range(priorfile[, 4]), range(priorfile[, 5]),range(priorfile[, 6]),range(priorfile[, 7])),
	hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

write.table(abc.postSC$weights,"weight.sc.m1.homom.homon",quote=F,row.names=F,col.names=F)
write.table(abc.postSC$adj.values,"adjval.sc.m1.homom.homon",quote=F,row.names=F)
z<-summary(abc.postSC)
write.table(z,"posterior.sc.m1.tol.0.01",quote=F,col.names=T,row.names=F)

pdf(file="sc.m1.homo.homo.tol0.1.pdf",width=12,height=8)
par(mfrow=c(4,4))
prior.tetha1     <- density( (priorfile[,1] ) )
posterior.tetha1 <- density( abc.postSC$adj.values[,1] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta1 <- expression(theta[1])
#theta0 <- expression(theta[0])
plot( prior.tetha1,
      xlim=c(0,20),
      ylim=c(0,0.5), 
      ylab="probability density", 
      xlab=theta1,
      lwd=2, 
      main="theta1")  #or main= "Effective mutation rate (lf)")
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
      ylim=c(0,0.5), 
      ylab="probability density", 
      xlab=theta,
      lwd=2, 
      main="theta2") #main="Effective mutation rate (lp)")
lines(posterior.tetha2, col="blue",lwd=2)
text(-1,8,"IM",cex=1)

prior.tethaA     <- density( priorfile[,3] )
posterior.tethaA <- density( abc.postSC$adj.values[,3] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta <- expression(theta[A])
plot( prior.tethaA,
      xlim=c(0,20),
      ylim=c(0,0.5), 
      ylab="probability density", 
      xlab=theta,
      lwd=2, 
      main="thetaAncestral") # main="Effective mutation rate (ancestral population)")
lines(posterior.tethaA, col="blue",lwd=2)

prior.T     <- density( priorfile[,4] )
posterior.T <- density( abc.postSC$adj.values[,4] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,0.5), 
      ylab="probability density", 
      xlab=print_tau,
      lwd=2, 
      main="Divergence time (4Ngenerations)")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,5] )
posterior.T <- density( abc.postSC$adj.values[,5] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="Tsc",
      lwd=2, 
      main="Time of SC")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,6] )
posterior.T <- density( abc.postSC$adj.values[,6] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.07), 
      ylab="probability density", 
      xlab="M1",
      lwd=2, 
      main="M1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,7] )
posterior.T <- density( abc.postSC$adj.values[,7] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.07), 
      ylab="probability density", 
      xlab="M2",
      lwd=2, 
      main="M2")
lines(posterior.T, col="blue",lwd=2)

dev.off()

################################################## Hetero N Hetero M##################################################
nlinesFul=as.numeric(strsplit(system("wc -l sc.heterom.heteron.ABC.stat.txt ", intern=T), " ")[[1]][1]) 
M_SC4=matrix(scan("sc.heterom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 

means3 <- colMeans(M_SC4, na.rm=TRUE)
for (j in 1:ncol(M_SC4)){
     M_SC4[is.na(M_SC4[, j]), j] <- means3[j]
 }

M_SC4b=M_SC4[,-c(1:3,12,13,18:25)]
#Hetero Hetero
priorfile=matrix(scan("sc.heterom.heteron.priorfile.txt"), byrow=T, nrow=nlinesFul) 
colnames(priorfile)=c("N1","N2","Na", "Tsplit","Tsc","shape1Ne","shape2Ne",        
						"propNtrlNe1","propNtrlNe2","M1","M2",
						"shape1M1",
						"shape2M1",      
						"shape1M2",      
						"shape2M2",      
						"propNtrlM1",     
						"propNtrlM2")

abc.postSC <- abc(target  = obs, param = priorfile, sumstat = M_SC4b, tol = 1000/1e6, transf=c(rep("logit",17) ), 
	logit.bounds = rbind(range(priorfile[, 1]),range(priorfile[, 2]), range(priorfile[, 3]), range(priorfile[, 4]), range(priorfile[, 5]),range(priorfile[, 6]),range(priorfile[, 7]),range(priorfile[, 8]),range(priorfile[, 9]),
	range(priorfile[, 10]), range(priorfile[, 11]), range(priorfile[, 12]), range(priorfile[, 13]), range(priorfile[, 14]), range(priorfile[, 15]), range(priorfile[, 16]), range(priorfile[, 17])),
	hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

write.table(abc.postSC$weights,"weight.sc.m4.heterom.heteron",quote=F,row.names=F,col.names=F)
write.table(abc.postSC$adj.values,"adjval.sc.m4.heteroom.heteron",quote=F,row.names=F)
z<-summary(abc.postSC)
write.table(z,"posterior.sc.m4.tol.0.01",quote=F,col.names=T,row.names=F)

pdf(file="sc.m4.het.het.tol0.1.pdf",width=12,height=8)
par(mfrow=c(4,4))
prior.tetha1     <- density( (priorfile[,1] ) )
posterior.tetha1 <- density( abc.postSC$adj.values[,1] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta1 <- expression(theta[1])
#theta0 <- expression(theta[0])
plot( prior.tetha1,
      xlim=c(0,20),
      ylim=c(0,0.5), 
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
      ylim=c(0,0.5), 
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
      ylim=c(0,0.5), 
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
      ylim=c(0,0.5), 
      ylab="probability density", 
      xlab=print_tau,
      lwd=2, 
      main="Divergence time (4Ngenerations)")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,5] )
posterior.T <- density( abc.postSC$adj.values[,5] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="Tsc",
      lwd=2, 
      main="Time of SC")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,10] )
posterior.T <- density( abc.postSC$adj.values[,10] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.07), 
      ylab="probability density", 
      xlab="M1",
      lwd=2, 
      main="M1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,11] )
posterior.T <- density( abc.postSC$adj.values[,11] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.07), 
      ylab="probability density", 
      xlab="M2",
      lwd=2, 
      main="M2")
lines(posterior.T, col="blue",lwd=2)

#						Homo M Hetero N
########################################################################################################################################
nlinesFul=as.numeric(strsplit(system("wc -l sc.homom.heteron.ABC.stat.txt ", intern=T), " ")[[1]][1])
M_SC3=matrix(scan("sc.homom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 

means1 <- colMeans(M_SC3, na.rm=TRUE)
for (j in 1:ncol(M_SC3)){
     M_SC3[is.na(M_SC3[, j]), j] <- means1[j]
 }
M_SC3b=M_SC3[,-c(1:3,12,13,18:25)]

#Hetero Hetero
priorfile=matrix(scan("sc.homom.heteron.priorfile.txt"), byrow=T, nrow=nlinesFul) 
colnames(priorfile)=c("N1","N2","Na",      
						"Tsplit",
						"Tsc",  
						"shape1Ne",        
						"shape2Ne",        
						"propNtrlNe1",     
						"propNtrlNe2",    
						"M1",
						"M2")

abc.postSC <- abc(target  = obs, param = priorfile, sumstat = M_SC3b, tol = 1000/1e6, transf=c(rep("logit",11) ), 
	logit.bounds = rbind(range(priorfile[, 1]),range(priorfile[, 2]), range(priorfile[, 3]), range(priorfile[, 4]), range(priorfile[, 5]),range(priorfile[, 6]),range(priorfile[, 7]),range(priorfile[, 8]),range(priorfile[, 9]),
	range(priorfile[, 10]), range(priorfile[, 11])),
	hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

write.table(abc.postSC$weights,"weight.sc.m3.homom.heteron",quote=F,row.names=F,col.names=F)
write.table(abc.postSC$adj.values,"adjval.sc.m3.homom.heteron",quote=F,row.names=F)
z<-summary(abc.postSC)
write.table(z,"posterior.sc.m3.homom.heteron.tol.0.01",quote=F,col.names=T,row.names=F)

pdf(file="sc.m3.homom.heteron.tol0.1.pdf",width=12,height=8)
par(mfrow=c(4,4))
prior.tetha1     <- density( (priorfile[,1] ) )
posterior.tetha1 <- density( abc.postSC$adj.values[,1] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta1 <- expression(theta[1])
#theta0 <- expression(theta[0])
plot( prior.tetha1,
      xlim=c(0,20),
      ylim=c(0,0.5), 
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
      ylim=c(0,0.5), 
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
      ylim=c(0,0.5), 
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
      ylim=c(0,0.5), 
      ylab="probability density", 
      xlab=print_tau,
      lwd=2, 
      main="Divergence time (4Ngenerations)")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,5] )
posterior.T <- density( abc.postSC$adj.values[,5] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="Tsc",
      lwd=2, 
      main="Time of SC")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,10] )
posterior.T <- density( abc.postSC$adj.values[,10] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.07), 
      ylab="probability density", 
      xlab="M1",
      lwd=2, 
      main="M1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,11] )
posterior.T <- density( abc.postSC$adj.values[,11] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.07), 
      ylab="probability density", 
      xlab="M2",
      lwd=2, 
      main="M2")
lines(posterior.T, col="blue",lwd=2)

dev.off()

#						Homo N Hetero M
########################################################################################################################################
M_SC1=matrix(scan("sc.heterom.homon.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 

means1 <- colMeans(M_SC1, na.rm=TRUE)
for (j in 1:ncol(M_SC1)){
     M_SC1[is.na(M_SC1[, j]), j] <- means1[j]
 }
M_SC1b=M_SC1[,-c(1:3,12,13,18:25)]

priorfile=matrix(scan("sc.heterom.homon.priorfile.txt"), byrow=T, nrow=nlinesFul) 
colnames(priorfile)=c("N1","N2","Na",      
						"Tsplit",
						"Tsc",  
						"M1",
						"M2",
						"shape1M1",
						"shape2M1",      
						"shape1M2",      
						"shape2M2",      
						"propNtrlM1",     
						"propNtrlM2")

abc.postSC <- abc(target  = obs, param = priorfile, sumstat = M_SC1b, tol = 1000/1e6, transf=c(rep("logit",13) ), 
	logit.bounds = rbind(range(priorfile[, 1]),range(priorfile[, 2]), range(priorfile[, 3]), range(priorfile[, 4]), range(priorfile[, 5]),range(priorfile[, 6]),range(priorfile[, 7]),range(priorfile[, 8]),range(priorfile[, 9]),
	range(priorfile[, 10]), range(priorfile[, 11]), range(priorfile[, 12]), range(priorfile[, 13])),
	hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

write.table(abc.postSC$weights,"weight.sc.m2.homom.homon",quote=F,row.names=F,col.names=F)
write.table(abc.postSC$adj.values,"adjval.sc.m2.homom.homon",quote=F,row.names=F)
z<-summary(abc.postSC)
write.table(z,"posterior.sc.m2.tol.0.01",quote=F,col.names=T,row.names=F)

##################################################"
pdf(file="sc.m2.heterom.homon.tol0.1.pdf",width=12,height=8)
par(mfrow=c(4,4))
prior.tetha1     <- density( (priorfile[,1] ) )
posterior.tetha1 <- density( abc.postSC$adj.values[,1] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
theta1 <- expression(theta[1])
#theta0 <- expression(theta[0])
plot( prior.tetha1,
      xlim=c(0,20),
      ylim=c(0,0.5), 
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
      ylim=c(0,0.5), 
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
      ylim=c(0,0.5), 
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
      ylim=c(0,0.5), 
      ylab="probability density", 
      xlab=print_tau,
      lwd=2, 
      main="Divergence time (4Ngenerations)")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,5] )
posterior.T <- density( abc.postSC$adj.values[,5] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="Tsc",
      lwd=2, 
      main="Time of SC")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,6] )
posterior.T <- density( abc.postSC$adj.values[,6] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.07), 
      ylab="probability density", 
      xlab="M1",
      lwd=2, 
      main="M1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,7] )
posterior.T <- density( abc.postSC$adj.values[,7] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.07), 
      ylab="probability density", 
      xlab="M2",
      lwd=2, 
      main="M2")
lines(posterior.T, col="blue",lwd=2)

dev.off()
