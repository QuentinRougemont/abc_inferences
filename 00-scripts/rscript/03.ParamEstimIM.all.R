#ParamEstimIM.R
source('00-scripts/rscript/cv4abc.R')
if("abc" %in% rownames(installed.packages()) == FALSE)
    {install.packages("abc") }

library(abc)

argv<-commandArgs(TRUE)
model  <-argv[1] #identification of the model either 1,2,3,4
#m1 = homom homon
#m2 = heterom homon
#m3 = homom heteron
#m4 = heterom heteron
model <- 1
#load data
target=as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))
target=target[-c(1:3,12,13,18:25)]
obs=target[1:19]

nlinesFul=as.numeric(strsplit(system
    ("wc -l 00.data/im.homom.homon.ABC.stat.txt ", intern=T), " ")[[1]][1])
if(model==1){
M1=matrix(scan("00.data/im.homom.homon.ABC.stat.txt"), 
    byrow=T, nrow=nlinesFul) 
}else if (model == 2 ){
M1=matrix(scan("00.data/im.heterom.homon.ABC.stat.txt"), 
    byrow=T, nrow=nlinesFul) 
}else if (model == 3 ){
M1=matrix(scan("00.data/im.homom.heteron.ABC.stat.txt"), 
    byrow=T, nrow=nlinesFul) 
}else if (model == 4 ){
M1=matrix(scan("00.data/im.heterom.heteron.ABC.stat.txt"), 
    byrow=T, nrow=nlinesFul) 
} 
#replace NA by means
clean <-function(model){
    means <-colMeans(model, na.rm=T)
    for (j in 1:ncol(model)){
         model[is.na(model[, j]), j] <- means[j]
    }
}
f <- function(x){ 
m <- mean(x, na.rm = TRUE) 
x[is.na(x)] <- m 
x 
} 
M1 <- apply(M1, 2, f) 
M1b=M1[,-c(1:3,12,13,18:25)]

#read priorfile 
if(model==1){
priorfile=matrix(scan("00.data/im.homom.homon.priorfile.txt"), 
    byrow=T, nrow=nlinesFul) 
}else if (model == 2 ){
priorfile=matrix(scan("00.data/im.heterom.homon.priorfile.txt"), 
    byrow=T, nrow=nlinesFul) 
}else if (model == 3 ){
priorfile=matrix(scan("00.data/im.homom.heteron.priorfile.txt"), 
    byrow=T, nrow=nlinesFul) 
}else if (model == 4 ){
priorfile=matrix(scan("00.data/im.heterom.heteron.priorfile.txt"), 
    byrow=T, nrow=nlinesFul) 
} 

if(model==1){
colnames(priorfile)=c("N1","N2","Na","Tsplit","M1","M2")
}else if (model == 2 ){
colnames(priorfile)=c("N1","N2","Na", "Tsplit",
                    "M1", "M2", 
                    "shape1M1","shape2M1", "shape1M2", "shape2M2", 
                    "propNtrlM1","propNtrlM2")
}else if (model == 3 ){
colnames(priorfile)=c("N1","N2","Na", "Tsplit",
                    "shape1Ne", "shape2Ne", "propNtrlNe1", "propNtrlNe2", 
                    "M1", "M2")
}else if (model == 4 ){
colnames(priorfile)=c("N1","N2","Na", "Tsplit",
                    "shape1Ne", "shape2Ne", "propNtrlNe1","propNtrlNe2", 
                    "M1", "M2", 
                    "shape1M1","shape2M1", "shape1M2", "shape2M2", 
                    "propNtrlM1","propNtrlM2")
} 

log.bound=matrix(NA,nrow=ncol(priorfile), ncol=2)
for(i in 1:ncol(priorfile))
{
    log.bound[i,]=(range(priorfile[,i]))
}

#parameter estimates
abc.post <- abc(target  = obs, param = priorfile, 
    sumstat = M1b, 
    tol = 1000/1e6, 
    transf=c(rep("logit",ncol(priorfile) ) ), 
    logit.bounds = log.bound,
    hcorr=T, 
    method  = "neuralnet", 
    numnet=50, sizenet=15)

if(model==1){
sink("postadjval.im.homo.homo.tol0.1.txt")
print(abc.post$adj.values)
sink()
}else if (model == 2 ){
sink("postadjval.im.heterom.homon.tol0.1.txt")
print(abc.post$adj.values)
sink()
}else if (model == 3 ){
sink("postadjval.im.homom.heteron.tol0.1.txt")
print(abc.post$adj.values)
sink()
}else if (model == 4 ){
sink("postadjval.im.heteron.heteron.tol0.1.txt")
print(abc.post$adj.values)
sink()
} 

if(model==1){
sink("postadjval.im.homo.homo.tol0.1.txt")
print(summary(abc.post))
sink()
}else if (model == 2 ){
sink("postadjval.im.heterom.homon.tol0.1.txt")
print(summary(abc.post))
sink()
}else if (model == 3 ){
sink("postadjval.im.homom.heteron.tol0.1.txt")
print(summary(abc.post))
sink()
}else if (model == 4 ){
sink("postadjval.im.heteron.heteron.tol0.1.txt")
print(summary(abc.post))
sink()
} 

pdf(file="im.homo.homo.tol0.1.pdf",width=12,height=8)
par(mfrow=c(4,4))
prior.tetha1     <- density( (priorfile[,1] ) )
posterior.tetha1 <- density( abc.post$adj.values[,1] , 
    weights=abc.post$weights/sum(abc.post$weights) )
theta1 <- expression(theta[1])
#theta0 <- expression(theta[0])
plot( prior.tetha1,
      xlim=c(0,20),
      ylim=c(0,1), 
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
posterior.tetha2 <- density( abc.post$adj.values[,2] , 
    weights=abc.post$weights/sum(abc.post$weights) )
theta <- expression(theta[2])
plot( prior.tetha2,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta,
      lwd=2, 
      main="theta2") #main="Effective mutation rate (lp)")
lines(posterior.tetha2, col="blue",lwd=2)
text(-1,8,"IM",cex=1)

prior.tethaA     <- density( priorfile[,3] )
posterior.tethaA <- density( abc.post$adj.values[,3] , 
    weights=abc.post$weights/sum(abc.post$weights) )
theta <- expression(theta[A])
plot( prior.tethaA,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta,
      lwd=2, 
      main="thetaAncestral") # main="Effective mutation rate (ancestral population)")
lines(posterior.tethaA, col="blue",lwd=2)

prior.T     <- density( priorfile[,4] )
posterior.T <- density( abc.post$adj.values[,4] , 
    weights=abc.post$weights/sum(abc.post$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=print_tau,
      lwd=2, 
      main="Divergence time (4Ngenerations)")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,5] )
posterior.T <- density( abc.post$adj.values[,5] , 
    weights=abc.post$weights/sum(abc.post$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="M1",
      lwd=2, 
      main="M1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,6] )
posterior.T <- density( abc.post$adj.values[,6] , 
    weights=abc.post$weights/sum(abc.post$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.10), 
      ylab="probability density", 
      xlab="M2",
      lwd=2, 
      main="M2")
lines(posterior.T, col="blue",lwd=2)
dev.off()

################################################# Hetero N Hetero M##################################################
##################################################
######################################################################################################################
pdf(file="im.het.het.tol0.1.pdf",width=12,height=8)
par(mfrow=c(4,4))
prior.tetha1     <- density( (priorfile[,1] ) )
posterior.tetha1 <- density( abc.postSC$adj.values[,1] , 
    weights=abc.postSC$weights/sum(abc.postSC$weights) )
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

prior.T    <- density( priorfile[,5] )
posterior.T <- density( abc.postSC$adj.values[,5] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="M1",
      lwd=2, 
      main="M1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,10] )
posterior.T <- density( abc.postSC$adj.values[,10] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.10), 
      ylab="probability density", 
      xlab="M2",
      lwd=2, 
      main="M2")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,11] )
posterior.T <- density( abc.postSC$adj.values[,11] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.1), 
      ylab="probability density", 
      xlab="propNtrlNe1",
      lwd=2, 
      main="propNtrlNe1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,8] )
posterior.T <- density( abc.postSC$adj.values[,8] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,6000),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="propNtrlNe2",
      lwd=2, 
      main="propNtrlNe2")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,9] )
posterior.T <- density( abc.postSC$adj.values[,9] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,6000),
      ylim=c(0,0.20), 
      ylab="probability density", 
      xlab="shape1M1 ",
      lwd=2, 
      main="shape1M1  ")
lines(posterior.T, col="blue",lwd=2)


prior.T    <- density( priorfile[,12] )
posterior.T <- density( abc.postSC$adj.values[,12] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,20),
      ylim=c(0,0.2), 
      ylab="probability density", 
      xlab="shape2M1 ",
      lwd=2, 
      main="shape2M1 ")
lines(posterior.T, col="blue",lwd=2)


prior.T    <- density( priorfile[,13] )
posterior.T <- density( abc.postSC$adj.values[,13] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,20),
      ylim=c(0,0.2), 
      ylab="probability density", 
      xlab="shape1M2 ",
      lwd=2, 
      main="shape1M2 ")
lines(posterior.T, col="blue",lwd=2)


prior.T    <- density( priorfile[,14] )
posterior.T <- density( abc.postSC$adj.values[,14] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,20),
      ylim=c(0,0.2), 
      ylab="probability density", 
      xlab="shape2M2 ",
      lwd=2, 
      main="shape2M2 ")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,15] )
posterior.T <- density( abc.postSC$adj.values[,15] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,200),
      ylim=c(0,0.2), 
      ylab="probability density", 
      xlab="propNtrlM1",
      lwd=2, 
      main="propNtrlM1 ")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,16] )
posterior.T <- density( abc.postSC$adj.values[,16] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,6000),
      ylim=c(0,0.2), 
      ylab="probability density", 
      xlab="propNtrlM2 ",
      lwd=2, 
      main="propNtrlM2 ")
lines(posterior.T, col="blue",lwd=2)
dev.off()
########################################################################
#
#						Homo M Hetero N
#
########################################################################
pdf(file="im.homom.hetn.tol0.1.pdf",width=12,height=8)
par(mfrow=c(4,4))
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

prior.T    <- density( priorfile[,5] )
posterior.T <- density( abc.postSC$adj.values[,5] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="M1",
      lwd=2, 
      main="M1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,10] )
posterior.T <- density( abc.postSC$adj.values[,10] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.10), 
      ylab="probability density", 
      xlab="M2",
      lwd=2, 
      main="M2")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,11] )
posterior.T <- density( abc.postSC$adj.values[,11] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.1), 
      ylab="probability density", 
      xlab="propNtrlNe1",
      lwd=2, 
      main="propNtrlNe1")
lines(posterior.T, col="blue",lwd=2)


prior.T    <- density( priorfile[,8] )
posterior.T <- density( abc.postSC$adj.values[,8] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,6000),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="propNtrlNe2",
      lwd=2, 
      main="propNtrlNe2")
lines(posterior.T, col="blue",lwd=2)

dev.off()
#
#		Homo N Hetero M
#
##################################################"
pdf(file="im.heterom.homon.tol0.1.pdf",width=12,height=8)
par(mfrow=c(4,4))
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

prior.T    <- density( priorfile[,5] )
posterior.T <- density( abc.postSC$adj.values[,5] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,30),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="M1",
      lwd=2, 
      main="M1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,6] )
posterior.T <- density( abc.postSC$adj.values[,6] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.10), 
      ylab="probability density", 
      xlab="M2",
      lwd=2, 
      main="M2")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,7] )
posterior.T <- density( abc.postSC$adj.values[,7] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,40),
      ylim=c(0,0.1), 
      ylab="probability density", 
      xlab="shape1M1",
      lwd=2, 
      main="shape1M1")
lines(posterior.T, col="blue",lwd=2)


prior.T    <- density( priorfile[,8] )
posterior.T <- density( abc.postSC$adj.values[,8] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,200),
      ylim=c(0,0.6), 
      ylab="probability density", 
      xlab="shape2M1",
      lwd=2, 
      main="shape2M1")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,9] )
posterior.T <- density( abc.postSC$adj.values[,9] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,200),
      ylim=c(0,0.20), 
      ylab="probability density", 
      xlab="shape1M2 ",
      lwd=2, 
      main="shape1M2 ")
lines(posterior.T, col="blue",lwd=2)

prior.T    <- density( priorfile[,10] )
posterior.T <- density( abc.postSC$adj.values[,10] , weights=abc.postSC$weights/sum(abc.postSC$weights) )
print_tau <- expression(tau)
plot( prior.T,
      xlim=c(0,200),
      ylim=c(0,0.2), 
      ylab="probability density", 
      xlab="shape2M2 ",
      lwd=2, 
      main="shape2M2 ")
lines(posterior.T, col="blue",lwd=2)
dev.off()
