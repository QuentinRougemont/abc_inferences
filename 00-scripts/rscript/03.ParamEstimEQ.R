#script to estimates parameters under EQ
#coalescent sims are stored in folder 00.data

source('00-scripts/rscript/cv4abc.R')
if("abc" %in% rownames(installed.packages()) == FALSE)
    {install.packages("abc") }
library(abc)

#load data
target=as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))
target=target[-c(1:3,12,13,18:25)]
obs=target[1:19]
nlinesFul=as.numeric(strsplit(system("wc -l 00.data/eq.homom.homon.ABC.stat.txt ", intern=T), " ")[[1]][1])

M1=matrix(scan("00.data/eq.homom.homon.ABC.stat.txt"), byrow=T, nrow=nlinesFul) 
#replace NA by means
clean <-function(model){
    means <-colMeans(model, na.rm=T)
    for (j in 1:ncol(model)){
         model[is.na(model[, j]), j] <- means[j]
    }
    #return(model)
}
f <- function(x){ 
m <- mean(x, na.rm = TRUE) 
x[is.na(x)] <- m 
x 
} 
M1 <- apply(M1, 2, f) 
M1b=M1[,-c(1:3,12,13,18:25)]

#load priorfile
priorfile=matrix(scan("00.data/eq.homom.homon.priorfile.txt"), 
    byrow=T, 
    nrow=nlinesFul) 
colnames(priorfile)=c("N1","N2","M1","M2")

log.bound=matrix(NA,nrow=ncol(priorfile), ncol=2)
for(i in 1:ncol(priorfile))
{
    log.bound[i,]=(range(priorfile[,i]))
}

abc.post <- abc(target  = obs, param = priorfile, 
    sumstat = M1b, 
    tol = 1000/1e6, 
    transf=c(rep("logit",ncol(priorfile) ) ), 
    logit.bounds = log.bound,
    hcorr=T, 
    method  = "neuralnet", 
    numnet=50, sizenet=15)

sink("postadjval.eq.homo.homo.tol0.1.txt")
print(abc.post$adj.values)
sink()

#Compute mean, mode, 95%HPD
sink("summary.post.eq.homo.homo.tol0.1.txt")
print(summary(abc.post))
sink()

pdf(file="eq.homo.homo.tol0.1.pdf",width=12,height=8)
par(mfrow=c(2,2))
prior.tetha1     <- density( (priorfile[,1] ) )
post.tetha1 <- density( abc.post$adj.values[,1] , 
    weights=abc.post$weights/sum(abc.post$weights) )
theta1 <- expression(theta[1])
#theta0 <- expression(theta[0])
plot( prior.tetha1,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta1,
      lwd=2, 
      main="theta1=4*N1*µ/4Nref*µ")  #or main= "Effective mutation rate ")
lines(post.tetha1, col="blue",lwd=2)
legend(
 legend = c("prior distribution","post"),
col = c("black","blue"),
lty = c(1,1), lwd = c(1,1),
x = "topright",
cex = 1,
bty ="n")
prior.tetha2     <- density( priorfile[,2] )
post.tetha2 <- density( abc.post$adj.values[,2] , 
    weights=abc.post$weights/sum(abc.post$weights) )
theta <- expression(theta[2])
plot( prior.tetha2,
      xlim=c(0,20),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab=theta,
      lwd=2, 
      main="theta2=4*N2*µ/4Nref*µ") #main="Effective mutation rate (lp)")
lines(post.tetha2, col="blue",lwd=2)
text(-1,8,"EQ",cex=1)
prior.M1     <- density( priorfile[,3] )
post.M1 <- density( abc.post$adj.values[,3] , 
    weights=abc.post$weights/sum(abc.post$weights) )
plot( prior.M1,
      xlim=c(0,40),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab="M1<-2",
      lwd=2, 
      main="Migration rate M1<-2")
lines(post.M1, col="blue",lwd=2)

prior.M2     <- density( priorfile[,4] )
post.M2 <- density( abc.post$adj.values[,4] , 
    weights=abc.post$weights/sum(abc.post$weights) )
plot( prior.M2,
      xlim=c(0,40),
      ylim=c(0,1), 
      ylab="probability density", 
      xlab="M2<-1",
      lwd=2, 
      main="Migration rate M2<-1")
lines(post.M2, col="blue",lwd=2)
dev.off()
