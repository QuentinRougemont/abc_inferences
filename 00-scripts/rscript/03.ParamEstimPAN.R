source('00-scripts/rscript/cv4abc.R')
if("abc" %in% rownames(installed.packages()) == FALSE)
    {install.packages("abc") }

library(abc)

argv<-commandArgs(TRUE)
model  <-argv[1] #identification of the model either 1,2,3,4
model <- 1
#load data
target=as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))
target=target[-c(1:3,12,13,18:25)]
obs=target[1:19]

nlinesFul=as.numeric(strsplit(system
    ("wc -l 00.data/pan.homom.homon.ABC.stat.txt ", intern=T), " ")[[1]][1])
if(model==1){
M1=matrix(scan("00.data/pan.homom.homon.ABC.stat.txt"), 
    byrow=T, nrow=nlinesFul) 
}else if (model == 2 ){
M1=matrix(scan("00.data/pan.homom.heteron.ABC.stat.txt"), 
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
priorfile=matrix(scan("00.data/pan.homom.homon.priorfile.txt"), 
    byrow=T, nrow=nlinesFul) 
}else if (model == 2 ){
priorfile=matrix(scan("00.data/pan.heterom.heteron.priorfile.txt"), 
    byrow=T, nrow=nlinesFul) 
} 

if(model==1){
colnames(priorfile)=c("N1")
}else if (model == 2 ){
colnames(priorfile)=c("N1",
                    "shape1Ne", "shape2Ne", "propNtrlNe1","propNtrlNe2")
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
sink("postadjval.pan.homo.tol0.1.txt")
print(abc.post$adj.values)
sink()
}else if (model == 2 ){
sink("postadjval.pan.heteron.tol0.1.txt")
print(abc.post$adj.values)
sink()
} 

if(model==1){
sink("postadjval.pan.homo.tol0.1.txt")
print(summary(abc.post))
sink()
}else if (model == 2 ){
sink("postadjval.pan.heteron.tol0.1.txt")
print(summary(abc.post))
sink()
} 
pdf(file="pan.homo.homo.tol0.1.pdf",width=12,height=8)
par(mfrow=c(1,2))
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

dev.off()
