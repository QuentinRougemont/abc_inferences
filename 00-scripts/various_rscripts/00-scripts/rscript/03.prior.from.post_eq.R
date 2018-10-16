library(abc)

target=as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))
target=target[-c(1:3,12,13,18:25)]
obs=target[1:19]

nlinesFul <- as.numeric(strsplit(system("wc -l 00.data/eq.homom.homon.ABC.stat.txt ", intern=T), " ")[[1]][1])
M1 <- matrix(scan("00.data/eq.homom.homon.ABC.stat.txt"), byrow=T, nrow=nlinesFul)
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

abc.post <- abc(target  = obs, 
                param = priorfile, 
                sumstat = M1b, 
                tol = 1000/1e6, 
                transf=c(rep("logit",ncol(priorfile)) ),
         logit.bounds = log.bound,
         hcorr=T, method  = "neuralnet", numnet=50, sizenet=15)

simpost <- 5000
newsamp <- sample(1:(dim(abc.post$adj)[1]),size=simpost,replace=T,prob=abc.post$weights) 
newsamp <- abc.post$adj[newsamp,]
write.table(newsamp,"results/eq1.prior.from.post", quote=F, row.names=F, col.names=F)

