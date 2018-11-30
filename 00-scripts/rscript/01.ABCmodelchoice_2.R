source('/home/quentin/script/rscript/cv4abc.R')
if("data.table" %in% rownames(installed.packages()) == FALSE) 
{install.packages("qqman", repos="https://cloud.r-project.org") 
    print("installing packages data.table..." ) }

library(data.table)


target=as.numeric(read.table("OBS.ABC.stat.txt",skip=2,h=F))

SI1=fread("zcat 00.data/si.homom.homon.ABC.stat.txt.gz"    )
SI3=fread("zcat 00.data/si.homom.heteron.ABC.stat.txt.gz"  )

f <- function(x){
    m <- mean(x, na.rm = TRUE)
    x[is.na(x)] <- m
    x
}

M_PA1 <-fread("zcat 00.data/pan.homom.homon.ABC.stat.txt.gz"    ) #les 1e6 simuls sont mergés ensemble au préalable avec du bash et sed.
M_IM1 <-fread("zcat 00.data/im.homom.homon.ABC.stat.txt.gz"    ) #les 1e6 simuls sont mergés ensemble au préalable avec du bash et sed.
M_BO1 <-fread("zcat 00.data/bot.homom.homon.ABC.stat.txt.gz"    ) #les 1e6 sSIuls sont mergés ensemble au préalable avec du bash et sed.
M_PA1 <- apply(M_PA1, 2, f)
M_IM1 <- apply(M_IM1, 2, f)
M_BO1 <- apply(M_BO1, 2, f)

M_PA1b=M_PA1[,-c(1:3,12,13,18:25)]
M_IM1b=M_IM1[,-c(1:3,12,13,18:25)]
M_BO1b=M_BO1[,-c(1:3,12,13,18:25)]

target=target[-c(1:3,12,13,18:25)]
obs=matrix(rep(target[1:19],10), byrow=T, nrow=10)

#Comparaison All:
x=as.factor(c(rep("A",nlinesFul),rep("B",nlinesFul),rep("C",nlinesFul))),

z=rbind(M_BO1b,M_IM1b,M_PA1b)
res5=model_selection_abc_nnet(target=obs, 
                              x=x, 
                              sumstat=z,
                              tol=1000/(length(x)), 
                              noweight=F, 
                              rejmethod=F, 
                              nb.nnet=50, 
                              size.nnet=15, 
                              output="ALL.MODELS.29.06.16")

mean_all=apply(res5,2,mean)
sd_all=apply(res5,2,sd)
write.table(mean_all,"mean_all",quote=F,row.names=c("BO1","IM1","PA1b"),col.names=F)
