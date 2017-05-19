#!/usr/bin/Rscript

# Global variables
argv    <- commandArgs(TRUE)
lowbnd  <- as.numeric(argv[1])
numlines <- 50 #100
uperbnd <- lowbnd + numlines - 1
target  <- argv[2] 

# Source ABC function
source('01.script/cv4abc.R')

# Read model data
nlinesFul=as.numeric(strsplit(system("wc -l  00.data/im.homom.homon.ABC.stat.txt      ",  intern=T), " ")[[1]][1])

IM.1=matrix(scan("00.data/im.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul)
IM.2=matrix(scan("00.data/im.heterom.homon.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul)
IM.3=matrix(scan("00.data/im.homom.heteron.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul)
IM.4=matrix(scan("00.data/im.heterom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul)
AM.1=matrix(scan("00.data/am.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul)
AM.2=matrix(scan("00.data/am.heterom.homon.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul)
AM.3=matrix(scan("00.data/am.homom.heteron.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul)
AM.4=matrix(scan("00.data/am.heterom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul)
SC.1=matrix(scan("00.data/sc.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul)
SC.2=matrix(scan("00.data/sc.heterom.homon.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul)
SC.3=matrix(scan("00.data/sc.homom.heteron.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul)
SC.4=matrix(scan("00.data/sc.heterom.heteron.ABC.stat.txt"), byrow=T, nrow=nlinesFul)
SI.1=matrix(scan("00.data/si.homom.homon.ABC.stat.txt"    ), byrow=T, nrow=nlinesFul)
SI.3=matrix(scan("00.data/si.homom.heteron.ABC.stat.txt"  ), byrow=T, nrow=nlinesFul) 

# Replace missing data
for(i in 1:ncol(IM.1)){
  IM.1[which(IM.1[,i]=="NaN"),i]=mean(IM.1[,i], na.rm=T)
  IM.2[which(IM.2[,i]=="NaN"),i]=mean(IM.2[,i], na.rm=T)
  IM.3[which(IM.3[,i]=="NaN"),i]=mean(IM.3[,i], na.rm=T)
  IM.4[which(IM.4[,i]=="NaN"),i]=mean(IM.4[,i], na.rm=T)
  AM.1[which(AM.1[,i]=="NaN"),i]=mean(AM.1[,i], na.rm=T)
  AM.2[which(AM.2[,i]=="NaN"),i]=mean(AM.2[,i], na.rm=T)
  AM.3[which(AM.3[,i]=="NaN"),i]=mean(AM.3[,i], na.rm=T)
  AM.4[which(AM.4[,i]=="NaN"),i]=mean(AM.4[,i], na.rm=T)
  SC.1[which(SC.1[,i]=="NaN"),i]=mean(SC.1[,i], na.rm=T)
  SC.2[which(SC.2[,i]=="NaN"),i]=mean(SC.2[,i], na.rm=T)
  SC.3[which(SC.3[,i]=="NaN"),i]=mean(SC.3[,i], na.rm=T)
  SC.4[which(SC.4[,i]=="NaN"),i]=mean(SC.4[,i], na.rm=T)
  SI.1[which(SI.1[,i]=="NaN"),i]=mean(SI.1[,i], na.rm=T)
  SI.3[which(SI.3[,i]=="NaN"),i]=mean(SI.3[,i], na.rm=T)
}

# Extract working data (lower and upper bound and wanted columnt)
IM.1b=IM.1[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
IM.2b=IM.2[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
IM.3b=IM.3[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
IM.4b=IM.4[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
SC.1b=SC.1[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
SC.2b=SC.2[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
SC.3b=SC.3[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
SC.4b=SC.4[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
AM.1b=AM.1[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
AM.2b=AM.2[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
AM.3b=AM.3[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
AM.4b=AM.4[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
SI.1b=SI.1[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]
SI.3b=SI.3[-c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]

nlineful2 <- nlinesFul-numlines

# Comparaison All
x <- as.factor(rep(c(1:14),each=nlineful2))
z = rbind(SI.1b,SI.3b,IM.1b,IM.2b,IM.3b,IM.4b,AM.1b,AM.2b,AM.3b,AM.4b,SC.1b,SC.2b,SC.3b,SC.4b)

targetdata = eval(parse(text = target))
obs = targetdata[c(lowbnd:uperbnd),-c(1:3,12,13,18:25)]

res = model_selection_abc_nnet(target=obs, x=x, sumstat=z, tol=3500/(length(x)), noweight=F, rejmethod=F, nb.nnet=50, size.nnet=15, output=paste("03.results/",target,"/robustess_",target, "_",lowbnd,sep=""))
