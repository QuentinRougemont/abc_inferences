#!/usr/bin/Rscript

#Script by QR - 15-09-16
#Miniscript to reshape a classical ped file into a genotypic matrix (usefull for several programms)
#Input files: 1) ped file 
#Output: genotypic matrix (AC, TG, AA, CC, etc. with inds in row, mk in cols)

argv <- commandArgs(TRUE) 

dat<-argv[1] #ped file 

dat2<-as.matrix(read.table(dat,h=F)) 
dat3<-dat2[,-c(1:6)] 

dat3[dat3 == "0"] <- NA

start <- seq(1, by = 2, length = ncol(dat3) / 2)
sdf <- sapply(start,function(i, dat3) paste(as.character(dat3[,i]),as.character(dat3[,i+1]), sep="") ,dat3 = dat3) 

write.table(cbind(dat2[,c(1:2)],as.data.frame(sdf)),"genotypic.matrix",quote=F,col.names=F,row.names=F)

system("sed -i 's/NANA/NA/g' genotypic.matrix*", wait=FALSE)
