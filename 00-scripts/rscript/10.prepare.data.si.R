##!/usr/bin/Rscript

argv   <- commandArgs(TRUE)        #pas utile pour l'instant 
input  <- argv[1]                 #pas utile pour l'instant #Nom dataframe
thresh <- as.numeric(argv[2])

geno  <- read.table("geno.tmp")

header <- read.table(input)[1,] #Argv2
header <- read.table("ABC.tmp")[1,] #Argv2

header <- header[-1]
p1 <-unique(t(header))[1,]
p2 <-unique(t(header))[2,]
sp1 <- 2*length(which(header==p1)) #sp.1
sp2 <- 2*length(which(header==p2)) #sp.2

geno<-as.matrix(geno)

for(i in 1:ncol(geno)){
geno[which(geno[,i]=="-"),i]<-NA
}

if(!is.matrix(geno) | !mode(geno)=="character") { stop("geno must be of 'matrix' class and 'character' mode.\n") }

m <- data.frame( N=rowSums(!is.na(geno)) ) 

allele <- apply(substr(geno,1,1),1,unique) 
if( is.matrix(allele) ) { allele <- lapply(apply(allele,1,as.list),as.character) }
allele <- lapply(allele,sort)
m$nballele <- sapply(allele,length)
if( any(m$nballele>2) ) { stop("geno contains more than two alleles.\n") }

ref = NULL

m$ancestral <- NA
inds <- which(m$nballele>0)
m$ancestral[inds] <- sapply(allele[inds],'[[',1)
m$derived <- NA
inds <- which(m$nballele>1)
m$derived[inds] <- sapply(allele[inds],'[[',2)

mono=(which(m$nballele==1))

geno2  <- read.table("geno.tmp")
geno.2 <- geno2[-c(mono),] #virer dans geno les markeurs monomorphes

ancestral.state=as.matrix(geno.2[,1])
geno.3<-matrix(ncol=ncol(geno.2),nrow=nrow(geno.2))

for(i in 1:nrow(geno.2)){
      for(j in 1:ncol(geno.2)){
      if(geno.2[i,j]==ancestral.state[i,])
        {
        geno.3[i,j] <- 0
        }
        else
        	if(geno.2[i,j]=="NA")
        	{
        	geno.3[i,j]<-NA
        	}
        	else
        	geno.3[i,j]<- 1
      }

}

for (i in 1:ncol(geno.3)){      #a remplacer par un sapply ou équivalent
geno.3[which(geno.3[,i]=="-"),i]<-NA
  }

geno.5<-as.data.frame(t(geno.3)) 
tmp <- matrix(ncol=ncol(geno.5), nrow=1)
for(i in 1:ncol(geno.5)){

    tmp[,i]=sum(is.na(geno.5[,i]))

}

m<-m[-c(mono),]
#on veux virer ceux avec trop de missing data:
miss_rate <- 0.2 #A passer en argument dans le script
rate <-m[,1] * miss_rate
miss <- which(tmp>=rate) #attention le seuil de 8 sera à adapter #trouver critere, a mettre en argument à passer au script R dans le bash

if (length(miss)>0)
{
geno.6 <- geno.5[,-miss]
} else {
geno.6 <- geno.5
}

geno.30 <- read.table("ABC", skip=2)
geno.30 <- geno.30[-mono,]

if (length(miss)>0)
{
geno.40 <- geno.30[,-miss]
} else {
geno.40 <- geno.30
}

write.table(geno.40,"matrix.filter.snps.txt",quote=F,row.names=F,col.names=F)

fst   <- read.table("../03.fst.he.data/empirical.stats2.txt",T)[,8]

q   <- quantile(fst,thresh/100)

fst <- which(fst>q)

zoro  <- geno.40[fst,]

eder2 <- read.table("ABC.tmp")

write.table(rbind(eder2,as.matrix(zoro)),file=paste(p1,p2,".",thresh,".txt",sep=""),quote=F,row.names=F,col.names=F, sep="\t")
