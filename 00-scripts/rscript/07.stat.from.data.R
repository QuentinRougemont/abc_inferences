##!/usr/bin/Rscript

argv    <- commandArgs(TRUE)        #pas utile pour l'instant 
input_1 <- argv[1]                 #pas utile pour l'instant #Nom dataframe
input_2 <- argv[2]

input_2 <- "geno"
input_1 <- "header"
geno <- read.table(input_2)
header <- read.table(input_1)[1,] #Argv2

#function for heterozygotie #Copyright: E. pararadis. Pegas package
Het <- function(x)
{
    n <- length(x)
    f <- table(x)/n
    p <- sum(f^2)
    Het <- n * (1 - p) / (n - 1)
  return(Het)
}

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
m$nballele = sapply(allele,length)
if( any(m$nballele>2) ) { stop("geno contains more than two alleles.\n") }

ref=NULL

m$ancestral = NA
inds <- which(m$nballele>0)
m$ancestral[inds] <- sapply(allele[inds],'[[',1)
m$derived = NA
inds <- which(m$nballele>1)
m$derived[inds] <- sapply(allele[inds],'[[',2)

mono=(which(m$nballele==1))

geno2 <-read.table("geno.tmp")
#geno.2<-geno2[-c(mono),] #virer dans geno les markeurs monomorphes
geno.2 <- geno2
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

m <- m
#m<-m[-c(mono),]
#on veux virer ceux avec trop de missing data:
miss_rate <- 0.2 #A passer en argument dans le script
rate <-m[,1] * miss_rate
miss <- which(tmp>=rate) 

if (length(miss)>0)
{
geno.6 <- geno.5[,-miss]
} else {
geno.6 <- geno.5
}

geno.6<-as.list(geno.6)
x <- geno.6 
x <- na.omit(geno.6)

#preparing bpfile:
for (i in 1:ncol(geno)){      #a remplacer par un sapply ou équivalent
geno[which(geno[,i]=="-"),i]<-NA
  }
m$Nsp1 <-rowSums(  !is.na(as.matrix(geno[,1:sp1]) )  )
m$Nsp2 <--rowSums(  !is.na(as.matrix(geno[,(sp1+1):(sp1+sp2)]) )  ) # !!!!

#m$Nsp1 <-rowSums(  !is.na(as.matrix(geno[-mono,1:sp1]) )  )
#m$Nsp2 <--rowSums(  !is.na(as.matrix(geno[-mono,(sp1+1):(sp1+sp2)]) )  ) # !!!!
bp <- t(cbind(m$Nsp1,-m$Nsp2)) #du coup je reverse en mettant un "-" mais ça devrait pas !!!!
#bp <- bp.tmp[,-c(mono)]

if(length(miss)>0)
{
bp <- bp[,c(miss)]
} else {
bp <- bp 
}

#il n"y a plus qu'a faire la boucle ou équivalent
sample_size_pop1 <- as.matrix(bp[1,])
sample_size_pop2 <- as.matrix(bp[2,])
sample_size_total <- as.matrix(sample_size_pop1+sample_size_pop2)

stat_data <- "empirical.stats2.txt" #la matrice finale
write.table( cbind(
                  "H1",
                  "H2",
                  "H_total",
                  "He1",
                  "He2",
                  "He_total",
                  "fis",
                  "fst"
                 ),
              file=stat_data,sep=" ",
              quote=F,col.names=F,row.names=F,append=F)

for(i in 1:length(x))
    {
    h1     <- Het(as.factor(as.numeric(x[[i]][1:sample_size_pop1[i]]))) #Rq : on peut démontrer que équivalent à PiA dans le cas single SNPs comme ici
    h2     <- Het(as.factor(as.numeric(x[[i]][(sample_size_pop1[i]+1):sample_size_total[i]]))) #PiB
    htot   <- Het(as.factor(as.numeric(x[[i]]))) #PiT
    gst    <- (htot - mean(cbind(h1,h2)))/htot #équivaut donc à 1-Pi12/PiTot
    he_1   <- 1-sum((table(as.factor(as.numeric(x[[i]][1:sample_size_pop1[i]])))/sum(table(as.factor(as.numeric(x[[i]][1:sample_size_pop1[i]]))),na.rm=T))^2)
    he_2   <- 1-sum((table(as.factor(as.numeric(x[[i]][(sample_size_pop1[i]+1):sample_size_total[i]])))/sum(table(as.factor(as.numeric(x[[i]][(sample_size_pop1[i]+1):sample_size_total[i]]))),na.rm=T))^2)
    he_tot <-1-sum((table(as.factor(as.numeric(x[[i]])))/sum(table(as.factor(as.numeric(x[[i]]))),na.rm=T))^2)
    fis    <- 1-(htot/he_tot)
    write.table(cbind(h1,h2,htot,he_1,he_2,he_tot,fis,gst),stat_data,quote=F,col.names=F,row.names=F,append=T)
}
