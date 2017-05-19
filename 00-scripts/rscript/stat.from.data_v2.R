##!/usr/bin/Rscript
##!/usr/local/bioinfo/src/R/R-3.2.2/bin/Rscript #pour genotoul
#fst and other stats .from.data.R
#here missing data in genotypic table should be coded as '-' #coulde make it more generable later

#z<-paste("#!/bin/bash","\n", "sed -e '1,2d' EmtCRO.ABC.tmp |cut -f2-20000 | sed -r 's/([^ ])/\1 /g' | sed -e 's/\t/ /g'  >> pouet ") #Il y a un pb au niveau de la regexp (2nd sed)
#write.table(z,"reshape.sh",quote=F,col.names=F,row.names=F)

#cmd<-"bash ../bin/reshape.sh" 
#system(cmd, wait = T)

argv <- commandArgs(TRUE)        #pas utile pour l'instant 
input <- argv[1]                 #pas utile pour l'instant #Nom dataframe
#header <- argv[2]                  #pas utile pour l'instant #Non genotypic matrice (header only are read)
#miss_rate <- as.numeric(argv[4]) #pas utile pour l'instant #Seuil de missing

#geno <- read.table(argv[1])
geno <- read.table("pouet")

#function for heterozygotie
Het <- function(x)
{
    n <- length(x)
    f <- table(x)/n
    p <- sum(f^2)
    Het <- n * (1 - p) / (n - 1)
  return(Het)
}

header <- read.table(input)[1,] #Argv2
header <- header[-1]
p1 <-unique(t(header))[1,]
p2 <-unique(t(header))[2,]
sp1 <- 2*length(which(header==p1)) #sp.1
sp2 <- 2*length(which(header==p2)) #sp.2

allele <- apply(geno,1,unique)
allele.sp1 <- apply(geno[,1:sp1],1,unique)
allele.sp2 <- apply(geno[,(sp1+1):(sp1+sp2)],1,unique)

m <- data.frame( N=rowSums(!is.na(as.matrix(geno)) ) ) #numbers of individuals successfully genottyped at each marqueurs
m$numallele    <- sapply(allele,length) #number of alleles
m$numallelesp1 <- sapply(allele.sp1,length) #number of alleles
m$numallelesp2 <- sapply(allele.sp2,length) #number of alleles

#length(which(m$numallele==1))
mono<-(which(m$numallele==1)) #list of monomorphic markers

geno.2<-geno[-c(mono),] #virer dans geno les markeurs monomorphes
#on défini l'état ancestral à partir de ce qu'on obtiens dans la premiere colonne. 
ancestral.state=as.matrix(geno.2[,1])
geno.3<-matrix(ncol=ncol(geno.2),nrow=nrow(geno.2))

for(i in 1:nrow(geno.2)){
      for(j in 1:ncol(geno.2)){
      if(geno.2[i,j]==ancestral.state[i,])
        {
        geno.3[i,j] <- 0
        }
        else
        	if(geno.2[i,j]=="-")
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

#on veux virer ceux avec trop de missing data:
miss_rate <- 0.1 #A passer en argument dans le script
rate <-m[1,1] * miss_rate
miss <- which(tmp>=rate) #attention le seuil de 8 sera à adapter #trouver critere, a mettre en argument à passer au script R dans le bash
geno.6 <- geno.5[,-miss] #On a 4537 SNPs #Contre 4534 dans spinput
geno.6<-as.list(geno.6)
x <- geno.6 
#x <- x[-(which(sapply(x,is.na)))] #,arr.ind=TRUE))] #4531 on perd 6 SNPs!!!!! 
x <- na.omit(geno.6)

#preparing bpfile:
for (i in 1:ncol(geno)){      #a remplacer par un sapply ou équivalent
geno[which(geno[,i]=="-"),i]<-NA
  }

m$Nsp1 <-rowSums(  !is.na(as.matrix(geno[,1:sp1]) )  )
m$Nsp2 <--rowSums(  !is.na(as.matrix(geno[,(sp1+1):(sp1+sp2)]) )  ) #Ici ça devrait être positif !!!!
bp.tmp <- t(cbind(m$Nsp1,-m$Nsp2)) #du coup je reverse en mettant un "-" mais ça devrait pas !!!!
bp <- bp.tmp[,-c(mono)]
bp <- bp[,-c(miss)]


#il n"y a plus qu'a faire la boucle ou équivalent
sample_size_pop1 <- as.matrix(bp[1,])
sample_size_pop2 <- as.matrix(bp[2,])
sample_size_total <- as.matrix(sample_size_pop1+sample_size_pop2)

stat_data <- "empirical.stats.txt" #la matrice finale
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
    h1 <- Het(as.factor(as.numeric(x[[i]][1:sample_size_pop1[i]]))) #Rq : on peut démontrer que équivalent à PiA dans le cas single SNPs comme ici
    h2 <- Het(as.factor(as.numeric(x[[i]][(sample_size_pop1[i]+1):sample_size_total[i]]))) #PiB
    htot <- Het(as.factor(as.numeric(x[[i]]))) #PiT
    gst <- (htot - mean(cbind(h1,h2)))/htot #équivaut donc à 1-Pi12/PiTot
    he_1 <- 1-sum((table(as.factor(as.numeric(x[[i]][1:sample_size_pop1[i]])))/sum(table(as.factor(as.numeric(x[[i]][1:sample_size_pop1[i]]))),na.rm=T))^2)
    he_2 <- 1-sum((table(as.factor(as.numeric(x[[i]][(sample_size_pop1[i]+1):sample_size_total[i]])))/sum(table(as.factor(as.numeric(x[[i]][(sample_size_pop1[i]+1):sample_size_total[i]]))),na.rm=T))^2)
    he_tot <-1-sum((table(as.factor(as.numeric(x[[i]])))/sum(table(as.factor(as.numeric(x[[i]]))),na.rm=T))^2)
    fis <- 1-htot/he_tot
    write.table(cbind(h1,h2,htot,he_1,he_2,he_tot,fis,gst),stat_data,quote=F,col.names=F,row.names=F,append=T)
}
