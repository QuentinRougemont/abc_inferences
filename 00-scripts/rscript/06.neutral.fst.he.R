argv <- commandArgs(TRUE)        #pas utile pour l'instant
input <- argv[1]                 #pas utile pour l'instant #Nom dataframe
nsp1  <- as.numeric(argv[2])
nsp2  <- as.numeric(argv[3]) 
#miss_rate <- as.numeric(argv[4]) #pas utile pour l'instant #Seuil de missing

geno <- as.numeric(t(read.table(input)))

#function for heterozygotie
Het <- function(x)
{
    n <- length(x)
    f <- table(x)/n
    p <- sum(f^2)
    Het <- n * (1 - p) / (n - 1)
  return(Het)
}

x<-seq_along(geno)
ntot <- nsp1 + nsp2 
x1 <- split(geno, ceiling(x/ntot))

sample_size_pop1 <-nsp1 
sample_size_pop2 <- nsp2
sample_size_total <- as.matrix(sample_size_pop1+sample_size_pop2)

stat_data <- "neutral.stats.txt" #la matrice finale
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

x<-x1
for(i in 1:length(x))
    {
    h1 <- Het(as.factor(as.numeric(x[[i]][1:sample_size_pop1[1]]))) #Rq : on peut démontrer que équivalent à PiA dans le cas single SNPs comme ici
    h2 <- Het(as.factor(as.numeric(x[[i]][(sample_size_pop1[1]+1):sample_size_total[1]]))) #PiB
    htot <- Het(as.factor(as.numeric(x[[i]]))) #PiT
    gst <- (htot - mean(cbind(h1,h2)))/htot #équivaut donc à 1-Pi12/PiTot
    he_1 <- 1-sum((table(as.factor(as.numeric(x[[i]][1:sample_size_pop1[1]])))/sum(table(as.factor(as.numeric(x[[i]][1:sample_size_pop1[1]]))),na.rm=T))^2)
    he_2 <- 1-sum((table(as.factor(as.numeric(x[[i]][(sample_size_pop1[1]+1):sample_size_total[1]])))/sum(table(as.factor(as.numeric(x[[i]][(sample_size_pop1[1]+1):sample_size_total[1]]))),na.rm=T))^2)
    he_tot <-1-sum((table(as.factor(as.numeric(x[[i]])))/sum(table(as.factor(as.numeric(x[[i]]))),na.rm=T))^2)
    fis <- 1-htot/he_tot
    write.table(cbind(h1,h2,htot,he_1,he_2,he_tot,fis,gst),stat_data,quote=F,col.names=F,row.names=F,append=T)
}

