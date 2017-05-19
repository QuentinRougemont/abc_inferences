#!/usr/bin/env Rscript
#QR - 12-09-16
#Input: #A genotypic matrix containing individuals in rows and loci in colomns (AG,GG,GC,TC,...), NA allowed
# Output:
#1. allele frequency normalized covariance table (see Bradburg et al. 2016)
#2. mean sample sizes
#3. allele count for each inds and loci 
#4. sample size  for each inds and loci

argv <- commandArgs(TRUE) 

dat<-argv[1] 
 
dat0=read.table(dat,h=F)

dat1<-dat0[,-c(1:2)]
genot<-t(as.matrix(dat1)) #set class to matrix

m <- data.frame( N=rowSums(!is.na(genot)) ) #numbers of individuals successfully genottyped at each m

alleles <- apply(cbind(substr(genot,1,1),substr(genot,2,2)),1,unique) #grep the different alleles at each loci
if( is.matrix(alleles) ) { alleles <- lapply(apply(alleles,1,as.list),as.character) } #transform the matrice as a class list and as characters
alleles <- lapply(alleles,sort) #sort alleles, alphabetical order; ignore NA
m$numAlleles = sapply(alleles,length) #
if( any(m$numAlleles>2) ) { stop("genot contains more than two alleles.\n") } #check that no m have more than 2 different alleles

m$A1 <- NA
inds <- which(m$numAlleles>0) #ceux sans NA
m$A1[inds] <- sapply(alleles[inds],'[[',1)
m$A2 <- NA
inds <- which(m$numAlleles>1)
m$A2[inds] <- sapply(alleles[inds],'[[',2)

ref <- m$A1
alt <- NA
inds <- which(ref==m$A1); alt[inds] <- m$A2[inds]
inds <- which(ref==m$A2); alt[inds] <- m$A1[inds]

if( any(ref!=m$A1 & ref!=m$A2) ) { warning("ref allele not present in genot for some m. Conversions for these m cannot be performed and will be coerced to NA.\n") }

m$G2 = paste(ref,ref,sep="")	#2 copies of ref
m$G1.1 = paste(ref,alt,sep="")	#1 copy of ref, ref allele coded first
m$G1.2 = paste(alt,ref,sep="")	#1 copy of ref, reversed coding
m$G0 = paste(alt,alt,sep="")	#0 copy of ref
m$G2[is.na(ref)] <- NA #imputing Missing
m$G1.1[is.na(alt)] <- NA
m$G1.2[is.na(alt)] <- NA
m$G0[is.na(alt)] <- NA

genot.mat <- matrix( NA, ncol=ncol(genot), nrow=nrow(genot), dimnames=dimnames(genot) ) #create final genotypic matrix
genot.mat[genot==m$G2] <- 2 #impute values
genot.mat[genot==m$G1.1 | genot==m$G1.2] <- 1
genot.mat[genot==m$G0] <- 0

#Get Sample Size
SampleSize <- matrix( NA, ncol=ncol(genot), nrow=nrow(genot), dimnames=dimnames(genot) )
SampleSize[genot.mat==2] <- 2
SampleSize[genot.mat==1] <- 2
SampleSize[genot.mat==0] <- 2
SampleSize[is.na(SampleSize)] <-0
SampSizful=t(SampleSize)
write.table(cbind((dat0[,1:2]),SampSizful), "samp.siz.full.data",quote=F,row.names=F,col.names=F) #write table of sample size

samp.siz.1<-cbind(dat0[,1:2],SampSizful)

genot.mat[is.na(genot.mat)] <-0

write.table(cbind((dat0[,c(1:2)]),t(genot.mat)),'allel.count.full.data',quote=F,row.names=F) #write table of allelic count

#system(" cat loci.list tmp >> allele.count.salm.non.admix ",wait=F) #write table of allele count
genot.1<-cbind(dat0[,c(1:2)],t(genot.mat))

#Create population sample size matrix
sample.sizes<- mapply(rowsum, as.data.frame(samp.siz.1[,c(3:ncol(samp.siz.1))]), as.data.frame(samp.siz.1[,c(1)]))
dimnames(sample.sizes)<- list(levels(factor(samp.siz.1[,c(1)]))) #, 1:nrow(samp.siz.1[,c(1)]))

#Create the vector of mean sample size
apply(sample.sizes,1,mean, na.rm=T) -> z
write.table(z,"mean.sample.sizes",quote=F,col.names=F)

#Create pop. allel freq.list 
tmp1<- mapply(rowsum, as.data.frame(genot.1[,c(3:ncol(genot.1))]), as.data.frame(genot.1[,c(1)]))
dimnames(tmp1)<- list(levels(factor(genot.1[,c(1)]))) #, 1:nrow(samp.siz.1[,c(1)]))
frequencies<-tmp1/sample.sizes

#Now calculate the weigthed mean allele frequencies at each locus (f* in Bradburg et al. 2016)
#for learning only
calculate.mean.sample.freq <- function(frequencies,sample.sizes){
frequencies[is.na(frequencies)] <- 0
weighted.sample.frequencies <- colSums(frequencies * sample.sizes) /
colSums(sample.sizes)
return(weighted.sample.frequencies)
}
tmp3<-calculate.mean.sample.freq(frequencies,sample.sizes)

#Check for monomorphic loci :
if( any(tmp3==1) ) { warning("Warning, non polymoprphic m still present in the dataset, covariance won't be computed \n please clean your input file") }

#to_remove=(which(tmp3==1))
#dat3=dat1[,-c(to_remove)]

#Now normalize the variance accross locus: (f*/sqrt(f(1-f*)
# First, define a normalizing function:
normalize.allele.freqs <- function(frequencies,sample.sizes){
mean.freqs <- calculate.mean.sample.freq(frequencies,sample.sizes)
mean.freqs.mat <- matrix(mean.freqs,
nrow=nrow(frequencies),
ncol=ncol(frequencies),
byrow=TRUE)
normalized.freqs <- frequencies/sqrt(mean.freqs.mat * (1-mean.freqs.mat))
return(normalized.freqs)
}

# Application to data:
normalized.allele.freqs <- normalize.allele.freqs(frequencies,sample.sizes)
nrow(normalized.allele.freqs) #check number of rows
#Now compute sample covariance of the normalized allele frequencies
sample.cov <- cov(t(normalized.allele.freqs),use="pairwise.complete.obs")
write.table(sample.cov,"sample.cov", quote=F)
