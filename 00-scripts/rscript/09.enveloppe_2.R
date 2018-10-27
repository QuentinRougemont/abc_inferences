#!/usr/bin/Rscript

argv <- commandArgs(TRUE)

neut <- argv[1]
emp  <- argv[2]

a <- read.table(neut,T) 
b <- read.table(emp,T)
a <- read.table("neutral.stats2.txt",T)
b <- read.table("statistics_with_WC_Fst", T)

pdf(file="fst_he_simulated_data2.pdf")                
plot(a[,6],a[,9],xlab="Expected Heterozygosity", ylab = "simulated Fst", pch="*")
dev.off()

x=0 # initiate x 
y=0 # initiate y
probs=c(0.005,0.025,0.25,0.5,0.75,0.975,0.995,0.99999)
#z = matrix(nrow=length(seq(0.025, 0.5,0.025)), 
#            ncol=length(probs)+1) #NULL
z = NULL
while(x < 0.5) { #loop while incrémentée de 0.025, pour un pas de 0.025 en 0.025 
          x <- x+0.025 # borne max
  print(x)
  print(y)
  v=NULL
        for (i in 1:nrow(a)){ # pour toutes les valeurs 
              if (a[i,6] < x ){
              if (a[i,6] > y){
              v = rbind(v, a[i,c(6,9)]) # value with name
               }
             } 
           }
   # imprimer les résultats des quantiles pour chaque milieu de classe (x-0.025, si l'incrémentation est de 0.025)
  z = rbind(z, c(x,quantile(v[2,], probs=probs,na.rm=T)))
  y <- y+0.025 # borne min
}
z <- as.data.frame(z)

pdf(file="neutral_enveloppe_with_outliers.pdf",12,10)
plot(b[,4]~b[,11],col='black',pch="*",
     xlab="Expected Heterozygosity", 
     ylab= "Simulated Fst",
     ylim=c(0:1))
points(z[,8]~z[,1],type="l", 
       lty = 2, col = "blue", 
       ylim=c(0:1), lwd="2")
points(z[,9]~z[,1],type='l', col="red")
dev.off()

#zz<-subset(b, b[,6]>0.025 & b[,6]<0.050 & b[,8]>0.01979980 & b[,8]<0.020)

outli<-"outlier.count"
write.table(cbind(
    "SNP_ID",
    "POS",
    "CHROM",
    "WCFST",
    "CHROM",
	"H1",
	"H2",
	"H_Tot",
	"He1",
	"He2",
	"HeTot",
	"Fis",
	"Fst"),
	file=outli,sep=" ",
	row.names=F,col.names=F,quote=F)

for(i in 1:nrow(z))
	{
	tmp <- subset(b, b[,11]>z[i,1] & b[,11]<z[i+1,1] & b[,4]>z[i,9] & b[,4]<z[i+1,9])
	write.table(tmp,outli,append=T,col.names=F,row.names=F, quote=F)
}

tmp2 <- subset(b, b[,11] >z[19,1] & b[,4]>z[20,9])
write.table(tmp2,outli, append =T,col.names=F,
	row.names=F,quote=F)

write.table(z,"quantile.distrib",col.names=c(
	"he",
	"q0.005",
	"q0.025",
	"q0.25",
	"q0.5",
	"q0.75",
	"q0.975",
	"q0.995",
	"q0.99999"),
	row.names=F,quote=F)
