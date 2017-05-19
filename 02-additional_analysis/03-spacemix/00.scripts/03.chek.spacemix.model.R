#!/usr/local/bioinfo/src/R/R-3.2.2/bin/Rscript

#Author= "Quentin Rougemont"
#purpose = "check validity spacemix results" (see https://github.com/gbradburd/SpaceMix)
#how to run = "03.chek.spacemix.model.R arg1 arg2 arg3 
#last uptade = "06.09.2016"

argv <- commandArgs(TRUE)

spmix.data<-argv[1] 
mcn.freq.list<-argv[2]
mcmc.output<-argv[3]

#SpaceMix Control and Plot 
library(SpaceMix)

load(spmix.data)
load(mcn.freq.list)
load(mcmc.output)

#check the data
#str(MCN.frequencies.list)
#str(spacemix.data)
#ls()

#Plot MCMC trace 
pdf(file="1.spacemix.log.pdf")
plot(Prob,xlab="MCMC iterations",ylab="value",
    main="Posterior probability trace plot",type='l')
    
dev.off()
#### Plot Trace of Nugget parameters
pdf(file="2.spacemix.nugget.pdf")    
matplot(t(nugget),type='l',
            xlab="MCMC iterations",ylab="Parameter value",
            main="Trace plot of nugget parameters")    
dev.off()

#Joint Marginal Plots
pdf(file="3.joint.marginal.pdf")
plot(a0,a1,xlab="a0",ylab="a1",
    main="Joint marginal of a0 and a1",pch=20,
    col=adjustcolor(rainbow(1000,start=4/6,end=6/6),0.3))
legend(x="bottomright",pch=19,cex=0.8,
        col=rainbow(1000,start=4/6,end=6/6)[c(1,500,1000)],
        legend=c("Sampled MCMC iteration 1",
                 "Sampled MCMC iteration 500",
                 "Sampled MCMC iteration 1000"))
dev.off()

#Acceptance Rate
pdf(file="4.accept.rate.pdf")
plot(accept_rates$a0_accept_rate,
        xlab="MCMC iterations",ylab="Acceptance rate",
        main="Acceptance rate of a0",type='l',
        ylim=c(0,1))
    abline(h=0.44,col="gray",lty=2)
dev.off()

#Acceptatnce Rates of the nugget
pdf(file="5.accept.rate.nugget.pdf")
matplot(t(accept_rates$nugget_accept_rate),
            xlab="MCMC iterations",ylab="Acceptance rate",
            main="Acceptance rates of nuggets",type='l',
            ylim=c(0,0.7))
    abline(h=0.44,col="gray",lty=2)
dev.off()

######################Measurement of Model Adequacy ##############################
#calculate the sample covariance from the mean centered and normalized sample allele frequencies.
 sample.covariance <- cov(t(MCN.frequencies.list$mean.centered.normalized.sample.frequencies),
                                use="pairwise.complete.obs")
                                
k <- nrow(MCN.frequencies.list$mean.centered.normalized.sample.frequencies)
MC.matrix <- diag(k) - matrix(1/last.params$inv.mean.sample.sizes / 
                                    (sum(1/last.params$inv.mean.sample.sizes)),
                                        nrow=k,ncol=k,byrow=TRUE)
                                
MC.parametric.covariance <- (MC.matrix) %*% last.params$admixed.covariance %*%  t(MC.matrix) #

index.matrix <- upper.tri(sample.covariance,diag=TRUE)

pdf(file="6.model.adequacy.pdf")
plot(sample.covariance[index.matrix], 
    MC.parametric.covariance[index.matrix],
    col=adjustcolor("black",0.3),pch=20,
    xlab="Sample covariance",
    ylab="Parametric covariance",
    main="Model adequacy:\n matrix comparison")
    abline(0,1,col="red")                              
dev.off()

pdf(file="7.ibd.trend.pdf")
plot(last.params$D[1:k,1:k][index.matrix], 
        sample.covariance[index.matrix],
        pch=19,col="black",
        xlab="geogenetic distance",
        ylab="covariance",
        main="Model adequacy:\n IBD patterns")
        points(last.params$D[1:k,1:k][index.matrix], 
                MC.parametric.covariance[index.matrix],col="red",pch=20)
        legend(x="topright",pch=19,col=c(1,2),
                legend=c("observed","model estimate"))
                
dev.off()
