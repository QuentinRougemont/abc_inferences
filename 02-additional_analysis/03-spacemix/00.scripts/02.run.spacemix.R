#!/usr/bin/env Rscript

#Author= "Quentin Rougemont"
#purpose = "analyse spacemix results" (see https://github.com/gbradburd/SpaceMix)
#how to run = "02.run.spacemix.R sample.covariance.matrix sample.size x_y.coordinates fast.model.opts long.model.opts nloci 
#last uptade = "06.09.2016"

#load library
library(SpaceMix)


argv <- commandArgs(TRUE)

sample.cov<-argv[1] #pool catalog
samp.size<-argv[2]
coord<-argv[3]

sample.cov<-read.table(sample.cov,skip=1)
sample.cov<-sample.cov[,-c(1)]
sample.cov<-as.matrix(sample.cov)
sample.siz<-read.table(samp.size)[,2]
geo.coord <- as.matrix(read.table(coord)[,c(2:3)])

# Data option: allele counts and sample sizes
# Fast Model option: estimating geogenetic locations
# Long Model option: estimating geogenetic locations and 
#                    admixture source locations
# Spatial priors: default variance,
#                   observed geographic sampling locations
spm=run.spacemix.analysis(n.fast.reps = 5,
                        fast.MCMC.ngen = 1e6,
                        fast.model.option = "source_and_target",
                        long.model.option = "source_and_target",
                        data.type = "sample.covariance",
                        sample.frequencies=NULL,
                        mean.sample.sizes=sample.siz,
                        counts =NULL ,
                        sample.sizes = NULL,
                        sample.covariance=sample.cov,
                        target.spatial.prior.scale=NULL,
                        source.spatial.prior.scale=NULL,
                        spatial.prior.X.coordinates = geo.coord[,1],
                        spatial.prior.Y.coordinates = geo.coord[,2],
                        round.earth = TRUE,
                        long.run.initial.parameters=NULL,
                        k = nrow(sample.cov),
                        loci = 4656,
                        ngen = 1e6,
                        printfreq = 1e2,
                        samplefreq = 1e3,
                        mixing.diagn.freq = 50,
                        savefreq = 1e5,
                        directory=NULL,
                        prefix = "salmon")
