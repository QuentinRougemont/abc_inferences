model_selection_abc_nnet <- function(target,x,sumstat,tol,gwt,noweight=F,rejmethod=F,nb.nnet=10,size.nnet=5, output="output")
{
require(nnet)
normalise <- function(x,y){

if(mad(y) == 0)
return (x)
else
return (x/mad(y))
}

###weight decay parameter
repet<-floor(nb.nnet/4)+1
the_decay<-rep(c(10^(-4),10^(-3),10^(-2),10^(-1)),repet)[1:nb.nnet]

# target is the set of target summary stats
# x is the parameter vector (long vector of numbers from the simulations) and is the dependent variable for the regression. In the context of model selection, it is a number indicating which model we are working with
# sumstat is an array of simulated summary stats (i.e. independent variables).
# tol is the required proportion of points nearest the target values
# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
# if noweight=T, no Epanechnikov weights are calculated

# if rejmethod=T it doesn't bother with the neural networks, and just compute the smooth ratio of the acceptance ratios

# If rejmethod=F it returns a list with the following components:-

# $x regression adjusted values
# $vals - unadjusted values in rejection region (i.e. normal rejection)
# $wt - the regression weight (i.e. the Epanechnikov weight)
# $ss - the sumstats corresponding to these points
# $predmean - estimate of the posterior mean
# $fv - the fitted value from the regression


if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))

nss <- length(sumstat[1,])


# scale everything

    scaled.sumstat <- sumstat

    for(j in 1:nss){

    	scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
    }
    target.s.tot <- target

    for(j in 1:nss){

    	target.s.tot[,j] <- normalise(target[,j],sumstat[,j][gwt])
    }

	rm(sumstat)

#Cross validation
the_pred.tot=NULL
for(pseu in 1:nrow(target.s.tot)){
	target.s=as.numeric(target.s.tot[pseu,])
	
	# calc euclidean distance
	    sum1 <- 0
	    for(j in 1:nss){
	    	sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
	   }
	   dst <- sqrt(sum1)
	# includes the effect of gwt in the tolerance
	    dst[!gwt] <- floor(max(dst[gwt])+10)


	# wt1 defines the region we're interested in
	    abstol <- quantile(dst,tol)
	    wt1 <- dst <= abstol

	    if(rejmethod){
		regwt <- 1-dst[wt1]^2/abstol^2
		the_pred<-NULL
		for (i in 1:length(unique(x)))
		{
			aux<-unique(x)[i]
		  	the_pred<-c(the_pred,sum(regwt[x[wt1]==aux]))
		}
		the_pred<-the_pred/sum(the_pred)
	    }
	    else{
		  regwt <- 1-dst[wt1]^2/abstol^2
		if(noweight)
		  	regwt <- rep(1,length(regwt))

		#Fit a neural network for predicting the class
		ll<-NULL
	    	for (i in 1:nb.nnet)
	    	{
			fit1 <- nnet(scaled.sumstat[wt1,],class.ind(as.factor(x))[wt1,],weights=regwt,decay=the_decay[i],size=size.nnet,softmax=T,maxit=500,trace=F)
			ll<-c(ll,list(fit1))
		}
		vect_pred<-NULL
		for (i in 1:nb.nnet)
	    	{
		  	vect_pred<-cbind(vect_pred,as.numeric(predict(ll[[i]],data.frame(rbind(target.s)))))
		}
		the_pred<-apply(vect_pred,FUN=mean,MARGIN=1)
	    }
	    the_pred
	    print(the_pred)
	    the_pred.tot=rbind(the_pred.tot, the_pred)
		write.table(the_pred.tot, file=output, col.names=F, row.names=F)
	   }
	   return(the_pred.tot)
}




