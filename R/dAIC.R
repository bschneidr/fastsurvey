

logLik.svyglm<-function(object,...){
   warning("svyglm not fitted by maximum likelihood.")
   object$deviance/-2
}

AIC.svyglm<-function(object,...,k=2){
	if (length(list(...))){
		do.call(rbind,lapply(list(object,...),extractAIC,k=k))
    } else {
	   extractAIC(object,k=k)
    }
}
extractAIC.svyglm<-function(fit,scale,k=2,...){
	if (length(attr(terms(fit),"factors"))){
	    r<-regTermTest(fit, delete.response(formula(fit)), method="LRT")
	    deltabar<-mean(r$lambda)
	} else {
	    r<-list(lambda=0)
	    deltabar<-NaN
	}
	d<-fit$deviance
	c(eff.p=sum(r$lambda), AIC=d+k*sum(r$lambda),deltabar=deltabar)
}

extractAIC.svrepglm<-extractAIC.svyglm

BIC.svyglm<-function(object,...,maximal){
	if (length(list(...))){
		do.call(rbind,lapply(list(object,...),dBIC,modelM=maximal))
    } else {
	   dBIC(object,modelM=maximal)
    }
	
}
	
dBIC<-function(modela,modelM){
	pm<-modela$rank
	pM<-modelM$rank	

	if (any(!(names(coef(modela))%in% names(coef(modelM))))){
		stop("coefficients in model but not in maximal model")
		}
	index<-!(names(coef(modelM))%in% names(coef(modela)))
	n<-1+modela$df.null	
	if(any(index)){
		wald<-coef(modelM)[index]%*%solve(vcov(modelM)[index,index],coef(modelM)[index])
		detDelta<-det(solve(modelM$naive.cov[index,index,drop=FALSE],modelM$cov.unscaled[index,index,drop=FALSE]))
		dbar<-detDelta^(1/(pM-pm))
		nstar<-n/dbar	
	}else {
		wald<-0
		detDelta<-1
		dbar<-1
		nstar=NaN
		}
	c(p=pm, BIC=wald+pm*log(n)+log(detDelta)+deviance(modelM),neff=nstar)
    }

extractAIC.svrepcoxph<-function (fit, scale, k = 2, ...) .NotYetImplemented()
extractAIC.svycoxph<-function (fit, scale, k = 2, ...) 
{
    Delta<-solve(fit$inv.info, fit$var)
    deltabar <- mean(diag(Delta))
    d <- -2*fit$ll[1]
    c(eff.p = sum(diag(Delta)), AIC = d + k * sum(diag(Delta)), deltabar = deltabar)
}
AIC.svycoxph<-function (object, ..., k = 2) 
{
    if (length(list(...))) {
        do.call(rbind, lapply(list(object, ...), extractAIC, 
            k = k))
    }
    else {
        extractAIC(object, k = k)
    }
}
