
fisherinf<-function(beta, model, X=model.matrix(model)){
    design<-model$survey.design
    eta<-X %*%beta
    mu<-drop(family(model)$linkinv(eta))
    V<-family(model)$variance(mu)
    d<-drop(family(model)$mu.eta(eta))
    w<-weights(design,"sampling")
    D<-X*d
    t(D)%*%(w*D/V)
}



scores<-function(model, newbeta, X=model.matrix(model)){
    design<-model$survey.design
    out<-newbeta==0
    eta<-X %*%newbeta
    Xin<-X[, !out, drop=FALSE]
    Xout<-X[, out, drop=FALSE]
    Xout_in<-qr.resid(qr(Xin),Xout)
    
    mu<-drop(family(model)$linkinv(eta))
    V<-family(model)$variance(mu)
    d<-drop(family(model)$mu.eta(eta))
    Y<-model$y
    Uin<-Xin*((d/V)*(Y-mu))
    Uout<-Xout*((d/V)*(Y-mu))
    
    I_rss98<-fisherinf(newbeta, model, X) ## equation 3.3, Rao, Scott, and Skinner 1998
    Atwiddle<-I_rss98[out,!out]%*%solve(I_rss98[!out,!out]) #equation 3.2, ibid
    Uout_in<-Uout-Uin%*%t(Atwiddle)
    
    svytotal( Uout_in, design)
}


pseudoscore<-function(model,newbeta,ddf, X=model.matrix(model)){

    s<-scores(model, newbeta,X)
    np<-length(coef(s))
    Qtest<-coef(s)%*%solve(vcov(s),coef(s))
    if (is.null(ddf)) ddf<-model$df.resid
    if (ddf<=0) {
        warning("negative or zero ddf set to 1")
        ddf<-1
    }
    ndf<-qr(vcov(s))$rank
    if(ndf<np) warning(paste("df reduced from",np,"to",ndf,"due to linear dependencies"))
    if (is.finite(ddf)) {
        p<-pf(Qtest,ndf,ddf,lower.tail=FALSE)
    } else{
       p<- pchisq(Qtest,ndf,lower.tail=FALSE)
        }
    c(X2=Qtest, df=ndf, ddf=ddf, p=p)
}


workingscore<-function(model,newbeta,ddf, lrt.approximation, X=model.matrix(model)){

    s<-scores(model, newbeta,X)
    I<-fisherinf(newbeta,model,X)
    out<-newbeta==0
    V<-I[out,out,drop=FALSE]-I[out,!out,drop=FALSE]%*%solve(I[!out,!out,drop=FALSE],I[!out,out,drop=FALSE])
    phi<-coef(summary(model)$dispersion)
    mwt<-mean(weights(model$survey.design))
    Qtest<-coef(s)%*%solve(V,coef(s))/phi/mwt
    Vmis<-solve(V,vcov(s))/phi/mwt

    if (is.null(ddf))  ddf<-model$df.resid

    if (ddf<=0) {
        warning("negative or zero ddf set to 1")
        ddf<-1
    }

    
    if (is.finite(ddf)) {
        p<-pFsum(Qtest,a=eigen(Vmis)$values, df=rep(1,NROW(Vmis)),ddf=ddf,lower.tail=FALSE, method=lrt.approximation)
    } else {
       p<- pchisqsum(Qtest,a=eigen(Vmis)$values, df=rep(1,NROW(Vmis)),lower.tail=FALSE, method=lrt.approximation)
    }

    
    c("Rao-Scott X^2"=Qtest/sum(diag(Vmis)),ddf=ddf, p=p)
}

svyscoretest<-function(model, drop.terms=NULL, add.terms=NULL, method=c("working","pseudoscore","individual"),
                       ddf=NULL,lrt.approximation="satterthwaite", ...){
    UseMethod("svyscoretest", model)
}


scoretest_addterms.svyglm<-function(model,add.terms, method=method,ddf=ddf, lrt.approximation=lrt.approximation,...){
    Xin<-model.matrix(model)
    design<-model$survey.design
    bigformula<-eval(bquote(update(.(formula(model)), .~.+.(add.terms[[-1]]))))
    Xbig<-model.matrix(bigformula, model.frame(design))
    added<-!(colnames(Xbig) %in% colnames(Xin))
    newbeta<-rep(0,NCOL(Xbig))
    idx<-match(colnames(Xin), colnames(Xbig))
    newbeta[idx]=coef(model)  
    s<-scores(model, newbeta, Xbig)
    
    switch(method, 
           working=workingscore(model,newbeta, lrt.approximation=lrt.approximation,ddf=ddf, X=Xbig),
           pseudoscore=pseudoscore(model,newbeta,ddf=ddf, X=Xbig),
           individual=s)
    }


svyscoretest.svyglm<-function(model,drop.terms=NULL, add.terms=NULL, method=c("working","pseudoscore","individual"),
                              ddf=NULL,lrt.approximation="satterthwaite",...){
    if(!xor(is.null(drop.terms),is.null(add.terms))){
        stop("One of 'drop.terms' and 'add.terms' must be non-NULL")
    }
    method<-match.arg(method)
    
    canonicalOrder<-function(term){
        tt<-strsplit(term,":")
        tt<-lapply(tt,sort)
        sapply(tt,paste,collapse=":")
    }
    
    if(inherits(drop.terms,"formula")){
        test_intercept<-explicit1(drop.terms)
        drop.terms<-attr(terms(drop.terms),"term.labels")
    } else{
        test_intercept<-FALSE
        return(scoretest_addterms.svyglm(model,add.terms, method=method,ddf=ddf, lrt.approximation=lrt.approximation,...))
    }
    
    okbeta<-!is.na(coef(model,na.rm=FALSE)) ## na.rm for svyglm
    tt<-attr(terms(model),"term.labels")
    aa<-attr(model.matrix(model),"assign")[okbeta]

    
    index<-which(aa %in% match(canonicalOrder(drop.terms),canonicalOrder(tt)))
    if (any(is.na(index)))
        stop("Terms didn't match:",canonicalOrder(drop.terms),canonicalOrder(tt))
    
    if (test_intercept){
        if (attr(terms(model),"intercept"))
            index<-unique(c(1,index))
        else
            stop("model does not have an intercept")
    }
      if (test_intercept) {
            test.formula <- make.formula(c(1, drop.terms))[[2]]
        }
        else {
            test.formula <- make.formula(drop.terms)[[2]]
        }
        if (!("formula" %in% names(model$call))) 
            names(model$call)[[2]] <- "formula"
   
    model0 <- eval(bquote(update(.(model), . ~ . - (.(test.formula)))), 
                   environment(formula(model)))
    
    beta<-coef(model)*0
    beta[-index]<-coef(model0)

    switch(method, 
           working=workingscore(model,beta, lrt.approximation=lrt.approximation,ddf=ddf),
           pseudoscore=pseudoscore(model,beta,ddf=ddf),
           individual=scores(model, beta))

}
