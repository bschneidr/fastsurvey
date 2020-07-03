svyquantile<-function(x,design,quantiles,...) UseMethod("svyquantile", design)

svyquantile.survey.design <- function (x, design, quantiles, alpha = 0.05,
                                       interval.type = c("mean", "beta","xlogit", "asin","score"),
                                       na.rm = FALSE, 
                                       se = TRUE,
                                       qrule=c("math","school","shahvaish","hf1","hf2","hf3","hf4","hf7","hf8"),
                                       df = NULL, ...) {
    
   if (inherits(x, "formula")) 
        x <- model.frame(x, model.frame(design), na.action = na.pass)
    else if (typeof(x) %in% c("expression", "symbol")) 
        x <- eval(x, model.frame(design, na.action = na.pass))
    if (na.rm) {
        nas <- rowSums(is.na(x))
        design <- design[nas == 0, ]
        if (length(nas) > length(design$prob)) 
            x <- x[nas == 0, , drop = FALSE]
        else x[nas > 0, ] <- 0
    }

    if(is.null(df))
        df<-degf(design)

    if (is.character(qrule)){
        qrulename<-paste("qrule",match.arg(qrule),sep="_")
        qrule<-get(qrulename, mode="function")
    }

    qcrit<-if(df==Inf) qnorm else function(...) qt(...,df=df)
    

    interval.type<-match.arg(interval.type)
    if (interval.type=="score"){
        ci_fun<-ffullerCI
    } else {
        ci_fun<-function(...) woodruffCI(...,method=interval.type)
    }

    w<-weights(design)
    
    rvals<-lapply(x,
                  function(xi){ sapply(quantiles,
                                      function(p){
                                          qhat<-qrule(xi,w,p)
                                          ci<-ci_fun(xi,qhat,p,design,qrule,alpha,df)
                                          names(ci)<-c(signif(alpha/2,2),signif(1-alpha/2,2))
                                          list(quantile=qhat,ci=ci,
                                               se=diff(ci)/(2*qcrit(1-alpha/2))
                                               )
                                      }
                                      )
                  }
                  )
    
    class(rvals)<-"newsvyquantile"
    rvals
}

svyquantile.svyrep.design <- function (x, design, quantiles, alpha = 0.05,
                                       interval.type = c("mean", "beta","xlogit", "asin","quantile"),
                                       na.rm = FALSE, 
                                       se = TRUE,
                                       qrule=c("school","math","shahvaish","hf1","hf2","hf3","hf4","hf7","hf8"),
                                       df = NULL, return.replicates=FALSE,...) {
    interval <- match.arg(interval.type)
    
    if (design$type %in% c("JK1", "JKn") && interval == "quantile") 
        warning("Jackknife replicate weights may not give valid standard errors for quantiles")
    if (design$type %in% "other" && interval == "quantile") 
        warning("Not all replicate weight designs give valid standard errors for quantiles.")

    if (inherits(x, "formula")) 
        x <- model.frame(x, design$variables, na.action = if (na.rm) 
            na.pass
        else na.fail)
    else if (typeof(x) %in% c("expression", "symbol")) 
        x <- eval(x, design$variables)
    if (na.rm) {
        nas <- rowSums(is.na(x))
        design <- design[nas == 0, ]
        if (length(nas) > length(design$prob)) 
            x <- x[nas == 0, , drop = FALSE]
        else x[nas > 0, ] <- 0
    }

    if (is.character(qrule)){
        qrulename<-paste("qrule",match.arg(qrule),sep="_")
        qrule<-get(qrulename, mode="function")
    }
    
    if (is.null(df)) df <- degf(design)
    if (df == Inf)  qcrit <- qnorm else qcrit <- function(...) qt(..., df = df)

    w<-weights(design,"sampling")

    if (interval.type=="quantile"){
        ci_fun<-function(...) repCI(...,return.replicates=return.replicates)
    } else {
        if(return.replicates) warning("return.replicates=TRUE only implemented for interval.type='quantile'")
        ci<-woodruffCI
    }
    
    rvals<-lapply(x,
                  function(xi){ sapply(quantiles,
                                       function(p){
                                           qhat<-qrule(xi,w,p)
                                           ci<-ci_fun(xi,qhat,p,design,qrule,alpha,df)
                                           names(ci)<-c(signif(alpha/2,2),signif(1-alpha/2,2))
                                           list(quantile=qhat,ci=ci,
                                                se=diff(ci)/(2*qcrit(1-alpha/2))
                                                )
                                       }
                                       )
                  }
                  )
    
    rvals
}

SE.svyquantile<-function(object) {
    lapply(object, function(v) sapply(v, function(vp) vp$se))
    }

## Need to call qs() from qrule(), because definition of p varies by rule. <sigh>

qrule_math <-qrule_hf1 <- function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
        }
    qdata<-qs(x,w,p)
    if(qdata$wlow==0) qdata$qlow else qdata$qup
}

qrule_school<-qrule_hf2 <- function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
        }
    qdata<-qs(x,w,p)
    if(qdata$wlow==0) (qdata$qlow+qdata$qup)/2 else qdata$qup
}

qrule_hf3 <- function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    qdata<-qs(x,w,p)
    if((qdata$wlow==0)&& (qdata$ilow %%2 ==0)) qdata$qlow else qdata$qup
}

## modified version of hf6
qrule_shahvaish<-function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
    padj<-floor(p*n)/(n+1)
    qdata<-qs(x,w,padj)
    wj <- (qdata$wup+qdata$wlow)/sum(w)
    
    padj<- (floor(p*n)+0.5 -0.5*wj)/(n+1)
    qdata<-qs(x,w,padj)
    
    if(qdata$wlow==0) qdata$qlow else qdata$qup
}


qrule_hf4 <- function(x,w,p){
     if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
     }
     qdata<-qs(x,w,p)
     gamma<-with(qdata, wup/(wup+wlow))
     qdata$qlow*(1-gamma)+qdata$qup*gamma
}


qrule_hf7 <- function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
    j<-floor(p*n+1-p)
    padj<-j/n

    qdata<-qs(x,w,padj)
    
    gamma<-with(qdata, wup/(wup+wlow))
    qdata$qlow*(1-gamma)+qdata$qup*gamma
}

qrule_hf8 <- function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
    j<-floor(p*n+(p+1)/3)
    padj<-j/n
    
    qdata<-qs(x,w,padj)
    
    gamma<-with(qdata, wup/(wup+wlow))
    qdata$qlow*(1-gamma)+qdata$qup*gamma
}



qs <- function(x, w, p){
    ## already has missings removed, ties handled.
    n<-length(x)
    ii<-order(x)
    x<-x[ii]
    cumw<-cumsum(w[ii])

    pos<-max(which(cumw<=p*sum(w)))
    posnext<-if(pos==length(x)) pos else pos+1
    
    list(qlow=x[pos], qup=x[posnext], ilow=pos,iup=posnext, wlow=p-cumw[pos]/sum(w),wup=cumw[posnext]/sum(w)-p)
    
}
    

woodruffCI<-function(x, qhat,p, design, qrule,alpha,df,method=c("mean","beta","xlogit","asin")){

    method<-match.arg(method)
    m<-svymean(x<=qhat, design)
    names(m)<-"pmed"

    pconfint<-switch(method,
                    mean=as.vector(confint(m, 1, level = 1-alpha, df = df)),
                    xlogit= {xform <- svycontrast(m, quote(log(`pmed`/(1 - `pmed`))));
                        expit(as.vector(confint(xform, 1, level = 1-alpha,  df = df)))},
                    beta={n.eff <- coef(m) * (1 - coef(m))/vcov(m);
                        rval <- coef(m)[1]
                        n.eff <- n.eff * (qt(alpha/2, nrow(design) - 1)/qt(alpha/2, degf(design)))^2
                        c(qbeta(alpha/2, n.eff * rval, n.eff * (1 - rval) +  1), qbeta(1 - alpha/2, n.eff * rval + 1, n.eff *  (1 - rval)))
                    },
                    asin={xform <- svycontrast(m, quote(asin(sqrt(`pmed`))))
                        sin(as.vector(confint(xform, 1, level = 1-alpha,  df = df)))^2
                    }
                    )
    lower<-qrule(x,weights(design,"sampling"), pconfint[1])
    upper<-qrule(x,weights(design,"sampling"), pconfint[2])
    c(lower, upper)
    
}



ffullerCI<-function(x, qhat, p, design, qrule, alpha, df){
        qcrit<-if(df==Inf) qnorm else function(...) qt(...,df=df)
        U <- function(theta) {
            ((x > theta) - (1 - p))
        }
        scoretest <- function(theta, qlimit) {
            umean <- svymean(U(theta), design)
            umean/SE(umean) - qlimit
        }
        iqr <- IQR(x)
        lowerT <- min(x) + iqr/100
        upperT <- max(x) - iqr/100
        tol <- 1/(100 * sqrt(nrow(design)))
        
        qlow<- uniroot(scoretest, interval = c(lowerT, upperT), qlimit = qcrit(alpha/2, lower.tail = FALSE), tol = tol)$root
        qup<-uniroot(scoretest, interval = c(lowerT, upperT), qlimit = qcrit(alpha/2, lower.tail = TRUE), tol = tol)$root

        w<-weights(design)

        c(qrule(x, w, mean(x<=qlow)), qrule(x, w, mean(x<=qup)))
        
    }


repCI<-function(x, qhat, p, design, qrule, alpha, df,return.replicates){
    qcrit<-if(df==Inf) qnorm else function(...) qt(...,df=df)

    wrep<-weights(x,"analysis")
    reps<-apply(wrep,2, function(wi) qrule(x,wi,p))
    v<-with(design, svrVar(reps, scale=scale, rscales=rscales, mse=mse,coef=qhat))

    ci<- qhat+ c(-1,1)*sqrt(v)*qcrit(1-alpha/2)

    if (return.replicates) attr(ci,"replicates")<-reps

    ci
}
