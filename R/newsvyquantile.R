svyquantile<-function(x,design,quantiles,...) UseMethod("svyquantile", design)



qrule_math <-qrule_hf1 <- function(qdata){
    if(qdata$wlow==0) qdata$qlow else qdata$qup
}

qrule_hf2 <- function(qdata){
    if(qdata$wlow==0) (qdata$qlow+qdata$qup)/2 else qdata$qup
}

qrule_hf3 <- function(qdata){
    if((qdata$wlow==0)&& (qdata$ilow %%2 ==0)) qdata$qlow else qdata$qup
}




qs <- function(x, w, p){
    ## already has missings removed, ties handled.
    n<-length(x)
    ii<-order(x)
    x<-x[ii]
    cumw<-cumsum(w[ii])

    pos<-max(which(cumw<=p*sum(w)))
    posnext<-if(pos==length(x)) pos else pos+1
    
    list(qlow=x[pos], qup=x[posnext], ilow=pos,iup=posnext, wlow=cumw[pos]/sum(w)-p,wup=cumw[posnext]/sum(w)-p)
   
}
    

woodruff<-function(x, med, design, qrule,method=c("mean","beta","xlogit","asin"),level,df){

    method<-match.arg(method)
    m<-svymean(x<=med, design)

    pconfint<-switch(method,
                    mean=as.vector(confint(m, 1, level = level, df = df)),
                    xlogit= {xform <- svycontrast(m, quote(log(`1`/(1 - `1`))));
                        expit(as.vector(confint(xform, 1, level = level,  df = df)))},
                    beta={n.eff <- coef(m) * (1 - coef(m))/vcov(m);
                        rval <- coef(m)[1]
                        alpha <- 1 - level
                        n.eff <- n.eff * (qt(alpha/2, nrow(design) - 1)/qt(alpha/2, degf(design)))^2
                        c(qbeta(alpha/2, n.eff * rval, n.eff * (1 - rval) +  1), qbeta(1 - alpha/2, n.eff * rval + 1, n.eff *  (1 - rval)))
                    },
                    asin={xform <- svycontrast(m, quote(asin(sqrt(`1`))))
                        sin(as.vector(confint(xform, 1, level = level,  df = df)))^2
                    }
                    )
    lower<-qrule(qs(x,weights(design,"sampling"), pconfint[1]))
    upper<-qrule(qs(x,weights(design,"sampling"), pconfint[2]))
    c(lower, upper)
    
}
