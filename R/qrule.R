
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
    if (is.na(padj)) return(NaN)
    qdata<-qs(x,w,padj)
    
    if(qdata$wlow==0) qdata$qlow else qdata$qup
}

qrule_hf6<-function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
    padj<-floor(p*n)/(n+1)
    qdata<-qs(x,w,padj)
    gamma<-with(qdata, wlow/(wup+wlow))
    qdata$qlow*(1-gamma)+qdata$qup*gamma
}


qrule_hf4 <- function(x,w,p){
     if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
     }
     qdata<-qs(x,w,p)
     gamma<-with(qdata, wlow/(wup+wlow))
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
    
    gamma<-with(qdata, wlow/(wup+wlow))
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
    
    gamma<-with(qdata, wlow/(wup+wlow))
    qdata$qlow*(1-gamma)+qdata$qup*gamma
}

last<-function(a) {
    if (any(a)) max(which(a)) else -Inf
    }

qs <- function(x, w, p){
    ## already has missings removed, ties handled.
    n<-length(x)
    ii<-order(x)
    x<-x[ii]
    cumw<-cumsum(w[ii])

    pos<-last(cumw<=p*sum(w))
    posnext<-if(pos==length(x)) pos else pos+1
    
    list(qlow=x[pos], qup=x[posnext], ilow=pos,iup=posnext, wlow=p-cumw[pos]/sum(w),wup=cumw[posnext]/sum(w)-p)
    
}
