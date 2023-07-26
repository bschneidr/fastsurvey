svycontrast.svyvar<-function(stats, ...){
    s<-as.vector(as.matrix(stats))
    nms<-as.vector(outer(rownames(stats),colnames(stats),paste,sep=":"))
    v<-vcov(stats)
    names(s)<-nms
    dimnames(v)<-list(nms,nms)
    attr(s,"var")<-v
    attr(s,"statistic")<-"variance"
    class(s)<-"svystat"
    svycontrast(s,...)
}
