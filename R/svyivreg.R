
svyivreg<-function(formula, design, ...) UseMethod("svyivreg",design)

svyivreg.survey.design<-function(formula, design,...){

    .data<-model.frame(design)
    .data$.weights<-weights(design,"sampling")
    fitter<-get("ivreg", mode="function")
    .weights<-NULL ## make CMD check happy
    model<- fitter(formula, data=.data, weights=.weights)

    U<-estfun(model)/weights(design,"sampling")
    n<-NROW(U)
    infl<- U%*%model$cov.unscaled 
    v<-vcov(svytotal(infl,  design))
    
    model$invinf<-model$cov.unscaled
    model$cov.unscaled<-v
    model$df.residual<-degf(design)+1-length(coef(model))
    model$sigma<-model$sigma/sqrt(mean(weights(design,"sampling")))
    model$call<-sys.call(-1)
    class(model)<-c("svyivreg","ivreg")
    model
}


summary.svyivreg<-function(object, df = NULL, ...){
    class(object)<-"ivreg"
    summary(object, vcov.=NULL, df=df, diagnostics=FALSE,...)
}

vcov.svyivreg<-function(object,...) object$cov.unscaled

svyivreg.svyrep.design<-function(formula, design,return.replicates=FALSE,...){
    fitter<-get("ivreg", mode="function")
    .pweights<-NULL ## make CMD check happy

    withReplicates(design, return.replicates=return.replicates,
                   function(.weights, .data){
                       .data$.pweights<-.weights
                       m<-fitter(formula,data= .data, weights=.pweights)
                       coef(m)
                   })
                      
    }
