library(survey)
library(survival)
set.seed(2021-6-25)
test1 <- list(time=c(4,3,1,1,2,2,3),
              status=c(1,1,1,0,1,1,0),
              x1=as.factor(rbinom(7, 2, 0.5)),
              x=c(0,2,1,1,1,0,0))
# Fit a stratified model
mod_c <- coxph(Surv(time, status) ~ x1 + x, test1)
mod_d <- coxph(Surv(time, status) ~ x + x1, test1)
stopifnot(all.equal(regTermTest(mod_c, ~x1, df = Inf)[c("chisq","df","test.terms","p")],
regTermTest(mod_d, ~x1, df = Inf)[c("chisq","df","test.terms","p")]))

data(pbc, package="survival")

pbc$randomized<-with(pbc, !is.na(trt) & trt>0)
biasmodel<-glm(randomized~age*edema,data=pbc,family=binomial)
pbc$randprob<-fitted(biasmodel)
if (is.null(pbc$albumin)) pbc$albumin<-pbc$alb ##pre2.9.0

dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))
library(splines)
model<-svycoxph(formula = Surv(time, status > 0) ~ bili + protime + albumin+ns(bili,4)[,1:3], design = dpbc)
test<-regTermTest(model, ~ns(bili,4)[,1:3],method="LRT")
stopifnot(all.equal(test$chisq, 47.314, tolerance=1e-4))
stopifnot(all.equal(test$lambda, c(1.4764260, 1.0109836, 0.6923415),tolerance=1e-4))
