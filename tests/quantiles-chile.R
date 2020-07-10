library(survey)
datos <- readRDS("datos_ejemplo.rds")

design <- svydesign(id = ~id_directorio, strata = ~estrato, weights = ~f_pers, check.strata = TRUE, data = datos)
set.seed(234262762)
repdesign <- as.svrepdesign(design, type = "subbootstrap", replicates=20)
options(survey.lonely.psu="remove")

values<-datos$ing_t_p[datos$CL_GRUPO_OCU_08=="ISCO08_6"]

f0<-svyquantile(~ing_t_p, subset(design,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5), qrule="math")
f0.5<-svyquantile(~ing_t_p, subset(design,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5), qrule="school")
f1<-svyquantile(~ing_t_p, subset(design,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5),qrule="hf2")

all.equal(c(values[1],mean(values),values[2]), c(f0,f0.5,f1))

suppressWarnings({
rf0<-svyquantile(~ing_t_p, subset(repdesign,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5), qrule="math")
rf0.5<-svyquantile(~ing_t_p, subset(repdesign,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5), qrule="school")
rf1<-svyquantile(~ing_t_p, subset(repdesign,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5),qrule="hf2")
})

all.equal(c(values[1],mean(values),values[2]), c(rf0,rf0.5,rf1))
