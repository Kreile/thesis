rm(list = ls())
PATH_HOME = path.expand("~") # user home
PATH = file.path(PATH_HOME, 'Data/PubBias')
PATH2 = file.path(PATH_HOME, 'PubBias')
FILE = 'cochrane_2018-06-09.csv'
PATH_DATA = file.path(PATH, 'data')
PATH_CODE = file.path(PATH2, 'code')
PATH_RESULTS = file.path(PATH2, 'results')
PATH_FIGURES = file.path(PATH_RESULTS, 'figures')

file_results = "pb.RData"

source(file.path(PATH_CODE, 'PubBias_functions.R'))

file.dat <- "data.RData"
if (file.exists(file.path(PATH_RESULTS, file.dat))) {
  load(file.path(PATH_RESULTS, file.dat))
} else {
  data = pb.readData(path = PATH_DATA, file = FILE)
  tmp = pb.clean(data)
  data = tmp[[1]]
  aliases = tmp[[2]]
  save(data, file =  file.path(PATH_RESULTS, file.dat))
}

load(file.path(PATH_RESULTS, "meta_analyses_summary_complete.RData"))
load(file.path(PATH_RESULTS, "meta.bin.RData"))
load(file.path(PATH_RESULTS, "meta.cont.RData"))
load(file.path(PATH_RESULTS, "meta.surv.RData"))
load(file.path(PATH_RESULTS, "data_used_for_analysis.RData"))





#Metafor package part

mtd <- data.ext2 %>% filter(meta.id == 101436)
mtd$n <- mtd$total1 + mtd$total2
mtd$var.log.OR <- mtd$se.log.OR^2
mtd.c <- data.ext2 %>% filter(meta.id == 122943)

mtd <- escalc(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2, 
       measure = "RTET", data = mtd, append = T, var.names = c("rtet", "var.rtet"))
mtd <- escalc(ri = rtet, ni = n, measure = "ZCOR", data = mtd, append = T, 
              var.names = c("zcor", "var.zcor"))
mtd <- escalc(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2, 
              measure = "OR", data = mtd, append = T, var.names = c("OR", "var.OR"))
mtd <- escalc(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2, 
              measure = "PBIT", data = mtd, append = T, var.names = c("pbit", "var.pbit"))
mtd <- escalc(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2, 
              measure = "OR2DN", data = mtd, append = T, var.names = c("or2dn", "var.or2dn"))
mtd <- escalc(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2, 
              measure = "OR2DL", data = mtd, append = T, var.names = c("or2dl", "var.or2dl"))
mtd <- escalc(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2, 
              measure = "OR", data = mtd, append = T, var.names = c("OR", "var.OR"))
mtd <- escalc(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2, 
              measure = "PHI", data = mtd, append = T, var.names = c("phi", "var.phi"))
#Conversion of a probit SMD according to "Introduction to meta-analysis". Generally, no correspondance btw. methods..
mtd <- mtd %>% mutate(a = ((total1 + total2)^2)/(total1*total2),
                      newr  = pbit/sqrt((pbit^2) + a),
                      newr.var = ((a^2)*var.pbit)/(((pbit^2)+a)^3))

#Compare forestplots: CI's and relations should remain equal after conversions. 
#Odds ratios:
forest(rma.uni(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2,  
               measure = "OR", method = "FE", data = mtd))
forest(rma.uni(yi = log.OR, vi = var.log.OR, measure = "RR", method = "FE", 
               data = mtd)) #Own calculation, little different

# meta.cor.bin <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, 
#                         method = "Inverse", method.tau = "REML", data = mtd, sm = "RR")
# forest(meta.cor.bin)

#Correlations
forest(rma.uni(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2,  
               measure = "RR", method = "FE", data = mtd)) #For comparison
forest(rma.uni(yi = phi, vi = var.phi, measure = "PHI", method = "FE", 
               data = mtd)) #PHI correlation, similar
forest(rma.uni(yi = fishersz, vi = fishersz.variance, measure = "ZCOR", method = "FE", 
               data = mtd)) #Own fishers-z calculation, completely wrong
forest(rma.uni(yi = rtet, vi = var.rtet, method = "FE", data = mtd, 
               measure = "RTET")) #Tetrachronic correlation, also quite different
forest(rma.uni(yi = zcor, vi = var.zcor, method = "FE", measure = "ZCOR", 
               data = mtd)) #fishers corr from tetrachcronic correlation, same as before
forest(rma.uni(yi = newr, vi = newr.var, method = "FE", 
               data = mtd)) #By hand calculation of correlation from PBIT SMD, very similar to RTET
forest(rma.uni(ri = correlation, ni = n, measure = "COR", method = "FE", 
               data = mtd)) #By hand calculation of correlation from lor and md as in "Intro.."


#Standardized mean differences - All quite good:
forest(rma.uni(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2,  
               measure = "RR", method = "FE", data = mtd))
forest(rma.uni(yi = pbit, vi = var.pbit, measure = "PBIT", method = "FE", data = mtd))
forest(rma.uni(yi = or2dn, vi = var.or2dn, 
               measure = "OR2DN", method = "FE", data = mtd))
forest(rma.uni(yi = or2dl, vi = var.or2dl, 
               measure = "OR2DL", method = "FE", data = mtd))















#Compare "Rosenstein. et. al. to metafor functions:
mtd <- data.ext2 %>% filter(meta.id == 92458)
mtd <- escalc(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2, 
              measure = "OR", data = mtd, append = T, var.names = c("lor", "var.lor"))

mtd <- mtd %>% 
  mutate(smd2 = lor * (sqrt(3)/pi),
         var.smd2 = var.lor * (3/(pi^2)),
         a = ((total1 + total2)^2)/(total2*total1),
         corr2 = smd2/sqrt((smd2^2) + a),
         var.corr2 = ((a^2)*var.smd2)/((smd2^2) + a)^3
       # fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
       # fishersz.variance = 1/(total1 + total2 - 3),
       # pval.type = 2*(1- pnorm(abs(log.OR/se.log.OR)))
       )

forest(rma.uni(yi = lor, vi = var.lor,
               measure = "OR", method = "FE", data = mtd))
forest(rma.uni(yi = smd2, vi = var.smd2, measure = "SMD",
               method = "FE", data = mtd))
forest(rma.uni(yi = smd.ordl, vi = var.smd.ordl, measure = "SMD",
               method = "FE", data = mtd))
forest(rma.uni(yi = corr2, vi = var.corr2, measure = "COR",
               method = "FE", data = mtd))
forest(rma.uni(yi = cor.phi, vi = var.cor.phi, measure = "COR",
               method = "FE", data = mtd))
forest(rma.uni(yi = cor.phi, vi = var.cor.phi, measure = "COR",
               method = "FE", data = mtd))

# forest(rma.uni(yi = corr2, vi = var.corr2, measure = "COR", 
#                method = "FE", data = mtd))
# forest(rma.uni(yi = cor.phi, vi = var.cor.phi, measure = "PHI", 
#                method = "FE", data = mtd))
# forest(rma.uni(yi = smd.pbit, vi = var.smd.pbit, measure = "PBIT", 
#                method = "FE", data = mtd))


