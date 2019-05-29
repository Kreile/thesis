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

file.dat <- "data.processed.RData"
if (file.exists(file.path(PATH_RESULTS, file.dat))) {
	load(file.path(PATH_RESULTS, file.dat))
} else {
	data.ext2 = pb.process2(data)
	save(data.ext2, file =  file.path(PATH_RESULTS, file.dat))
}

file.bin <- "pb.bin.RData"
if (file.exists(file.path(PATH_RESULTS, file.bin))) {
	load(file.path(PATH_RESULTS, file.bin))
} else {
	meta.bin <- meta.bin.complete(data.ext, min.study.number = 10, sig.level = 0.1, sm = "RR")
	save(meta.bin, file =  file.path(PATH_RESULTS, file.bin))
}

file.cont <- "pb.cont.RData"
if (file.exists(file.path(PATH_RESULTS, file.cont))) {
	load(file.path(PATH_RESULTS, file.cont))
} else {
	meta.cont <- meta.cont.complete(data.ext, min.study.number = 10, sig.level = 0.1)
	save(meta.cont, file =  file.path(PATH_RESULTS, file.cont))
}

file.meta <- "meta.RData"
if (file.exists(file.path(PATH_RESULTS, file.meta))) {
	load(file.path(PATH_RESULTS, file.meta))
} else {
	meta <- pb.meta.merge(meta.bin, meta.cont)
	save(meta, file =  file.path(PATH_RESULTS, file.meta))
}

load(file.path(PATH_RESULTS, "meta.complete.RData"))


#Applying test and adjustment criteria:
metac.bin <- meta.bin %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac.cont <- meta.cont %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac <- meta %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)





#Due to issues with correlations in binary cases, try to use metafor package:
#New Data process function

#Import data
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

if (file.exists(file.path(PATH_RESULTS, "data.RData"))) {
  load(file.path(PATH_RESULTS, "data.RData"))
} else {
  data = pb.readData(path = PATH_DATA, file = FILE)
  tmp = pb.clean(data)
  data = tmp[[1]]
  aliases = tmp[[2]]
  save(data, file =  file.path(PATH_RESULTS, "data.RData"))
}

if (file.exists(file.path(PATH_RESULTS, "data.ext.RData"))) {
  load(file.path(PATH_RESULTS, "data.ext.RData"))
} else {
  data.ext <- pb.process(data)
  save(data.ext, file =  file.path(PATH_RESULTS, "data.ext.RData"))
}



require(metafor)

#Metafor package part

mtd <- data.ext %>% filter(meta.id == 101436)
mtd$n <- mtd$total1 + mtd$total2
mtd$var.log.OR <- mtd$se.log.OR^2
mtd.c <- data.ext %>% filter(meta.id == 122943)

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



