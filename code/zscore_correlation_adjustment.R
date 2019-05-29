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





tmp <- data.ext2 %>% group_by(meta.id) %>% mutate(n = n()) %>% filter(n > 9) %>% 
	filter(outcome.type != "rate" & !is.na(outcome.type)) # %>% filter(file.nr != 3014 & file.nr != 208)

# tmp <- tmp %>% filter(meta.id == 8361 | meta.id == 8363)
# tmp <- tmp %>% filter(meta.id == 378)

#Meta-analysis and ajustment part:
# meta.id.vector <- unique(tmp$meta.id)
meta.id.vector <- metac$meta.id
meta.id.vector <- meta.id.vector[-which(meta.id.vector == 159329)]
meta.analyses <- list()
meta.analyses.limitmeta <- list()
meta.analyses.copas <- list()
counter <- 0

for(u in meta.id.vector){
	counter <- counter + 1
	print(c(u,counter))
	meta.analyses[[counter]] <- mcor(data = tmp[tmp$meta.id == u,])
	meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
}


#List extraction:
correlation.meta.analysis.estimates <- cbind(
	meta.id = meta.id.vector,
	est.cor.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.fixed})),
	est.cor.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.random})),
	se.est.cor.fixef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.fixed})),
	se.est.cor.ranef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.random})),
	
	est.cor.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
	se.est.cor.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
	
	est.cor.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
	se.est.cor.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})))

meta.wcor <- merge(metac, zscore.meta.analysis.estimates, by = c("meta.id"))

#Meta-analysis and ajustment part:
# meta.id.vector <- unique(tmp$meta.id)
meta.id.vector <- metac$meta.id
meta.id.vector <- meta.id.vector[-which(meta.id.vector == 159329)]
meta.analyses <- list()
meta.analyses.limitmeta <- list()
meta.analyses.copas <- list()
counter <- 0

for(u in meta.id.vector){
	counter <- counter + 1
	print(c(u,counter))
	meta.analyses[[counter]] <- metagen(TE = z, seTE = sqrt(var.z), studlab = study.name, tmp[tmp$meta.id == u,])
	meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
}


#List extraction:
zscore.meta.analysis.estimates <- cbind(
	meta.id = meta.id.vector,
	est.z.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.fixed})),
	est.z.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.random})),
	se.est.z.fixef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.fixed})),
	se.est.z.ranef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.random})),
	
	est.z.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
	se.est.z.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
	
	est.z.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
	se.est.z.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})))




# correlation.meta.analysis.estimates <- data.frame(correlation.meta.analysis.estimates) %>% 
# 	mutate(est.z.fixef = 0.5 * log( (1 + est.cor.fixef)/(1 - est.cor.fixef) , base = exp(1)),
# 				 est.z.ranef = 0.5 * log( (1 + est.cor.ranef)/(1 - est.cor.ranef) , base = exp(1)),
# 				 est.z.reg = 0.5 * log( (1 + est.cor.reg)/(1 - est.cor.reg), base = exp(1)),
# 				 est.z.copas = 0.5 * log( (1 + est.cor.copas)/(1 - est.cor.copas) , base = exp(1)),
# 				 se.est.z.fixef = ifelse(total1 + total2 > 4))

meta.wcor <- merge(metac, zscore.meta.analysis.estimates, by = c("meta.id"))

save(meta.wcor, file =  file.path(PATH_RESULTS, "meta.z.RData"))
load(file.path(PATH_RESULTS, "meta.complete.RData"))









# #Examplary m.a. for test purpose: 134944 (bin), 158355 (cont)
# 
# bin <- data.ext2 %>% filter(meta.id == 134944)
# cont <- data.ext2 %>% filter(meta.id == 158355)
# bin2 <- data.ext2 %>% filter(meta.id == 16060)
# 
# 
# forest(meta.id(134944))
# forest(metacor(cor = cor.phi, n = total1 + total2, data = bin))
# forest(metacor(cor = cor.fisher, n = total1 + total2, data = bin))
# 
# 
# forest(rma.uni(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2,  
#                measure = "RR", method = "FE", data = bin)) #For comparison
# forest(rma.uni(yi = smd.pbit, vi = var.smd.pbit, measure = "PBIT", method = "FE", 
#                data = bin)) 
# forest(rma.uni(yi = cor.phi, vi = var.cor.phi, measure = "PHI", method = "FE", 
#                data = bin)) 
# forest(rma.uni(yi = cor.fisher, vi = var.cor.fisher, measure = "ZCOR", method = "FE", 
#                data = bin)) 
# 
# forest(rma.uni(yi = effect, sei = se, measure = "MD", method = "FE", data = cont))
# forest(rma.uni(m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, n1i = total1, n2i = total2, 
#                measure = "MD", method = "FE", data = cont))
# forest(rma.uni(yi = smd, vi = var.smd, measure = "SMD", method = "FE", data = cont))
# forest(rma.uni(yi = cor.pearson, vi = var.cor.pearson, measure = "COR", method = "FE", data = cont))
# forest(rma.uni(yi = cor.fisher, vi = var.cor.fisher, measure = "ZCOR", method = "FE", data = cont))



