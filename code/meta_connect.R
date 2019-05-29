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

data.ext2 <- pb.process2(data)

meta.surv <- meta.surv.complete2(data = data.ext2, sig.level = 0.1, min.study.number = 10)
meta.bin <- meta.bin.complete2(data = data.ext2, sig.level = 0.1, min.study.number = 10, sm1 = "RR")
meta.cont <- meta.cont.complete2(data = data.ext2, sig.level = 0.1, min.study.number = 10)

meta <- bind_rows(meta.bin, meta.cont, meta.surv)

metac.bin <- meta.bin %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac.cont <- meta.cont %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac <- meta %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)



#Adding adjusted correlation and z-score treatment effects:
tmp <- data.ext2

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

meta.wcor <- merge(metac, correlation.meta.analysis.estimates, by = c("meta.id"))


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


meta.wcor <- merge(meta.wcor, zscore.meta.analysis.estimates, by = c("meta.id"))

save(meta.wcor, file =  file.path(PATH_RESULTS, "meta.complete.RData"))
load(file.path(PATH_RESULTS, "meta.complete.RData"))






