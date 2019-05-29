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





comp <- metac %>% select(meta.id, comparison.name) #%>% distinct(comparison.name)

comp$isin <- mapply(FUN = grepl, x = comp$comparison.name, pattern = "sensitivity")
length(which(comp$isin == TRUE))
comp$isin2 <- mapply(FUN = grepl, x = comp$comparison.name, pattern = "Sensitivity")
length(which(comp$isin2 == TRUE))

meta.wout.sensitivity <- metac %>% filter(!grepl(pattern = "sensitivity", comparison.name) | !grepl(pattern = "Sensitivity", comparison.name))


comp <- metac %>% select(meta.id, outcome.name) #%>% distinct(outcome.name)

comp$isin <- mapply(FUN = grepl, x = comp$outcome.name, pattern = "sensitivity")
length(which(comp$isin == TRUE))
comp$isin2 <- mapply(FUN = grepl, x = comp$outcome.name, pattern = "Sensitivity")
length(which(comp$isin2 == TRUE))

meta.wout.sensitivity <- metac %>% filter(!grepl(pattern = "sensitivity", outcome.name) | !grepl(pattern = "Sensitivity", outcome.name))


comp <- metac %>% select(meta.id, sungroup.name) #%>% distinct(sungroup.name)

comp$isin <- mapply(FUN = grepl, x = comp$sungroup.name, pattern = "sensitivity")
length(which(comp$isin == TRUE))
comp$isin2 <- mapply(FUN = grepl, x = comp$sungroup.name, pattern = "Sensitivity")
length(which(comp$isin2 == TRUE))

meta.wout.sensitivity <- metac %>% filter(!grepl(pattern = "sensitivity", sungroup.name) | !grepl(pattern = "Sensitivity", sungroup.name))
