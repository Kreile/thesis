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

load(file.path(PATH_RESULTS, file = "data.processed.RData"))

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


data.ext2 <- pb.process2(data)


#Meta filtering: 
metac.bin <- meta.bin %>% filter(n.sig.type.bin > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac.cont <- meta.cont %>% filter(n.sig.type.cont > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac <- meta %>% filter(n.sig.type > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)


#Pb test stats agreement
require(GGally)
topairsplot.pb.stats.bin <- metac %>% select(stat.egger, stat.harbord, stat.peter, 
																						stat.schwarzer, stat.rucker, harbord.test)
#topairsplot.pb.stats.bin <- abs(topairsplot.pb.stats.bin)
ggpairs(topairsplot.pb.stats.bin, columns = 1:5, aes(colour = factor(harbord.test), alpha = 0.5))
ggpairs(topairsplot.pb.stats.bin, columns = 1:5, aes(alpha = 0.5))

topairsplot.pb.stats.cont <- metac %>% select(stat.egger, stat.thomson, stat.begg,
																							egger.test)
#topairsplot.pb.stats.bin <- abs(topairsplot.pb.stats.bin)
ggpairs(topairsplot.pb.stats.cont, columns = 1:3, aes(colour = factor(egger.test)))
ggpairs(topairsplot.pb.stats.cont, columns = 1:3, aes(alpha = 0.7))

#$Complete dataset pairsplot:
topairsplot <- metac %>% select(comparison.nr, pval.Q, pval.thomson, n.sig.type,	missing.copas,
																 G.squared, n, var.samplesize,mean.publication.year
																)
topairsplot <- topairsplot %>% filter(var.samplesize < 100000 & comparison.nr < 20)

ggpairs(topairsplot, aes(alpha = 0.1))



topairsplot.pb.stats.bin <- meta %>% select(stat.egger, stat.harbord, stat.peter, 
																						stat.schwarzer, stat.rucker, harbord.test)
ggpairs(topairsplot.pb.stats, aes(colour = factor(harbord.test)))


meta.bin %>% ungroup() %>% 
	select(stat.egger, stat.harbord, stat.peter, stat.schwarzer, stat.rucker) %>% 
	gather(key = "test.type", value = "test.statistic") %>% 
	ggplot()

topairsplot.pb.stats.cont <- meta %>% select(stat.egger, stat.thomson, stat.begg)
pairs(topairsplot.pb.stats.cont, lower.panel = NULL, gap = 0)


meta.bin %>% ungroup() %>% 
	select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
	gather(key = "test.type", value = "null.hypothesis")
