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

load(file.path(PATH_RESULTS, file = "mly.RData"))
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

# metac %>% filter(meta.id == 30578) %>% select(est.trimfill.fixef, est.fixef, est.copas, est.reg.ranef)
# 
# print(meta.bin.complete(tmp, 10, sig.level = 0.1, sm = "RR") %>% select(doi, se.min, se.max, total.events,
# 																							logest.trimfill.fixef, est.fixef, logest.copas, logest.reg.ranef))
# 
# print(meta.bin.complete(tmp, 10, sig.level = 0.1, sm = "RR") %>% select(doi, se.min, se.max, total.events, logest.trimfill.fixef,
# 																							est.fixef, logest.copas, logest.reg.ranef))


########################################################################################################################
########################################################################################################################
#Trimfill lor comparison
########################################################################################################################
########################################################################################################################


#Scatterplots of metac-analysis and corrected log(abs(abs estimates
trimfill.fixef.sc <- metac %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.trimfill.fixef)))) + geom_point() + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	xlab("Fixed effects log estimate") + ylab("Trimfill log estimate") + ggtitle("Fixed effects to Trimfill")

trimfill.ranef.sc <- metac %>% ggplot(aes(x = log(abs((est.ranef))), y = log(abs(est.trimfill.fixef)))) + geom_point() + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	xlab("Random effects log estimate") + ylab("Trimfill log estimate") + ggtitle("Random effects to Trimfill")



########################################################################################################################
########################################################################################################################
#Copas lor comparison
########################################################################################################################
########################################################################################################################

#Scatterplots of metac-analysis and corrected log estimates
copas.fixef.sc <- metac %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.copas)))) + geom_point() + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	xlab("Fixed effects log estimate") + ylab("Copas log estimate") + ggtitle("Fixed effects to Copas") 

copas.ranef.sc <- metac %>% ggplot(aes(x = log(abs((est.ranef))), y = log(abs(est.copas)))) + geom_point() + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE, show.legend = T) + 
	xlab("Random effects log estimate") + ylab("Copas log estimate") + ggtitle("Random effects to Copas") 


########################################################################################################################
########################################################################################################################
#Regression lor comparison
########################################################################################################################
########################################################################################################################

#Scatterplots of metac-analysis and corrected log estimates
reg.fixef.sc <- metac %>% ggplot(aes(x = log(abs(log(abs((est.fixef))))), y = log(abs(log(abs(est.reg.ranef)))))) + geom_point() + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	xlab("Fixed effects log estimate") + ylab("Regression log estimate") + ggtitle("Fixed effects to Regression")

reg.ranef.sc <- metac %>% ggplot(aes(x = log(abs((est.ranef))), y = log(abs(est.reg.ranef)))) + geom_point() + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	xlab("Random effects log estimate") + ylab("Regression log estimate") + ggtitle("Random effects to Regression")





########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

#Scatterplots of test- statistics:
trimfill.fixef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.trimfill.fixef))) + geom_point() +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Fixed effects and trimfill z statistics") + ylab("Trimfill  statistic") + xlab("Fixed effects  statistic")

trimfill.ranef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.ranef), y = abs(zval.trimfill.fixef))) + geom_point() +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Random effects and trimfill z statistics") + ylab("Trimfill  statistic") + xlab("Random effects  statistic")

#Scatterplots of test- statistics:
copas.fixef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.copas))) + geom_point() +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Fixed effects and Copas z statistics") + ylab("Copas  statistic") + xlab("Fixed effects  statistic")

copas.ranef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.ranef), y = abs(zval.copas))) + geom_point() +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Random effects and Copas z statistics") + ylab("Copas  statistic") + xlab("Random effects  statistic")

#Scatterplots of test- statistics:
reg.fixef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.reg.ranef))) + geom_point() +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Fixed effects and Regression z statistics") + ylab("Regression  statistic") + xlab("Fixed effects  statistic")

reg.ranef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.ranef), y = abs(zval.reg.ranef))) + geom_point() +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Random effects and Regression z statistics") + ylab("Regression  statistic") + xlab("Random effects  statistic")




#####


p.missing.copas <- metac %>% ggplot(aes(x = missing.copas/n)) + geom_histogram(bins = 20) + 
	xlim(0, 0.5) + ylim(0, 200) + xlab("Proportion of missing studies") + ggtitle("Copas selection model") + theme_bw()
p.missing.trim <- metac %>% ggplot(aes(x = missing.trim)) + geom_histogram(bins = 20) + 
	xlim(0, 0.5) + ylim(0, 200) + xlab("Proportion of missing studies") + ggtitle("Trim-and-fill") + theme_bw()
grid.arrange(p.missing.copas, p.missing.trim, ncol = 2) 


