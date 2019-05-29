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

#Publication Bias Test Agreement:

meta %>% ggplot(aes(x = log(pval.fixef), y = log(pval.reg.ranef))) + geom_point(size = .5)
meta %>% ggplot(aes(x = (pval.fixef), y = (pval.reg.ranef))) + geom_point(size = .5)



meta.bin %>% mutate(n.sig = peter.test + rucker.test + egger.test + harbord.test + schwarzer.test) %>% 
  group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
  ggplot(aes(y = nn, x = n.sig)) + geom_col() 

meta.cont %>% mutate(n.sig = egger.test + thomson.test + begg.test) %>% 
  group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
  ggplot(aes(y = nn, x = n.sig)) + geom_col() 


########################################################################################################################
########################################################################################################################
#Trimfill lor comparison
########################################################################################################################
########################################################################################################################

#BINARY:

#Summaries:
metac.bin %>% mutate(xx = 1 - abs(logest.trimfill.fixef/log(est.fixef))) %>% filter(xx > 0) %>% ungroup() %>% 
  summarise(biased = sum(peter.test), larger.metac.effects = n())

metac.bin %>% mutate(xx = 1 - abs(logest.trimfill.fixef/log(est.fixef))) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(peter.test), smaller.metac.effects = n())

metac.bin %>% mutate(xx = 1 - abs(logest.trimfill.fixef/log(est.ranef))) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(peter.test), smaller.metac.effects.random = n())

#Histogramms of changes +- 100%
metac.bin %>% mutate(xx = 1 - abs(logest.trimfill.fixef/log(est.fixef))) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Trimfill")

metac.bin %>% mutate(xx = 1 - abs(logest.trimfill.fixef/log(est.fixef))) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Trimfill")

metac.bin %>% mutate(xx = 1 - abs(logest.trimfill.fixef/log(est.ranef))) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Trimfill")

metac.bin %>% mutate(xx = 1 - abs(logest.trimfill.fixef/log(est.ranef))) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Trimfill")

metac.bin %>% mutate(xx = 1 - abs(logest.trimfill.fixef/log(est.fixef))) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Trimfill")

metac.bin %>% mutate(xx = 1 - abs(logest.trimfill.fixef/log(est.ranef))) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Trimfill")

#Scatterplots of metac-analysis and corrected estimates
metac.bin %>% ggplot(aes(x = abs(log(est.fixef)), y = abs(logest.trimfill.fixef))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Fixed effects estimate") + ylab("Trimfill estimate") + ggtitle("Fixed effects to Trimfill")

metac.bin %>% ggplot(aes(x = abs(log(est.ranef)), y = abs(logest.trimfill.fixef))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Random effects estimate") + ylab("Trimfill estimate") + ggtitle("Random effects to Trimfill")

#Scatterplots of test-statistics:
metac.bin %>% ggplot(aes(x = zval.fixef, y = logest.trimfill.fixef/se.logest.trimfill.fixef)) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Fixed effects and trimfill test statistics") + ylab("Trimfill statistic") + xlab("Fixed effects statistic")

metac.bin %>% ggplot(aes(x = zval.ranef, y = logest.trimfill.fixef/se.logest.trimfill.fixef)) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Random effects and trimfill test statistics") + ylab("Trimfill statistic") + xlab("Random effects statistic")



#CONTINUOUS:

#Summaries:
meta.cont %>% mutate(xx = 1 - abs(est.trimfill.fixef.cont/est.fixef.cont)) %>% filter(xx > 0) %>% ungroup() %>% 
  summarise(biased = sum(egger.test), larger.metac.effects = n())

meta.cont %>% mutate(xx = 1 - abs(est.trimfill.fixef.cont/est.fixef.cont)) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(egger.test), smaller.metac.effects = n())

meta.cont %>% mutate(xx = 1 - abs(est.trimfill.fixef.cont/est.ranef.cont)) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(egger.test), smaller.metac.effects.random = n())

#Histogramms of changes +- 100%
meta.cont %>% mutate(xx = 1 - abs(est.trimfill.fixef.cont/est.fixef.cont)) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Trimfill")

meta.cont %>% mutate(xx = 1 - abs(est.trimfill.fixef.cont/est.fixef.cont)) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Trimfill")

meta.cont %>% mutate(xx = 1 - abs(est.trimfill.fixef.cont/est.ranef.cont)) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Trimfill")

meta.cont %>% mutate(xx = 1 - abs(est.trimfill.fixef.cont/est.ranef.cont)) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Trimfill")

meta.cont %>% mutate(xx = 1 - abs(est.trimfill.fixef.cont/est.fixef.cont)) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Trimfill")

meta.cont %>% mutate(xx = 1 - abs(est.trimfill.fixef.cont/est.ranef.cont)) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Trimfill")

#Scatterplots of metac-analysis and corrected estimates
meta.cont %>% ggplot(aes(x = abs(est.fixef.cont), y = abs(est.trimfill.fixef.cont))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Fixed effects estimate") + ylab("Trimfill estimate") + ggtitle("Fixed effects to Trimfill")

meta.cont %>% ggplot(aes(x = abs(est.ranef.cont), y = abs(est.trimfill.fixef.cont))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Random effects estimate") + ylab("Trimfill estimate") + ggtitle("Random effects to Trimfill")

#Scatterplots of test-statistics:
meta.cont %>% ggplot(aes(x = abs(zval.fixef.cont), y = abs(est.trimfill.fixef.cont/se.est.trimfill.fixef.cont))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Fixed effects and trimfill test statistics") + ylab("Trimfill statistic") + xlab("Fixed effects statistic")

meta.cont %>% ggplot(aes(x = abs(zval.ranef.cont), y = abs(est.trimfill.fixef.cont/se.est.trimfill.fixef.cont))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Random effects and trimfill test statistics") + ylab("Trimfill statistic") + xlab("Random effects statistic")

########################################################################################################################
########################################################################################################################
#Copas lor comparison
########################################################################################################################
########################################################################################################################

#Summaries
metac.bin %>% mutate(xx = 1 - abs(logest.copas/log(est.fixef))) %>% filter(xx > 0) %>% ungroup() %>% 
  summarise(biased = sum(peter.test), larger.metac.effects = n())

metac.bin %>% mutate(xx = 1 - abs(logest.copas/log(est.fixef))) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(peter.test), smaller.metac.effects = n())

metac.bin %>% mutate(xx = 1 - abs(logest.copas/log(est.ranef))) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(peter.test), smaller.metac.effects.random = n())

#Histogramms of changes +- 100%
metac.bin %>% mutate(xx = 1 - abs(logest.copas/log(est.fixef))) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Copas")

metac.bin %>% mutate(xx = 1 - abs(logest.copas/log(est.fixef))) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Copas")

metac.bin %>% mutate(xx = 1 - abs(logest.copas/log(est.ranef))) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Copas")

metac.bin %>% mutate(xx = 1 - abs(logest.copas/log(est.ranef))) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Copas")

metac.bin %>% mutate(xx = 1 - abs(logest.copas/log(est.fixef))) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Copas")

metac.bin %>% mutate(xx = 1 - abs(logest.copas/log(est.ranef))) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Copas")

#Scatterplots of metac-analysis and corrected estimates
metac.bin %>% ggplot(aes(x = abs(log(est.fixef)), y = abs(logest.copas))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Fixed effects estimate") + ylab("Copas estimate") + ggtitle("Fixed effects to Copas")

metac.bin %>% ggplot(aes(x = abs(log(est.ranef)), y = abs(logest.copas))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Random effects estimate") + ylab("Copas estimate") + ggtitle("Random effects to Copas")

#Scatterplots of test-statistics:
metac.bin %>% ggplot(aes(x = abs(zval.fixef), y = abs(logest.copas/se.logest.copas))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Fixed effects and copas test statistics") + ylab("copas statistic") + xlab("Fixed effects statistic")

metac.bin %>% ggplot(aes(x = abs(zval.ranef), y = abs(logest.copas/se.logest.copas))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Random effects and copas test statistics") + ylab("copas statistic") + xlab("Random effects statistic")



#CONTINUOUS:

#Summaries:
meta.cont %>% mutate(xx = 1 - abs(est.copas.cont/est.fixef.cont)) %>% filter(xx > 0) %>% ungroup() %>% 
  summarise(biased = sum(egger.test), larger.metac.effects = n())

meta.cont %>% mutate(xx = 1 - abs(est.copas.cont/est.fixef.cont)) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(egger.test), smaller.metac.effects = n())

meta.cont %>% mutate(xx = 1 - abs(est.copas.cont/est.ranef.cont)) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(egger.test), smaller.metac.effects.random = n())

#Histogramms of changes +- 100%
meta.cont %>% mutate(xx = 1 - abs(est.copas.cont/est.fixef.cont)) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to copas")

meta.cont %>% mutate(xx = 1 - abs(est.copas.cont/est.fixef.cont)) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to copas")

meta.cont %>% mutate(xx = 1 - abs(est.copas.cont/est.ranef.cont)) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to copas")

meta.cont %>% mutate(xx = 1 - abs(est.copas.cont/est.ranef.cont)) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to copas")

meta.cont %>% mutate(xx = 1 - abs(est.copas.cont/est.fixef.cont)) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to copas")

meta.cont %>% mutate(xx = 1 - abs(est.copas.cont/est.ranef.cont)) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to copas")

#Scatterplots of metac-analysis and corrected estimates
meta.cont %>% ggplot(aes(x = abs(est.fixef.cont), y = abs(est.copas.cont))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Fixed effects estimate") + ylab("copas estimate") + ggtitle("Fixed effects to copas")

meta.cont %>% ggplot(aes(x = abs(est.ranef.cont), y = abs(est.copas.cont))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Random effects estimate") + ylab("copas estimate") + ggtitle("Random effects to copas")

#Scatterplots of test-statistics:
metac.cont %>% ggplot(aes(x = abs(zval.fixef.cont), y = abs(est.copas.cont/se.est.copas.cont))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Fixed effects and copas test statistics") + ylab("copas statistic") + xlab("Fixed effects statistic")

metac.cont %>% ggplot(aes(x = abs(zval.ranef.cont), y = abs(est.copas.cont/se.est.copas.cont))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Random effects and copas test statistics") + ylab("copas statistic") + xlab("Random effects statistic")

########################################################################################################################
########################################################################################################################
#Regression lor comparison
########################################################################################################################
########################################################################################################################

#Summaries
metac.bin %>% mutate(xx = 1 - abs(logest.reg.ranef/log(est.fixef))) %>% filter(xx > 0) %>% ungroup() %>% 
  summarise(biased = sum(peter.test), larger.metac.effects = n())

metac.bin %>% mutate(xx = 1 - abs(logest.reg.ranef/log(est.fixef))) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(peter.test), smaller.metac.effects = n())

metac.bin %>% mutate(xx = 1 - abs(logest.reg.ranef/log(est.ranef))) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(peter.test), smaller.metac.effects.random = n())

#Histogramms of changes +- 100%
metac.bin %>% mutate(xx = 1 - abs(logest.reg.ranef/log(est.fixef))) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Regression")

metac.bin %>% mutate(xx = 1 - abs(logest.reg.ranef/log(est.fixef))) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Regression")

metac.bin %>% mutate(xx = 1 - abs(logest.reg.ranef/log(est.ranef))) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Regression")

metac.bin %>% mutate(xx = 1 - abs(logest.reg.ranef/log(est.ranef))) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(peter.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Regression")

metac.bin %>% mutate(xx = 1 - abs(logest.reg.ranef/log(est.fixef))) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Regression")

metac.bin %>% mutate(xx = 1 - abs(logest.reg.ranef/log(est.ranef))) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Regression")

#Scatterplots of metac-analysis and corrected estimates
metac.bin %>% ggplot(aes(x = abs(log(est.fixef)), y = abs(logest.reg.ranef))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Fixed effects estimate") + ylab("Regression estimate") + ggtitle("Fixed effects to Regression")

metac.bin %>% ggplot(aes(x = abs(log(est.ranef)), y = abs(logest.reg.ranef))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Random effects estimate") + ylab("Regression estimate") + ggtitle("Random effects to Regression")

#Scatterplots of test-statistics:
metac.bin %>% ggplot(aes(x = abs(zval.fixef), y = abs(logest.reg.ranef/se.logest.reg.ranef))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Fixed effects and reg.ranef test statistics") + ylab("reg.ranef statistic") + xlab("Fixed effects statistic")

metac.bin %>% ggplot(aes(x = abs(zval.ranef), y = abs(logest.reg.ranef/se.logest.reg.ranef))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Random effects and reg.ranef test statistics") + ylab("reg.ranef statistic") + xlab("Random effects statistic")



#CONTINUOUS:

#Summaries:
meta.cont %>% mutate(xx = 1 - abs(est.reg.ranef.cont/est.fixef.cont)) %>% filter(xx > 0) %>% ungroup() %>% 
  summarise(biased = sum(egger.test), larger.metac.effects = n())

meta.cont %>% mutate(xx = 1 - abs(est.reg.ranef.cont/est.fixef.cont)) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(egger.test), smaller.metac.effects = n())

meta.cont %>% mutate(xx = 1 - abs(est.reg.ranef.cont/est.ranef.cont)) %>% filter(xx < 0) %>% ungroup() %>% 
  summarise(biased = sum(egger.test), smaller.metac.effects.random = n())

#Histogramms of changes +- 100%
meta.cont %>% mutate(xx = 1 - abs(est.reg.ranef.cont/est.fixef.cont)) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Regression")

meta.cont %>% mutate(xx = 1 - abs(est.reg.ranef.cont/est.fixef.cont)) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Regression")

meta.cont %>% mutate(xx = 1 - abs(est.reg.ranef.cont/est.ranef.cont)) %>% filter(xx > 0) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Regression")

meta.cont %>% mutate(xx = 1 - abs(est.reg.ranef.cont/est.ranef.cont)) %>% filter(xx > -1) %>% 
  ggplot(aes(x = xx, fill = factor(egger.test))) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Regression")

meta.cont %>% mutate(xx = 1 - abs(est.reg.ranef.cont/est.fixef.cont)) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Fixed effects to Regression")

meta.cont %>% mutate(xx = 1 - abs(est.reg.ranef.cont/est.ranef.cont)) %>% 
  ggplot(aes(x = xx)) + geom_histogram() + xlab("Fraction of metac-analysis estimate") + ggtitle("Random effects to Regression")

#Scatterplots of metac-analysis and corrected estimates
meta.cont %>% ggplot(aes(x = abs(est.fixef.cont), y = abs(est.reg.ranef.cont))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Fixed effects estimate") + ylab("Regression estimate") + ggtitle("Fixed effects to Regression")

meta.cont %>% ggplot(aes(x = abs(est.ranef.cont), y = abs(est.reg.ranef.cont))) + geom_point() + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  xlab("Random effects estimate") + ylab("Regression estimate") + ggtitle("Random effects to Regression")

#Scatterplots of test-statistics:
metac.cont %>% ggplot(aes(x = abs(zval.fixef.cont), y = abs(est.reg.ranef.cont/se.est.reg.ranef.cont))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Fixed effects and reg.ranef test statistics") + ylab("reg.ranef statistic") + xlab("Fixed effects statistic")

metac.cont %>% ggplot(aes(x = abs(zval.ranef.cont), y = abs(est.reg.ranef.cont/se.est.reg.ranef.cont))) + geom_point() +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm") + 
  ggtitle("Random effects and reg.ranef test statistics") + ylab("reg.ranef statistic") + xlab("Random effects statistic")


# metac.bin %>% mutate(effect.reduction = abs((log(est.fixef) - logest.trimfill.fixef))) %>% ungroup() %>% 
#   summarise(max(effect.reduction), min(effect.reduction), 
#             mean(effect.reduction), median(effect.reduction))
