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


# data = pb.readData(path = PATH_DATA, file = FILE)
# tmp = pb.clean(data)
# data = tmp[[1]]
# aliases = tmp[[2]]


meta.id <- function(meta.id.tolook){
	dt <- data.ext %>% filter(meta.id == meta.id.tolook)
	if(any(is.na(dt$events1)) | all(dt$events1 == 0)){
		meta.ex <- metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, data = dt)
	} else{
		meta.ex <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR", data = dt)
	}
	return(meta.ex)
}

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


file.cont <- "mly.cont.RData"
if (file.exists(file.path(PATH_RESULTS, file.cont))) {
	load(file.path(PATH_RESULTS, file.cont))
} else {
	data.cont <- mly.cont(data.ext, 0.05, min.study.number = 2)
	save(data.cont, file =  file.path(PATH_RESULTS, file.cont))
}

file.bin <- "mly.bin.RData"
if (file.exists(file.path(PATH_RESULTS, file.bin))) {
	load(file.path(PATH_RESULTS, file.bin))
} else {
	data.bin <- mly.bin(data.ext, 0.05, min.study.number = 2)
	save(data.bin, file =  file.path(PATH_RESULTS, file.bin))
}


require(biostatUZH)
require(tidyverse)
require(meta)
require(metasens)
require(gridExtra)


#Meta filtering: 
metac.bin <- meta.bin %>% filter(n.sig.type.bin > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac.cont <- meta.cont %>% filter(n.sig.type.cont > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac <- meta %>% filter(n.sig.type > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)

# metac %>% filter(meta.id == 30578) %>% select(est.trimfill.fixef, est.fixef, est.copas, est.reg.ranef)
# 
# print(meta.bin.complete(tmp, 10, sig.level = 0.1, sm = "RR") %>% select(doi, se.min, se.max, total.events,
# 																							logest.trimfill.fixef.bin, est.fixef.bin, logest.copas.bin, logest.reg.ranef.bin))
# 
# print(meta.bin.complete(tmp, 10, sig.level = 0.1, sm = "RR") %>% select(doi, se.min, se.max, total.events, logest.trimfill.fixef.bin,
# 																							est.fixef.bin, logest.copas.bin, logest.reg.ranef.bin))


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



#Scatterplots of test-log statistics:
metac %>% ggplot(aes(x = log(abs(zval.fixef), y = log(abs(zval.reg.ranef))) + geom_point() +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	ggtitle("Fixed effects and trimfill z statistics") + ylab("Trimfill log statistic") + xlab("Fixed effects log statistic")

metac %>% ggplot(aes(x = log(abs(zval.ranef), y = log(abs(zval.reg.ranef))) + geom_point() +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	ggtitle("Random effects and trimfill z statistics") + ylab("Trimfill log statistic") + xlab("Random effects log statistic")
	
	#Scatterplots of test-log statistics:
	metac %>% ggplot(aes(x = log(abs(zval.fixef), y = log(abs(zval.trimfill.fixef))) + geom_point() +
											 	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
											 	ggtitle("Fixed effects and trimfill z statistics") + ylab("Trimfill log statistic") + xlab("Fixed effects log statistic")
											 
											 metac %>% ggplot(aes(x = log(abs(zval.ranef), y = log(abs(zval.trimfill.fixef))) + geom_point() +
											 										 	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
											 										 	ggtitle("Random effects and trimfill z statistics") + ylab("Trimfill log statistic") + xlab("Random effects log statistic")
											 										 
											 										 #Scatterplots of test-log statistics:
											 										 metac %>% ggplot(aes(x = log(abs(zval.fixef), y = log(abs(zval.copas))) + geom_point() +
											 										 										 	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
											 										 										 	ggtitle("Fixed effects and trimfill z statistics") + ylab("Trimfill log statistic") + xlab("Fixed effects log statistic")
											 										 										 
											 										 										 metac %>% ggplot(aes(x = log(abs(zval.ranef), y = log(abs(zval.copas))) + geom_point() +
											 										 										 										 	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
											 										 										 										 	ggtitle("Random effects and trimfill z statistics") + ylab("Trimfill log statistic") + xlab("Random effects log statistic")
