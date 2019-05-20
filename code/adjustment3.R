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
		meta.ex <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "RR", data = dt, method = "Inverse")
	}
	return(meta.ex)
}

funnel.id <- function(meta.id.tolook){
	funnel(meta.id(meta.id.tolook))
}

trimfill.id <- function(meta.id.tolook){
	funnel(trimfill(meta.id(meta.id.tolook)))
}

limitmeta.id <- function(meta.id.tolook){
	funnel(limitmeta(meta.id(meta.id.tolook)))
}

copas.id <- function(meta.id.tolook){
	auto.copas(meta.id(meta.id.tolook), 0.1)
}

setting.id <- function(meta.id.tolook){
	return(meta %>% filter(meta.id == meta.id.tolook) %>% select(comparison.name, comparison.nr, outcome.name, sungroup.name))
}

est.id <- function(meta.id.tolook){
	return(meta %>% filter(meta.id == meta.id.tolook) %>% select(est.fixef, est.ranef, est.copas, est.trimfill.fixef, est.reg.ranef))
}

single.study.id <- function(meta.id.tolook){
	return(print(data.ext %>% filter(meta.id == meta.id.tolook) %>% select(study.name, effect, se, total1, total2, events1, events2, mean1, mean2) %>% 
							 	arrange(se), n = 200))
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
	meta.bin <- meta.bin.complete(data.ext, min.study.number = 10, sig.level = 0.1, sm1 = "RR")
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


effect.diff <- metac %>% mutate(
	est.fixef = ifelse(!is.na(est.fixef.bin), est.fixef - 1, est.fixef),
	est.ranef = ifelse(!is.na(est.fixef.bin), est.ranef - 1, est.ranef),
	est.trimfill.fixef = ifelse(!is.na(est.fixef.bin), est.trimfill.fixef - 1, est.trimfill.fixef),
	est.copas = ifelse(!is.na(est.fixef.bin), est.copas - 1, est.copas),
	est.reg.ranef = ifelse(!is.na(est.fixef.bin), est.reg.ranef - 1, est.reg.ranef),
	fixef.trimfill = ifelse(abs(est.fixef) > abs(est.trimfill.fixef), "Reduction", "Amplification"),
	fixef.copas = ifelse(abs(est.fixef) > abs(est.copas), "Reduction", "Amplification"),
	fixef.reg = ifelse(abs(est.fixef) > abs(est.reg.ranef), "Reduction", "Amplification"),
	fixef.trimfill.s = ifelse(sign(est.fixef) == sign(est.trimfill.fixef), "Unchanged", "Reversed"),
	fixef.copas.s = ifelse(sign(est.fixef) == sign(est.copas), "Unchanged", "Reversed"),
	fixef.reg.s = ifelse(sign(est.fixef) == sign(est.reg.ranef), "Unchanged", "Reversed")) 


trimfill.miss <- effect.diff %>% group_by(fixef.trimfill) %>% mutate(dif = abs(est.fixef - est.trimfill.fixef)) %>% 
	filter(dif > 2) %>% mutate(label = round(dif, 1), type = !is.na(est.fixef.bin)) %>% arrange(label) %>% select(label, fixef.trimfill.s, meta.id, type)

print(effect.diff %>% group_by(fixef.trimfill) %>% mutate(dif = abs(est.fixef - est.trimfill.fixef)) %>% 
	filter(dif > 2) %>% mutate(label = round(dif, 1), type = !is.na(est.fixef.bin)) %>% arrange(label) %>% 
	select(label, fixef.trimfill.s, meta.id, type, comparison.name, outcome.name), n = 200)

trimfill.label <- data.frame(label = paste("missing:", c(paste(trimfill.miss$label[trimfill.miss$fixef.trimfill == "Reduction"], collapse = ", "),
																												 paste(trimfill.miss$label[trimfill.miss$fixef.trimfill == "Amplification"], collapse = ", "))),
														fixef.trimfill = unique(trimfill.miss$fixef.trimfill), 
														fixef.trimfill.s  = "Reversed")


effect.diff %>% ggplot(aes(x = abs(est.fixef - est.trimfill.fixef), fill = factor(fixef.trimfill.s))) +
	geom_histogram() + facet_grid(factor(fixef.trimfill)~.) + xlim(0, 2)+ xlab("Effect difference") +
	geom_text(data = trimfill.label, mapping = aes(label = label, x = 1.25, y = 200), 
						color = "black") + scale_fill_discrete(name="Effect direction") + theme_bw()

print(effect.diff %>% group_by(fixef.trimfill) %>% mutate(dif = abs(est.fixef - est.trimfill.fixef)) %>% 
				filter(dif > 2) %>% mutate(label = round(dif, 1), type = !is.na(est.fixef.bin)) %>% arrange(label) %>% 
				select(label, fixef.trimfill.s, meta.id, type), n = 200)

#Odds ratio reduction by 3.9
funnel(trimfill(meta.id(165815)))


copas.miss <- effect.diff %>% filter(!is.na(est.copas)) %>% group_by(fixef.copas) %>% mutate(dif = abs(est.fixef - est.copas)) %>% 
	filter(dif > 2) %>% mutate(label = round(dif, 1), type = !is.na(est.fixef.bin)) %>% arrange(label) %>% select(label, fixef.copas.s, meta.id, type)

copas.label <- data.frame(label = paste("missing:", c(paste(copas.miss$label[copas.miss$fixef.copas == "Reduction"], collapse = ", "),
																												 paste(copas.miss$label[copas.miss$fixef.copas == "Amplification"], collapse = ", "))),
														 fixef.copas = unique(copas.miss$fixef.copas), 
														 fixef.copas.s  = "Reversed")


effect.diff %>% filter(!is.na(est.copas)) %>% ggplot(aes(x = abs(est.fixef - est.copas), fill = factor(fixef.copas.s))) +
	geom_histogram() + facet_grid(factor(fixef.copas)~.) + xlim(0, 2)+ xlab("Effect difference") +
	geom_text(data = copas.label, mapping = aes(label = label, x = 1.25, y = 100), 
						color = "black") + scale_fill_discrete(name="Effect direction") + theme_bw()

funnel(meta.id(56200))

print(effect.diff %>% group_by(fixef.copas) %>% mutate(dif = abs(est.fixef - est.copas)) %>% 
				filter(dif > 2) %>% filter(fixef.trimfill == "Amplification") %>% mutate(label = round(dif, 1), type = !is.na(est.fixef.bin)) %>% arrange(label) %>% 
				select(label,  meta.id, type, comparison.name, outcome.name), n = 200)

reg.ranef.label <- data.frame(label = paste("missing:", c(paste(reg.ranef.miss$label[reg.ranef.miss$fixef.reg == "Reduction"], collapse = ", "),
																											paste(reg.ranef.miss$label[reg.ranef.miss$fixef.reg == "Amplification"], collapse = ", "))),
													fixef.reg = unique(reg.ranef.miss$fixef.reg), 
													fixef.reg.s  = "Reversed")


effect.diff %>% ggplot(aes(x = abs(est.fixef - est.reg.ranef), fill = factor(fixef.reg.s))) +
	geom_histogram() + facet_grid(factor(fixef.reg)~.) + xlim(0, 3)+ xlab("Effect difference") +
	geom_text(data = reg.ranef.label, mapping = aes(label = label, x = 1.25, y = 100), 
						color = "black") + scale_fill_discrete(name="Effect direction") + theme_bw()

metac %>% filter(meta.id == 101436) %>% select(comparison.name, outcome.name, sungroup.name, est.fixef, est.copas, 
																							 est.trimfill.fixef, est.reg.ranef, est.copas)
par(mfrow = c(1,3))
funnel(meta.id(101436))
funnel(trimfill(meta.id(101436)))
funnel(limitmeta(meta.id(101436)))

print(effect.diff %>% group_by(fixef.reg) %>% mutate(dif = abs(est.fixef - est.reg.ranef)) %>% 
				filter(dif > 3) %>% mutate(label = round(dif, 1), type = !is.na(est.fixef.bin)) %>% 
				arrange(label) %>% select(label, type, meta.id, comparison.name, outcome.name), n = 200)


