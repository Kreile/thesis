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

#Applying test and adjustment criteria:
metac.bin <- meta.bin %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac.cont <- meta.cont %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac <- meta %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)



tmp <- data.ext2 %>% filter(outcome.type == "bin") %>% group_by(meta.id) %>% 
	mutate(n = n()) %>% filter(n >= 10) %>% filter(!all(events1 == 0) & !all(events2 == 0)) %>% 
	filter(meta.id != 94519)

# tmp <- data.ext2 %>% filter(meta.id == 378)

meta.id.vector <- unique(tmp$meta.id) #Optimization of likelihood fails..

#[-c(which(unique(tmp$meta.id) == 8765), which(unique(tmp$meta.id) == 17300),
																				 # which(unique(tmp$meta.id) == 17305), which(unique(tmp$meta.id) == 17310),
																				 # which(unique(tmp$meta.id) == 17316), which(unique(tmp$meta.id) == 17320),
																				 # which(unique(tmp$meta.id) == 17324), which(unique(tmp$meta.id) == 17341),
																				 # which(unique(tmp$meta.id) == 17347), which(unique(tmp$meta.id) == 17361),
																				 # which(unique(tmp$meta.id) == 17364), which(unique(tmp$meta.id) == 60945),
																				 # which(unique(tmp$meta.id) == 23813), which(unique(tmp$meta.id) == 60946),
																				 # which(unique(tmp$meta.id) == 60947), which(unique(tmp$meta.id) == 76838),
																				 # which(unique(tmp$meta.id) == 76846), which(unique(tmp$meta.id) == 76888),
																				 # which(unique(tmp$meta.id) == 79353), which(unique(tmp$meta.id) == 76888))]
#8765 hat keine events, 17300 bis 17364 haben keine Angaben (nur Odds ratios und standard errors, wären umzurechnen, alles vom gleichen Review)
# meta.id.vector <- meta.id.vector[which(unique(tmp$meta.id) == 94517):which(meta.id.vector == 126898)]
meta.analyses <- list()
meta.analyses.limitmeta <- list()
meta.analyses.copas <- list()
meta.analyses.trimfill <- list()
meta.tests.harbord <- list()
meta.tests.peter <- list()
counter <- 0
x <- c()

for(u in meta.id.vector){
	counter <- counter + 1
	print(c(u, counter))
	meta.analyses[[counter]] <- metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, 
																			studlab = study.name, sm = "RR",
																			method = "Inverse", data = tmp[tmp$meta.id == u,])
	# meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
	# meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
	# meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
	# meta.tests.harbord[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "score")
	# meta.tests.peter[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "peters")
	# meta.tests.harbord[[counter]] <- is.na(metabias(meta.analyses[[counter]], method.bias = "linreg")$p.value)
	x <- c(x, metabias(meta.analyses[[counter]], method.bias = "linreg")$p.value)
}

meta.analyses[[counter]] <- metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, 
																		studlab = study.name, sm = "RR",
																		method = "Inverse", data = tmp[tmp$meta.id == 80,])
meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
meta.tests.harbord[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "score")
meta.tests.peter[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "peters")



meta.bin.complete2 = function(data, sig.level, min.study.number, sm1) {
	meta <- data.ext2 %>% filter(outcome.type == "bin") %>% 
		group_by(meta.id) %>% 
		mutate(n = n()) %>% filter(n >= 10) %>% 
		filter(!all(events1 == 0) & !all(events2 == 0)) %>% 
		filter(meta.id != 94519) %>% #(Copas likelihood optimization issues there)
		summarize(doi = unique(doi),
							outcome.type = "bin",
							outcome.mesaure.meta = "log Risk Ratio",
							file.nr = unique(file.nr),
							comparison.nr = unique(comparison.nr),
							comparison.name = unique(comparison.name),
							outcome.name = unique(outcome.name),
							subgroup.name = unique(subgroup.name),
							outcome.measure.new = unique(outcome.measure.new),
							n = n(), #Number of studies in meta-analysis
							mean.samplesize = mean(total1 + total2, na.rm = T),
							total.samplesize = sum(total1 + total2),
							var.samplesize = var(total1 + total2, na.rm = T),
							total.events = sum(events1 + events2),
							mean.events = mean(events1 + events2),
							mean.publication.year = mean(study.year, na.rm = TRUE),
							first.publication.year = min(study.year, na.rm = T),
							n.sig.single = sum(sig.single, na.rm = T),
							NA.sig.single = sum(is.na(sig.single)),
							se.min = min(se, na.rm = T),
							se.max = max(se, na.rm = T),
							
							#Meta-analysis:
							est.fixef = #Pooled odds ratio estimate of fixed effects meta-analysis model
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1,
														method = "Inverse")$TE.fixed,
							se.est.fixef = #Pooled odds ratio estimate of fixed effects meta-analysis model
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1,
														method = "Inverse")$seTE.fixed,
							zval.fixef = #P-value
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1,
												method = "Inverse")$zval.fixed,
							pval.fixef = #P-value
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1,
												method = "Inverse")$pval.fixed,
							# sig.fixef = ifelse(pval.fixef > sig.level, 0, 1), #If significant
							se.est.ranef = #Pooled odds ratio estimate of random effects meta-analysis model
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1)$seTE.random,
							est.ranef = #Pooled odds ratio estimate of random effects meta-analysis model
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1)$TE.random,
							zval.ranef = 
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1)$zval.random,
							pval.ranef = 
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1)$pval.random,
							# sig.ranef = ifelse(pval.ranef > sig.level, 0, 1),
							Q = #Q between studies heterogeneity statistic
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1)$Q,
							pval.Q = #P-value for between study heterogeneity statistic = 0
								metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1)$pval.Q,
							# sig.Q = ifelse(pval.Q < sig.level, 1, 0),
							I2  = max(0, (Q - n + 1)/Q), #Proportion of variance of study estimates that is due to real variance
							
							#Publication bias correction:
							est.trimfill.fixef = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
								trimfill(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1,
																	method = "Inverse"))$TE.fixed,
							est.trimfill.ranef = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
								trimfill(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1,
																	method = "Inverse"))$TE.random,
							se.est.trimfill.ranef = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1,
																 method = "Inverse"))$seTE.random,
							se.est.trimfill.fixef = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1,
																 method = "Inverse"))$seTE.fixed,
							zval.trimfill.fixef = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1,
																 method = "Inverse"))$zval.fixed,
							missing.trim = #Estimated number of missing studies in trimfill model
								trimfill(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1,
																 method = "Inverse"))$k0 / n(),
							pval.trimfill = 2*(1-pnorm(est.trimfill.ranef/se.est.trimfill.ranef)),
							sig.trimfill = ifelse(pval.trimfill > sig.level, 0, 1),
							est.reg.ranef = #Adjusted log odds ratio estimate of regression model
								limitmeta(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1))$TE.adjust,
							se.est.reg.ranef = #Standard error of adjusted log odds ratio of regression model
								limitmeta(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1))$seTE.adjust,
							pval.reg.ranef = #Pvalue of adjusted log odds ratio estimate of regression model
								limitmeta(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1))$pval.adjust,
							Q.resid = #Residual between study heterogeneity (after adjustment with regression model)
								limitmeta(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1))$Q.resid,
							G.squared = 
								limitmeta(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1))$G.squared,
							# sig.reg.ranef = ifelse(pval.reg.ranef > sig.level, 0, 1),
							est.copas = #Adjusted log odds ratio estimate of copas model
								auto.copas(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1), sig.level = sig.level)[1],
							se.est.copas = #Standard error of adjusted log odds ratio of copas model
								auto.copas(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1), sig.level = sig.level)[2],
							missing.copas = #Estimated number of missing studies in copas model
								auto.copas(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm= sm1), sig.level = sig.level)[3]/n(),
							pval.copas = 2*(1 - pnorm(abs(est.copas/se.est.copas))), #Pvalue of adjusted log odds ratio estimate of regression model
							# sig.copas = ifelse(pval.copas > sig.level, 0, 1),
							
							#Significance change
							# sig.change.ranef.reg = sig.change(sig.before = sig.ranef, sig.after =  sig.reg.ranef), #Function to register if significance of estimate has changed after correction
							# sig.change.ranef.copas = sig.change(sig.before = sig.ranef, sig.after =  sig.copas),
							# sig.change.ranef.trimfill = sig.change(sig.before = sig.ranef, sig.after =  sig.trimfill),
							# sig.change.fixef.reg = sig.change(sig.before = sig.fixef, sig.after =  sig.reg.ranef),
							# sig.change.fixef.copas = sig.change(sig.before = sig.fixef, sig.after =  sig.copas),
							# sig.change.fixef.trimfill = sig.change(sig.before = sig.fixef, sig.after =  sig.trimfill),
							# 
							#Publication bias tests:
							pval.egger = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																				method = "linreg")$p.val, #P-value (of eggers test)
							egger.test = ifelse(pval.egger < sig.level, 1, 0),
							pval.schwarzer = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																						method = "count")$p.val, #Significance of p-value
							schwarzer.test = ifelse(pval.schwarzer < sig.level, 1, 0),
							pval.peter = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																				method = "peter")$p.val,
							peter.test = ifelse(pval.peter < sig.level, 1, 0),
							pval.harbord = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																					method = "score")$p.val,
							harbord.test = ifelse(pval.harbord < sig.level, 1, 0),
							pval.rucker = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = "ASD"), 
																				 method = "mm")$p.val,
							rucker.test = ifelse(pval.rucker < sig.level, 1, 0),
							stat.egger = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																				method = "linreg")$statistic, #Test statistic (of eggers test)
							stat.schwarzer = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																						method = "count")$statistic,
							stat.peter = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																				method = "peter")$statistic,
							stat.harbord = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																					method = "score")$statistic,
							stat.rucker = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = "ASD"), 
																				 method = "mm")$statistic)
	return(meta)
}

meta.bin.complete2(data.ext[1:10000,], sig.level = 0.1, min.study.number = 10, sm1 = "RR")

tmp <- data.ext2 %>% filter(outcome.type == "cont") %>% group_by(meta.id) %>% 
	mutate(n = n()) %>% filter(n >= 10) %>% filter(!all(mean1 == 0) & !all(mean2 == 0)) %>% 
	filter(!all(sd1 == 0) | !all(sd2 == 0)) %>% filter(meta.id < 42716 | meta.id > 42725) #single patient data..

meta.id.vector <- unique(tmp$meta.id)[which(unique(tmp$meta.id)==42715):length(unique(tmp$meta.id))]
meta.analyses <- list()
meta.analyses.limitmeta <- list()
meta.analyses.copas <- list()
meta.analyses.trimfill <- list()
meta.tests.linreg <- list()
meta.tests.begg <- list()
meta.tests.mm <- list()
counter <- 0

for(u in meta.id.vector){
	counter <- counter + 1
	print(c(u, counter))
	meta.analyses[[counter]] <- metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, 
																			 mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name,
																			 data = tmp[tmp$meta.id == u,])
	meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
	meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
	meta.tests.linreg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "linreg")
	meta.tests.begg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "rank")
	meta.tests.mm[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "mm")
}


# tp <- data.ext2 %>% filter(meta.id == 3545)
# nomean1.ind <- which(tp$outcome.type == "cont" & is.na(tp$mean1))
# nomean2.ind <- which(tp$outcome.type == "cont" & is.na(tp$mean2))
# noeffects.ind <- which(tp$outcome.type == "cont" & !is.na(tp$effect))
# nomeans.ind <- match(nomean1.ind, nomean2.ind)
# nomeans.ind <- match(noeffects.ind, nomeans.ind)
# nomeans.ind <- nomeans.ind[which(!is.na(nomeans.ind))]
# 
# tp[nomeans.ind,"mean1"] <- tp[nomeans.ind,"effect"]
# tp[nomeans.ind,"mean2"] <- 0



meta.cont.complete2 = function(data, sig.level, min.study.number) {
	meta <- data %>% filter(outcome.type == "cont") %>% group_by(meta.id) %>% 
		mutate(n = n()) %>% filter(n >= 10) %>% filter(!all(mean1 == 0) & !all(mean2 == 0)) %>% 
		filter(!all(sd1 == 0) | !all(sd2 == 0)) %>% filter(meta.id < 42716 | meta.id > 42725) %>% #single patient data..
		summarize(doi = unique(doi),
							outcome.type =="cont",
							outcome.measure.meta = "Std. Mean Difference",
							file.nr = unique(file.nr),
							comparison.nr = unique(comparison.nr),
							comparison.name = unique(comparison.name),
							outcome.name = unique(outcome.name),
							subgroup.name = unique(subgroup.name),
							outcome.measure.new = unique(outcome.measure.new),
							mean.samplesize = mean(total1 + total2, na.rm = T),
							total.samplesize = sum(total1 + total2),
							var.samplesize = var(total1 + total2, na.rm = T),
							total.events = sum(events1 + events2),
							mean.events = mean(events1 + events2),
							mean.publication.year = mean(study.year, na.rm = T),
							first.publication.year = min(study.year, na.rm = T),
							n.sig.single = sum(sig.single, na.rm = T),
							NA.sig.single = sum(is.na(sig.single)),
							se.min = min(se, na.rm = T),
							se.max = max(se, na.rm = T),
							
							#Meta-analysis:
							est.fixef = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$TE.fixed,
							se.est.fixef = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$seTE.fixed,
							zval.fixef = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$zval.fixed,
							pval.fixef = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$pval.fixed,
							# sig.fixef = ifelse(pval.fixef > sig.level, 0, 1),
							est.ranef = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$TE.random,
							se.est.ranef = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$seTE.random,
							zval.ranef = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$zval.random,
							pval.ranef = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$pval.random,
							# sig.ranef = ifelse(pval.ranef > sig.level, 0, 1),
							Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$Q,
							pval.Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$pval.Q,
							# sig.Q = ifelse(pval.Q < sig.level, 1, 0),
							I2  = max(0, (Q - n + 1)/Q),
							
							#Pubbias correction:
							est.trimfill.fixef = #Trimfill
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$TE.fixed,
							est.trimfill.ranef = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$TE.random,
							se.est.trimfill.ranef = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$seTE.random,
							se.est.trimfill.fixef = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$seTE.fixed,
							zval.trimfill.fixef = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$zval.fixed,
							pval.trimfill = 2*(1-pnorm(abs(est.trimfill.ranef/se.est.trimfill.ranef))),
							# sig.trimfill = ifelse(pval.trimfill > sig.level, 0, 1),
							missing.trim = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$k0 / n(),
							est.reg.ranef = #Regression
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$TE.adjust,
							se.est.reg.ranef = 
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$seTE.adjust,
							pval.reg.ranef = 
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$pval.adjust,
							Q.resid = 
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$Q.resid,
							G.squared = 
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$G.squared,
							# sig.reg.ranef = ifelse(pval.reg.ranef > sig.level, 0, 1),
							est.copas = #Copas
								auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), sig.level = sig.level)[1],
							se.est.copas = 
								auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), sig.level = sig.level)[2],
							missing.copas = 
								auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), sig.level = sig.level)[3]/n(),
							pval.copas = 2*(1 - pnorm(abs(est.copas/se.est.copas))),
							# sig.copas = ifelse(pval.copas > sig.level, 0, 1),
							
							#Significance change
							# sig.change.ranef.reg = sig.change(sig.before = sig.ranef, sig.after =  sig.reg.ranef),
							# sig.change.ranef.copas = sig.change(sig.before = sig.ranef, sig.after =  sig.copas),
							# sig.change.ranef.trimfill = sig.change(sig.before = sig.ranef, sig.after =  sig.trimfill),
							# sig.change.fixef.reg = sig.change(sig.before = sig.fixef, sig.after =  sig.reg.ranef),
							# sig.change.fixef.copas = sig.change(sig.before = sig.fixef, sig.after =  sig.copas),
							# sig.change.fixef.trimfill = sig.change(sig.before = sig.fixef, sig.after =  sig.trimfill),
							# 
							#Pubbias tests
							pval.egger = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																				 method = "linreg")$p.val,
							egger.test = ifelse(pval.egger < sig.level, 1, 0),
							pval.thomson = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																					 method = "mm")$p.val,
							thomson.test = ifelse(pval.thomson < sig.level, 1, 0),
							pval.begg = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																				method = "rank")$p.val,
							begg.test = ifelse(pval.begg < sig.level, 1, 0),
							stat.egger = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																				 method = "linreg")$statistic,
							stat.thomson = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																					 method = "mm")$statistic,
							stat.begg = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																				method = "rank")$statistic) 
	return(meta)
}



meta.surv.complete2 = function(data, sig.level, min.study.number) {
	meta <- data %>% filter(outcome.type == "surv") %>% group_by(meta.id) %>% 
		mutate(n = n()) %>% filter(n >= min.study.number) %>% 
		summarize(doi = unique(doi),
							outcome.type = "surv",
							outcome.measure.meta = "log Hazard Ratio",
							file.nr = unique(file.nr),
							comparison.nr = unique(comparison.nr),
							comparison.name = unique(comparison.name),
							outcome.name = unique(outcome.name),
							subgroup.name = unique(subgroup.name),
							outcome.measure.new = unique(outcome.measure.new),
							mean.samplesize = mean(total1 + total2, na.rm = T),
							total.samplesize = sum(total1 + total2),
							var.samplesize = var(total1 + total2, na.rm = T),
							total.events = sum(events1 + events2),
							mean.events = mean(events1 + events2),
							mean.publication.year = mean(study.year, na.rm = T),
							first.publication.year = min(study.year, na.rm = T),
							n.sig.single = sum(sig.single, na.rm = T),
							NA.sig.single = sum(is.na(sig.single)),
							se.min = min(se, na.rm = T),
							se.max = max(se, na.rm = T),
							
							#Meta-analysis:
							est.fixef = #Pooled odds ratio estimate of fixed effects meta-analysis model
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$TE.fixed,
							se.est.fixef = #Pooled odds ratio estimate of fixed effects meta-analysis model
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$seTE.fixed,
							zval.fixef = #P-value
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$zval.fixed,
							pval.fixef = #P-value
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$pval.fixed,
							# sig.fixef = ifelse(pval.fixef > sig.level, 0, 1), #If significant
							se.est.ranef = #Pooled odds ratio estimate of random effects meta-analysis model
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$seTE.random,
							est.ranef = #Pooled odds ratio estimate of random effects meta-analysis model
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$TE.random,
							zval.ranef = 
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$zval.random,
							pval.ranef = 
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$pval.random,
							# sig.ranef = ifelse(pval.ranef > sig.level, 0, 1),
							Q = #Q between studies heterogeneity statistic
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$Q,
							pval.Q = #P-value for between study heterogeneity statistic = 0
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$pval.Q,
							# sig.Q = ifelse(pval.Q < sig.level, 1, 0),
							I2  = max(0, (Q - n + 1)/Q), #Proportion of variance of study estimates that is due to real variance
							
							#Publication bias correction:
							est.trimfill.fixef = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$TE.fixed,
							est.trimfill.ranef = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$TE.random,
							se.est.trimfill.ranef = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$seTE.random,
							se.est.trimfill.fixef = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$seTE.fixed,
							zval.trimfill.fixef = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$zval.fixed,
							missing.trim = #Estimated number of missing studies in trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$k0 / n(),
							pval.trimfill = 2*(1-pnorm(est.trimfill.ranef/se.est.trimfill.ranef)),
							# sig.trimfill = ifelse(pval.trimfill > sig.level, 0, 1),
							est.reg.ranef = #Adjusted log odds ratio estimate of regression model
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$TE.adjust,
							se.est.reg.ranef = #Standard error of adjusted log odds ratio of regression model
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$seTE.adjust,
							pval.reg.ranef = #Pvalue of adjusted log odds ratio estimate of regression model
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$pval.adjust,
							Q.resid = #Residual between study heterogeneity (after adjustment with regression model)
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$Q.resid,
							G.squared = 
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$G.squared,
							# sig.reg.ranef = ifelse(pval.reg.ranef > sig.level, 0, 1),
							est.copas = #Adjusted log odds ratio estimate of copas model
								auto.copas(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), sig.level = sig.level)[1],
							se.est.copas = #Standard error of adjusted log odds ratio of copas model
								auto.copas(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), sig.level = sig.level)[2],
							missing.copas = #Estimated number of missing studies in copas model
								auto.copas(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), sig.level = sig.level)[3]/n(),
							pval.copas = 2*(1 - pnorm(abs(est.copas/se.est.copas))), #Pvalue of adjusted log odds ratio estimate of regression model
							# sig.copas = ifelse(pval.copas > sig.level, 0, 1),
							
							#Significance change
							# sig.change.ranef.reg = sig.change(sig.before = sig.ranef, sig.after =  sig.reg.ranef), #Function to register if significance of estimate has changed after correction
							# sig.change.ranef.copas = sig.change(sig.before = sig.ranef, sig.after =  sig.copas),
							# sig.change.ranef.trimfill = sig.change(sig.before = sig.ranef, sig.after =  sig.trimfill),
							# sig.change.fixef.reg = sig.change(sig.before = sig.fixef, sig.after =  sig.reg.ranef),
							# sig.change.fixef.copas = sig.change(sig.before = sig.fixef, sig.after =  sig.copas),
							# sig.change.fixef.trimfill = sig.change(sig.before = sig.fixef, sig.after =  sig.trimfill),
							
							#Pubbias tests
							pval.egger = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																				 method = "linreg")$p.val,
							egger.test = ifelse(pval.egger < sig.level, 1, 0),
							pval.thomson = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																					 method = "mm")$p.val,
							thomson.test = ifelse(pval.thomson < sig.level, 1, 0),
							pval.begg = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																				method = "rank")$p.val,
							begg.test = ifelse(pval.begg < sig.level, 1, 0),
							stat.egger = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																				 method = "linreg")$statistic,
							stat.thomson = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																					 method = "mm")$statistic,
							stat.begg = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																				method = "rank")$statistic) 
	return(meta)
}

meta.surv.complete2(data = data.ext2[1:50000, ], sig.level = 0.1, min.study.number = 10)
meta.bin.complete2(data.ext[1:40000,], sig.level = 0.1, min.study.number = 10, sm1 = "RR")
meta.cont.complete2(data = data.ext2[1:50000, ], sig.level = 0.1, min.study.number = 10)


meta.surv <- meta.surv.complete2(data = data.ext2, sig.level = 0.1, min.study.number = 10)
meta.bin <- meta.bin.complete2(data = data.ext2, sig.level = 0.1, min.study.number = 10, sm1 = "RR")
meta.cont <- meta.cont.complete2(data = data.ext2, sig.level = 0.1, min.study.number = 10)

meta <- bind_rows(meta.bin, meta.cont, meta.surv)

pb.meta.merge <- function(meta.binary, meta.continous){
	meta <- bind_rows(meta.bin, meta.cont) %>%
		mutate(pval.thomson = ifelse(!is.na(pval.rucker.bin), pval.rucker.bin, pval.thomson.cont),
					 thomson.test = ifelse(!is.na(thomson.test), thomson.test, rucker.test),
					 
					 est.fixef = ifelse(!is.na(est.fixef.bin), est.fixef.bin, est.fixef.cont),
					 zval.fixef = ifelse(!is.na(zval.fixef.bin), zval.fixef.bin, zval.fixef.cont),
					 
					 est.ranef = ifelse(!is.na(est.ranef.bin), est.ranef.bin, est.ranef.cont),
					 zval.ranef = ifelse(!is.na(zval.ranef.bin), zval.ranef.bin, zval.ranef.cont),
					 
					 est.trimfill.fixef = ifelse(!is.na(est.trimfill.fixef.bin), exp(logest.trimfill.fixef.bin), est.trimfill.fixef.cont),
					 se.est.trimfill.fixef = ifelse(!is.na(se.logest.trimfill.fixef.bin), se.logest.trimfill.fixef.bin, se.est.trimfill.fixef.cont),
					 
					 est.reg.ranef = ifelse(!is.na(logest.reg.ranef.bin), exp(logest.reg.ranef.bin), est.reg.ranef.cont),
					 se.est.reg.ranef = ifelse(!is.na(se.logest.reg.ranef.bin), se.logest.reg.ranef.bin, se.est.reg.ranef.cont),
					 zval.reg.ranef = ifelse(!is.na(logest.reg.ranef.bin), logest.reg.ranef.bin/se.logest.reg.ranef.bin,
					 												est.reg.ranef.cont/se.est.reg.ranef.cont),
					 
					 est.copas = ifelse(!is.na(logest.copas.bin), exp(logest.copas.bin), est.copas.cont),
					 se.est.copas = ifelse(!is.na(se.logest.copas.bin), se.logest.copas.bin, se.est.copas.cont),
					 zval.copas = ifelse(!is.na(logest.copas.bin), logest.copas.bin/se.logest.copas.bin,
					 										est.copas.cont/se.est.copas.cont),
					 
					 sig.ranef = ifelse(!is.na(sig.ranef.bin), sig.ranef.bin, sig.ranef.cont),
					 sig.fixef = ifelse(!is.na(sig.fixef.bin), sig.fixef.bin, sig.fixef.cont),
					 pval.fixef = ifelse(!is.na(pval.fixef.bin), pval.fixef.bin, pval.fixef.cont),
					 pval.ranef = ifelse(!is.na(pval.ranef.bin), pval.ranef.bin, pval.ranef.cont),
					 pval.copas = ifelse(!is.na(pval.copas.bin), pval.copas.bin, pval.copas.cont),
					 pval.reg.ranef = ifelse(!is.na(pval.reg.ranef.bin), pval.reg.ranef.bin, pval.reg.ranef.cont),
					 pval.trimfill = ifelse(!is.na(pval.trimfill.bin), pval.trimfill.bin, pval.trimfill.cont),
					 missing.trim = ifelse(!is.na(missing.trim.bin), missing.trim.bin, missing.trim.cont),
					 missing.copas = ifelse(!is.na(missing.copas.bin), missing.copas.bin, missing.copas.cont),
					 sig.copas = ifelse(!is.na(sig.copas.bin), sig.copas.bin, sig.copas.cont),
					 sig.trimfill = ifelse(!is.na(sig.trimfill.bin), sig.trimfill.bin, sig.trimfill.cont),
					 sig.reg.ranef = ifelse(!is.na(sig.reg.ranef.bin), sig.reg.ranef.bin, sig.reg.ranef.cont),
					 Q = ifelse(!is.na(Q.bin), Q.bin, Q.cont),
					 pval.Q = ifelse(!is.na(pval.Q.bin), pval.Q.bin, pval.Q.cont),
					 sig.Q = ifelse(!is.na(sig.Q.bin), sig.Q.bin, sig.Q.cont))
	return(meta)
}
