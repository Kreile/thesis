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



tmp <- data.ext2 %>% filter(outcome.type == "surv") %>% group_by(meta.id) %>% 
	mutate(n = n()) %>% filter(n >= 10) 

# tmp <- data.ext2 %>% filter(meta.id == 378)

meta.id.vector <- unique(tmp$meta.id)
meta.analyses <- list()
meta.analyses.limitmeta <- list()
meta.analyses.copas <- list()
meta.analyses.trimfill <- list()
meta.tests.harbord <- list()
meta.tests.peter <- list()
meta.tests.begg <- list()
counter <- 0


for(u in meta.id.vector){
	counter <- counter + 1
	print(c(u, counter))
	meta.analyses[[counter]] <- metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR", data = tmp[tmp$meta.id == u,])
	meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
	meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
	meta.tests.harbord[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "linreg")
	meta.tests.peter[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "mm")
	meta.tests.begg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "rank")
}

meta.analyses[[counter]] <- metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, 
																		studlab = study.name, sm = "RR",
																		method = "Inverse", data = tmp[tmp$meta.id == 80,])
meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
meta.tests.harbord[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "score")
meta.tests.peter[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "peters")



meta.surv.complete = function(data, sig.level, min.study.number) {
	meta.surv <- data %>% filter(outcome.type == "surv") %>% group_by(meta.id) %>% 
		mutate(n = n()) %>% filter(n >= min.study.number) %>% 
		summarize(doi = unique(doi),
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
							n.sig.single.surv = sum(sig.single, na.rm = T),
							NA.sig.single = sum(is.na(sig.single)),
							se.min = min(se, na.rm = T),
							se.max = max(se, na.rm = T),
							
							#Meta-analysis:
							est.fixef.surv = #Pooled odds ratio estimate of fixed effects meta-analysis model
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$TE.fixed,
							se.est.fixef.surv = #Pooled odds ratio estimate of fixed effects meta-analysis model
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$seTE.fixed,
							zval.fixef.surv = #P-value
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$zval.fixed,
							pval.fixef.surv = #P-value
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$pval.fixed,
							# sig.fixef.surv = ifelse(pval.fixef.surv > sig.level, 0, 1), #If significant
							se.est.ranef.surv = #Pooled odds ratio estimate of random effects meta-analysis model
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$seTE.random,
							est.ranef.surv = #Pooled odds ratio estimate of random effects meta-analysis model
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$TE.random,
							zval.ranef.surv = 
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$zval.random,
							pval.ranef.surv = 
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$pval.random,
							# sig.ranef.surv = ifelse(pval.ranef.surv > sig.level, 0, 1),
							Q.surv = #Q between studies heterogeneity statistic
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$Q,
							pval.Q.surv = #P-value for between study heterogeneity statistic = 0
								metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR")$pval.Q,
							# sig.Q.surv = ifelse(pval.Q.surv < sig.level, 1, 0),
							I2  = max(0, (Q.surv - n + 1)/Q.surv), #Proportion of variance of study estimates that is due to real variance
							
							#Publication bias correction:
							logest.trimfill.fixef.surv = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$TE.fixed,
							logest.trimfill.ranef.surv = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$TE.random,
							se.logest.trimfill.ranef.surv = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$seTE.random,
							se.logest.trimfill.fixef.surv = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$seTE.fixed,
							zval.trimfill.fixef = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$zval.fixed,
							missing.trim.surv = #Estimated number of missing studies in trimfill model
								trimfill(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$k0 / n(),
							pval.trimfill.surv = 2*(1-pnorm(logest.trimfill.ranef.surv/se.logest.trimfill.ranef.surv)),
							# sig.trimfill.surv = ifelse(pval.trimfill.surv > sig.level, 0, 1),
							logest.reg.ranef.surv = #Adjusted log odds ratio estimate of regression model
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$TE.adjust,
							se.logest.reg.ranef.surv = #Standard error of adjusted log odds ratio of regression model
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$seTE.adjust,
							pval.reg.ranef.surv = #Pvalue of adjusted log odds ratio estimate of regression model
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$pval.adjust,
							Q.resid = #Residual between study heterogeneity (after adjustment with regression model)
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$Q.resid,
							G.squared = 
								limitmeta(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"))$G.squared,
							# sig.reg.ranef.surv = ifelse(pval.reg.ranef.surv > sig.level, 0, 1),
							logest.copas.surv = #Adjusted log odds ratio estimate of copas model
								auto.copas(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), sig.level = sig.level)[1],
							se.logest.copas.surv = #Standard error of adjusted log odds ratio of copas model
								auto.copas(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), sig.level = sig.level)[2],
							missing.copas.surv = #Estimated number of missing studies in copas model
								auto.copas(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), sig.level = sig.level)[3]/n(),
							pval.copas.surv = 2*(1 - pnorm(abs(logest.copas.surv/se.logest.copas.surv))), #Pvalue of adjusted log odds ratio estimate of regression model
							# sig.copas.surv = ifelse(pval.copas.surv > sig.level, 0, 1),
							
							#Significance change
							# sig.change.ranef.reg = sig.change(sig.before = sig.ranef.surv, sig.after =  sig.reg.ranef.surv), #Function to register if significance of estimate has changed after correction
							# sig.change.ranef.copas = sig.change(sig.before = sig.ranef.surv, sig.after =  sig.copas.surv),
							# sig.change.ranef.trimfill = sig.change(sig.before = sig.ranef.surv, sig.after =  sig.trimfill.surv),
							# sig.change.fixef.reg = sig.change(sig.before = sig.fixef.surv, sig.after =  sig.reg.ranef.surv),
							# sig.change.fixef.copas = sig.change(sig.before = sig.fixef.surv, sig.after =  sig.copas.surv),
							# sig.change.fixef.trimfill = sig.change(sig.before = sig.fixef.surv, sig.after =  sig.trimfill.surv),
							
							#Pubbias tests
							pval.egger.surv = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																				 method = "linreg")$p.val,
							egger.test = ifelse(pval.egger.surv < sig.level, 1, 0),
							pval.thomson.surv = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																					 method = "mm")$p.val,
							thomson.test = ifelse(pval.thomson.surv < sig.level, 1, 0),
							pval.begg.surv = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																				method = "rank")$p.val,
							begg.test = ifelse(pval.begg.surv < sig.level, 1, 0),
							stat.egger.surv = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																				 method = "linreg")$statistic,
							stat.thomson.surv = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																					 method = "mm")$statistic,
							stat.begg.surv = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																				method = "rank")$statistic) 
	return(meta.surv)
}

meta.surv.complete(data = data.ext2[1:50000, ], sig.level = 0.1, min.study.number = 10)
