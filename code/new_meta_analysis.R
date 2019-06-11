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

# load(file.path(PATH_RESULTS, "meta_analyses_summary_bin_cont_surv.RData"))
# load(file.path(PATH_RESULTS, "meta_analyses_summary_complete.RData"))
# load(file.path(PATH_RESULTS, "data_used_for_summary.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_bin.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_cont.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_surv.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_cor.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_zscore.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_cum.RData"))

settings.meta(hakn=TRUE, method.tau="PM", method = "Inverse")

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#META-ANALYSIS CALCULATION:
#Binary outcomes:
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

tmp <- data.ext2 %>% filter(outcome.type == "bin") %>% group_by(meta.id) %>% 
	mutate(n = n()) %>% filter(n >= 10) %>% filter(!all(events1 == 0) & !all(events2 == 0)) 

meta.id.vector <- unique(tmp$meta.id)
meta.id.vector <- meta.id.vector[-c(which(meta.id.vector == 94519),
                                   which(meta.id.vector == 62301))] #Copas likelihood conv. problems..
meta.analyses <- list()
meta.analyses.asd <- list()
meta.analyses.limitmeta <- list()
meta.analyses.copas <- list()
meta.analyses.trimfill <- list()
meta.tests.harbord <- list()
meta.tests.peter <- list()
meta.tests.rucker <- list()
meta.tests.schwarzer <- list()
counter <- 0


# for(u in meta.id.vector){
# 	counter <- counter + 1
# 	print(c(u, counter))
# 	meta.analyses[[counter]] <- metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, 
# 																			studlab = study.name, sm = "RR",
# 																			method = "Inverse", data = tmp[tmp$meta.id == u,])
# 	meta.analyses.asd[[counter]] <- metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, 
# 																			studlab = study.name, sm = "ASD",
# 																			method = "Inverse", data = tmp[tmp$meta.id == u,])
# 	
# 	meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
# 	meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
# 	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
# 	
# 	meta.tests.harbord[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "score", k.min = 2)
# 	meta.tests.peter[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "peters", k.min = 2)
# 	meta.tests.schwarzer[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "count", k.min = 2)
# 	meta.tests.rucker[[counter]] <- metabias(meta.analyses.asd[[counter]], method.bias = "mm", k.min = 2)
# }
# 
# for(u in meta.id.vector){
# 	if(length(meta.analyses.copas[[counter]]) < 3){
# 		meta.analyses.copas[[counter]] <- c(NA, NA, NA)
# 	}
# 	
# }

meta.analyses <- bin.meta.list[[3]]
meta.analyses.limitmeta <- bin.meta.list[[5]]
meta.analyses.copas <- bin.meta.list[[6]]
meta.analyses.trimfill <- bin.meta.list[[7]]
meta.tests.harbord <- bin.meta.list[[8]]
meta.tests.peter <- bin.meta.list[[9]]
meta.tests.rucker <- bin.meta.list[[10]]
meta.tests.schwarzer <- bin.meta.list[[11]]

bin.meta.list <- list(tmp, meta.id.vector, meta.analyses,
											meta.analyses.asd,
											meta.analyses.limitmeta,
											meta.analyses.copas,
											meta.analyses.trimfill,
											meta.tests.harbord,
											meta.tests.peter,
											meta.tests.rucker,
											meta.tests.schwarzer)

#Extract from lists:
bin.results <- data.frame(
	meta.id = meta.id.vector,
	meta.es.measure = rep(times = length(meta.id.vector), "log risk ratio"),
	
	n.sig.single2 = unlist(lapply(meta.analyses, FUN = function(meta.analysis){length(which(meta.analysis$pval < 0.05))})),
	
	#Meta-analysis:
	est.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.fixed})),
	est.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.random})),
	
	zval.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$zval.fixed})),
	zval.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$zval.random})),
	
	pval.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$pval.fixed})),
	pval.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$pval.random})),
	
	se.est.fixef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.fixed})),
	se.est.ranef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.random})),
	
	Q = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$Q})),
	pval.Q = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$pval.Q})),
	#I2 = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$I2.w})),
	tau = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$tau})),
	#se.tau = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$se.tau})),
	sparse = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$sparse})),
	k = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$k})),
	
	#Adjustment:
	method.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$method.adjust})),
	est.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
	se.est.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
	zval.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$zval.adjust})),
	pval.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
	alpha.r = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
	beta.r = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
	Q.small = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$Q.small})),
	Q.resid = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$Q.resid})),
	G.squared = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$G.squared})),
	
	est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
	se.est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})),
	missing.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[3]})),
	
	#Tests;
	pval.harbord = unlist(lapply(meta.tests.harbord, FUN = function(meta.test){meta.test$p.value})),
	stat.harbord = unlist(lapply(meta.tests.harbord, FUN = function(meta.test){meta.test$statistic})),
	pval.peter = unlist(lapply(meta.tests.peter, FUN = function(meta.test){meta.test$p.value})),
	stat.peter = unlist(lapply(meta.tests.peter, FUN = function(meta.test){meta.test$statistic})),
	pval.rucker = unlist(lapply(meta.tests.rucker, FUN = function(meta.test){meta.test$p.value})),
	stat.rucker = unlist(lapply(meta.tests.rucker, FUN = function(meta.test){meta.test$statistic})),
	pval.schwarzer = unlist(lapply(meta.tests.schwarzer, FUN = function(meta.test){meta.test$p.value})),
	stat.schwarzer = unlist(lapply(meta.tests.schwarzer, FUN = function(meta.test){meta.test$statistic}))
	)



##########################################################################################
##########################################################################################
#Continuous outcomes:
##########################################################################################
##########################################################################################
tmp <- data.ext2 %>% filter(outcome.type == "cont") %>% group_by(meta.id) %>% 
	mutate(n = n()) %>% filter(n >= 10) %>% filter(!all(mean1 == 0) & !all(mean2 == 0)) %>% 
	filter(!all(sd1 == 0) | !all(sd2 == 0)) 

meta.id.vector <- unique(tmp$meta.id)
meta.id.vector <- meta.id.vector[-c(which(meta.id.vector > 42716 & meta.id.vector < 42725))] #single patient data..

# meta.analyses <- list()
# meta.analyses.limitmeta <- list()
# meta.analyses.copas <- list()
# meta.analyses.trimfill <- list()
# meta.tests.linreg <- list()
# meta.tests.begg <- list()
# meta.tests.mm <- list()
# counter <- 0
# 
# for(u in meta.id.vector){
# 	counter <- counter + 1
# 	print(c(u, counter))
# 	meta.analyses[[counter]] <- metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, 
# 																			 mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name,
# 																			 data = tmp[tmp$meta.id == u,])
# 	meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
# 	meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
# 	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
# 	
# 	meta.tests.linreg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "linreg")
# 	meta.tests.begg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "rank")
# 	meta.tests.mm[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "mm")
# }

meta.analyses <- cont.meta.list[[3]]
meta.analyses.limitmeta <- cont.meta.list[[4]]
meta.analyses.copas <- cont.meta.list[[5]]
meta.analyses.trimfill <- cont.meta.list[[6]]
meta.tests.linreg <- cont.meta.list[[7]]
meta.tests.begg <- cont.meta.list[[8]]
meta.tests.mm <- cont.meta.list[[9]]

cont.meta.list <- list(tmp, meta.id.vector, meta.analyses, 
											 meta.analyses.limitmeta,
											 meta.analyses.copas,
											 meta.analyses.trimfill,
											 meta.tests.linreg,
											 meta.tests.begg,
											 meta.tests.mm)

#Extract from lists:
cont.results <- data.frame(
	meta.id = meta.id.vector,
	meta.es.measure = rep("std. mean difference", times = length(meta.id.vector)),
	
	n.sig.single2 = unlist(lapply(meta.analyses, FUN = function(meta.analysis){length(which(meta.analysis$pval < 0.05))})),
	
	#Meta-Analysis
	est.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.fixed})),
	est.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.random})),
	
	zval.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$zval.fixed})),
	zval.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$zval.random})),
	
	pval.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$pval.fixed})),
	pval.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$pval.random})),
	
	se.est.fixef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.fixed})),
	se.est.ranef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.random})),
	
	Q = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$Q})),
	pval.Q = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$pval.Q})),
	#I2 = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$I2.w})),
	tau = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$tau})),
	#se.tau = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$se.tau})),
	#sparse = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$sparse})),
	k = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$k})),
	
	#Adjustment:
	method.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$method.adjust})),
	est.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
	se.est.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
	zval.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$zval.adjust})),
	pval.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
	alpha.r = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
	beta.r = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
	Q.small = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$Q.small})),
	Q.resid = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$Q.resid})),
	G.squared = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$G.squared})),
	
	est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
	se.est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})),
	missing.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[3]})),
	
	#Tests:
	pval.egger = unlist(lapply(meta.tests.linreg, FUN = function(meta.test){meta.test$p.value})),
	stat.egger = unlist(lapply(meta.tests.linreg, FUN = function(meta.test){meta.test$statistic})),
	pval.begg = unlist(lapply(meta.tests.begg, FUN = function(meta.test){meta.test$p.value})),
	stat.begg = unlist(lapply(meta.tests.begg, FUN = function(meta.test){meta.test$statistic})),
	pval.thompson = unlist(lapply(meta.tests.mm, FUN = function(meta.test){meta.test$p.value})),
	stat.thompson = unlist(lapply(meta.tests.mm, FUN = function(meta.test){meta.test$statistic}))
)


##########################################################################################
##########################################################################################
#Survival outcomes:
##########################################################################################
##########################################################################################

tmp <- data.ext2 %>% filter(outcome.type == "surv") %>% group_by(meta.id) %>% 
	mutate(n = n()) %>% filter(n >= 10) 

meta.id.vector <- unique(tmp$meta.id)
# meta.analyses <- list()
# meta.analyses.limitmeta <- list()
# meta.analyses.copas <- list()
# meta.analyses.trimfill <- list()
# meta.tests.linreg <- list()
# meta.tests.begg <- list()
# meta.tests.mm <- list()
# counter <- 0
# 
# 
# for(u in meta.id.vector){
# 	counter <- counter + 1
# 	print(c(u, counter))
# 	meta.analyses[[counter]] <- metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR", data = tmp[tmp$meta.id == u,])
# 
# 	meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
# 	meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
# 	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
# 	
# 	meta.tests.linreg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "linreg", k.min = 2)
# 	meta.tests.begg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "rank", k.min = 2)
# 	meta.tests.mm[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "mm", k.min = 2)
# }

meta.analyses <- surv.meta.list[[3]]
meta.analyses.limitmeta <- surv.meta.list[[4]]
meta.analyses.copas <- surv.meta.list[[5]]
meta.analyses.trimfill <- surv.meta.list[[6]]
meta.tests.linreg <- surv.meta.list[[7]]
meta.tests.begg <- surv.meta.list[[8]]
meta.tests.mm <- surv.meta.list[[9]]

surv.meta.list <- list(tmp, meta.id.vector, meta.analyses, 
											 meta.analyses.limitmeta,
											 meta.analyses.copas,
											 meta.analyses.trimfill,
											 meta.tests.linreg,
											 meta.tests.begg,
											 meta.tests.mm)

#Extract from lists:
surv.results <- data.frame(
	meta.id = meta.id.vector,
	meta.es.measure = rep(times = length(meta.id.vector), "log hazard ratio"),
	
	n.sig.single2 = unlist(lapply(meta.analyses, FUN = function(meta.analysis){length(which(meta.analysis$pval < 0.05))})),
	#Meta-Analysis
	est.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.fixed})),
	est.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.random})),
	
	zval.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$zval.fixed})),
	zval.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$zval.random})),
	
	pval.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$pval.fixed})),
	pval.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$pval.random})),
	
	se.est.fixef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.fixed})),
	se.est.ranef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.random})),
	
	Q = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$Q})),
	pval.Q = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$pval.Q})),
	# I2 = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$I2.w})),
	tau = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$tau})),
	# se.tau = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$se.tau})),
	# sparse = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$sparse})),
	k = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$k})),
	
	#Adjustment:
	method.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$method.adjust})),
	est.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
	se.est.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
	zval.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$zval.adjust})),
	pval.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
	alpha.r = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
	beta.r = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
	Q.small = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$Q.small})),
	Q.resid = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$Q.resid})),
	G.squared = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$G.squared})),
	
	est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
	se.est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})),
	missing.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[3]})),
	
	#Tests:
	pval.egger = unlist(lapply(meta.tests.linreg, FUN = function(meta.test){meta.test$p.value})),
	stat.egger = unlist(lapply(meta.tests.linreg, FUN = function(meta.test){meta.test$statistic})),
	pval.begg = unlist(lapply(meta.tests.begg, FUN = function(meta.test){meta.test$p.value})),
	stat.begg = unlist(lapply(meta.tests.begg, FUN = function(meta.test){meta.test$statistic})),
	pval.thompson = unlist(lapply(meta.tests.mm, FUN = function(meta.test){meta.test$p.value})),
	stat.thompson = unlist(lapply(meta.tests.mm, FUN = function(meta.test){meta.test$statistic}))
)

####################################################################################
####################################################################################
# Build complete dataset:
####################################################################################
####################################################################################

meta.results <- bind_rows(bin.results, cont.results, surv.results, .id = NULL)

meta.results <- meta.results %>% rowwise() %>% 
	mutate(I2  = max(0, (Q - k + 1)/Q)) #Proportion of variance of study estimates that is due to real variance)

meta.info <- data.ext2 %>%
	group_by(meta.id) %>% 
	mutate(n = n()) %>% filter(n >= 10) %>% 
	summarize(doi = unique(doi),
						outcome.type = unique(outcome.type),
						file.nr = unique(file.nr),
						comparison.nr = unique(comparison.nr),
						comparison.name = unique(comparison.name),
						outcome.name = unique(outcome.name),
						outcome.nr = unique(outcome.nr),
						subgroup.name = unique(subgroup.name),
						subgroup.nr = unique(subgroup.nr),
						outcome.measure.new = unique(outcome.measure.new),
						outcome.measure.new = unique(outcome.measure),
						n = n(), #Number of studies in meta-analysis
						mean.samplesize = mean(total1 + total2, na.rm = T),
						total.samplesize = sum(total1 + total2),
						var.samplesize = var(total1 + total2, na.rm = T),
						min.samplesize = min(min(total1, na.rm = T), min(total2, na.rm = T)),
						total.events = sum(events1 + events2),
						mean.events = mean(events1 + events2),
						mean.publication.year = mean(study.year, na.rm = TRUE),
						first.publication.year = min(study.year, na.rm = T),
						n.sig.single = sum(sig.single, na.rm = T),
						NA.sig.single = sum(is.na(sig.single)),
						se.min = min(se, na.rm = T),
						se.max = max(se, na.rm = T))

meta <- merge(meta.results, meta.info, by = c("meta.id"))
metac <- meta %>% filter(k > 9) %>% filter(n.sig.single > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)

# meta %>% filter(n.sig.single != n.sig.single2) %>% mutate(x = n.sig.single - n.sig.single2) %>% 
# 	select(meta.id, n.sig.single, n.sig.single2, outcome.type, outcome.measure.new, x)


##########################################################################################
##########################################################################################
#PEARSON CORRELATION META-ANALYSIS:
##########################################################################################
##########################################################################################



# tmp <- data.ext2 %>% group_by(meta.id) %>% mutate(n = n()) %>% filter(n > 9) %>% 
# 	filter(outcome.type != "rate" & !is.na(outcome.type)) # %>% filter(file.nr != 3014 & file.nr != 208)

# tmp <- tmp %>% filter(meta.id == 8361 | meta.id == 8363)
# tmp <- tmp %>% filter(meta.id == 378)

#Meta-analysis and ajustment part:
# meta.id.vector <- unique(tmp$meta.id)
meta.id.vector <- metac$meta.id[-which(metac$outcome.type == "surv")]
meta.id.vector <- meta.id.vector[-c(which(meta.id.vector == 157083),
																which(meta.id.vector == 159329),
																# which(meta.id.vector == 12459),
																# which(meta.id.vector == 12465),
																# which(meta.id.vector == 22087),
																which(meta.id.vector == 182298),
																which(meta.id.vector == 159329),
																which(meta.id.vector == 159329),
																which(meta.id.vector == 159329))] #157083 is one study and has total1 = 5 throughout.., 159329 gives issues with copas, 12459, 12465, 22087 (0 non NA cases)

meta.analyses <- list()
meta.analyses.limitmeta <- list()
meta.analyses.copas <- list()
counter <- 0

for(u in meta.id.vector){
	counter <- counter + 1
	print(c(u,counter))
	# meta.analyses[[counter]] <- metacor(cor = cor.pearson, n = total1 + total2, studlab = study.name, data = tmp[tmp$meta.id == u,])
	meta.analyses[[counter]] <- metagen(TE = cor.pearson, seTE = sqrt(var.cor.pearson), studlab = study.name, tmp[tmp$meta.id == u,])
	meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
}

# meta.analyses <- cor.meta.list[[3]]
# meta.analyses.limitmeta <- cor.meta.list[[4]]
# meta.analyses.copas <- cor.meta.list[[5]]

cor.meta.list <- list(tmp, meta.id.vector, meta.analyses, meta.analyses.limitmeta, meta.analyses.copas)

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

meta.f <- merge(metac, correlation.meta.analysis.estimates, by = c("meta.id"), all = T)

# correlation.meta.analysis.estimates <- data.frame(correlation.meta.analysis.estimates) %>% 
# 	mutate(est.z.fixef = 0.5 * log( (1 + est.cor.fixef)/(1 - est.cor.fixef) , base = exp(1)),
# 				 est.z.ranef = 0.5 * log( (1 + est.cor.ranef)/(1 - est.cor.ranef) , base = exp(1)),
# 				 est.z.reg = 0.5 * log( (1 + est.cor.reg)/(1 - est.cor.reg), base = exp(1)),
# 				 est.z.copas = 0.5 * log( (1 + est.cor.copas)/(1 - est.cor.copas) , base = exp(1)),
# 				 se.est.z.fixef = ifelse(total1 + total2 > 4))



##########################################################################################
##########################################################################################
#COHENS'd META-ANALYSIS:
##########################################################################################
##########################################################################################

metagen.bincont <- function(data){
  if (all(data$outcome.type == "bin")){
    meta <- metagen(TE = smd.ordl, seTE = sqrt(var.smd.ordl), studlab = study.name, tmp[tmp$meta.id == u,], sm = "SMD")
  } else{
    meta <- metagen(TE = smd, seTE = sqrt(var.smd), studlab = study.name, tmp[tmp$meta.id == u,], sm = "SMD")
  }
}

# meta.id.vector <- unique(tmp$meta.id)
meta.id.vector <- metac$meta.id[-which(metac$outcome.type == "surv")]
# meta.id.vector <- meta.id.vector[-c(which(meta.id.vector == 157083), which(meta.id.vector ==159329))]
meta.analyses <- list()
meta.analyses.limitmeta <- list()
meta.analyses.copas <- list()
counter <- 0

for(u in meta.id.vector){
  counter <- counter + 1
  print(c(u,counter))
  meta.analyses[[counter]] <- metagen.bincont(tmp[tmp$meta.id == u,])
  meta.analyses.limitmeta[[counter]] <- limitmeta(meta.analyses[[counter]])
  meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
}


# meta.analyses <- zscore.meta.list[[3]]
# meta.analyses.limitmeta <- zscore.meta.list[[4]]
# meta.analyses.copas <- zscore.meta.list[[5]]

cohensd.meta.list <- list(tmp, meta.id.vector, meta.analyses, meta.analyses.limitmeta, meta.analyses.copas)

#List extraction:
cohensd.meta.analysis.estimates <- cbind(
  meta.id = meta.id.vector,
  est.d.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.fixed})),
  est.d.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.random})),
  se.est.d.fixef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.fixed})),
  se.est.d.ranef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.random})),
  
  est.d.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
  se.est.d.reg = unlist(lapply(meta.analyses.limitmeta, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
  
  est.d.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
  se.est.d.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})))


meta.f <- merge(meta.f, cohensd.meta.analysis.estimates, by = c("meta.id"), all.x = T)


##########################################################################################
##########################################################################################
#Z-SCORE META-ANALYSIS:
##########################################################################################
##########################################################################################

# meta.id.vector <- unique(tmp$meta.id)
meta.id.vector <- metac$meta.id[-which(metac$outcome.type == "surv")]
meta.id.vector <- meta.id.vector[-c(which(meta.id.vector == 157083), which(meta.id.vector ==159329))]
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


# meta.analyses <- zscore.meta.list[[3]]
# meta.analyses.limitmeta <- zscore.meta.list[[4]]
# meta.analyses.copas <- zscore.meta.list[[5]]


zscore.meta.list <- list(tmp, meta.id.vector, meta.analyses, meta.analyses.limitmeta, meta.analyses.copas)

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


meta.f <- merge(meta.f, zscore.meta.analysis.estimates, by = c("meta.id"), all.x = T)
meta.bin <- meta.f %>% filter(outcome.type == "bin")
meta.cont <- meta.f %>% filter(outcome.type == "cont")
meta.surv <- meta.f %>% filter(outcome.type == "surv")





# #Cumulative meta-analysis
# meta.id.vector <- meta$meta.id
# meta.id.vector <- meta.id.vector[-c(which(meta.id.vector == 8395), which(meta.id.vector == 17017), 
#                                     which(meta.id.vector == 17019), which(meta.id.vector == 17019),
#                                     which(meta.id.vector == 26028),
#                                     which(meta.id.vector == 43408),
#                                     which(meta.id.vector == 82556 ),
#                                     which(meta.id.vector == 87619),
#                                     which(meta.id.vector == 87621),
#                                     which(meta.id.vector == 87636),
#                                     which(meta.id.vector == 183322),
#                                     which(meta.id.vector == 87621),
#                                     which(meta.id.vector == 87621),
#                                     which(meta.id.vector == 87621)
#                                     
#                                     )] #No z and var.z
# #[-c(which(meta.id.vector == 157083), which(meta.id.vector ==159329))]
# meta.analyses <- list()
# meta.analyses.cum <- list()
# counter <- 0
# 
# # tpp <- data.ext2 %>% filter(meta.id == meta$meta.id[778])
# # mt <- metagen(TE = z, seTE = sqrt(var.z), studlab = study.name, tpp)
# # xx <- metacum(mt, sortvar = study.year)
# 
# meta.analyses <- cum.meta.list[[3]]
# meta.analyses.cum <- cum.meta.list[[4]]
# 
# # for(u in meta.id.vector){
# # 	counter <- counter + 1
# # 	print(c(u,counter))
# # 	meta.analyses[[counter]] <- metagen(TE = z, seTE = sqrt(var.z), studlab = study.name, tmp[tmp$meta.id == u,])
# # 	meta.analyses.cum[[counter]] <- metacum(meta.analyses[[counter]], sortvar = study.year)
# # }
# # 
# # cum.meta.list <- list(meta.id.vector, tmp, meta.analyses, meta.analyses.cum)
# 
# z.cum <- c()
# se.z.cum <- c()
# study.years <- c()
# counter <- 0
# 
# for(u in meta.id.vector){
#   counter <- counter + 1
#   print(c(u,counter))
#   z.cum <- c(z.cum, meta.analyses.cum[[counter]]$TE)
#   se.z.cum <- c(se.z.cum, meta.analyses.cum[[counter]]$seTE)
#   if(length(meta.analyses.cum[[counter]]$TE) == length(sort(meta.analyses[[counter]]$data$study.year))){
#     study.years <- c(study.years, sort(meta.analyses[[counter]]$data$study.year))
#   } else{ #Complement with NA's
#     n.NAs <- length(meta.analyses.cum[[counter]]$TE) - length(sort(meta.analyses[[counter]]$data$study.year))
#     study.years <- c(study.years, c(sort(meta.analyses[[counter]]$data$study.year), rep(times = n.NAs, NA)))
#   }
#   
# }
# 
# 
# # study.years.centered <- scale(study.years, center = T, scale = F)
# # plot(y =abs(z.cum), study.years.centered, cex = .2)
# # # abline(loess(abs(z.cum)~study.years))
# # data.frame(study.years.centered, z.cum) %>% ggplot(aes(x = study.years.centered, y = abs(z.cum))) +
# #   geom_jitter(size = .5, alpha = .1) + geom_smooth() + geom_smooth(method = "lm", col = "red")
# # data.frame(study.years.centered, z.cum) %>% ggplot(aes(x = se.z.cum, y = abs(z.cum))) +
# #   geom_jitter(size = .5, alpha = .1) + geom_smooth() + geom_smooth(method = "lm", col = "red")
# # 
# # study.years.scaled <- scale(study.years, center = T, scale = T)
# # m1 <- lm(abs(z.cum)~study.years.scaled)
# # anova(m1)
# # m2 <- lm(abs(z.cum)~study.years.scaled*study.years)
# # anova(m1, m2)
# # m3 <- lm(abs(z.cum)~study.years.scaled*study.years + se.z.cum)
# # anova(m1, m2, m3)
# # plot(y =abs(z.cum), study.years, cex = .2)
# # abline(lm(abs(z.cum)~study.years))
# # 
# # study.years <- scale(study.years, center = T, scale = T)
# # plot(y =abs(z.cum), study.years, cex = .2)
# # abline(lm(abs(z.cum)~study.years))


save(meta, file =  file.path(PATH_RESULTS, "meta_analyses_summary_bin_cont_surv.RData"))
save(meta.f, file =  file.path(PATH_RESULTS, "meta_analyses_summary_complete.RData"))
save(meta.bin, file =  file.path(PATH_RESULTS, "meta.bin.RData"))
save(meta.cont, file =  file.path(PATH_RESULTS, "meta.cont.RData"))
save(meta.surv, file =  file.path(PATH_RESULTS, "meta.surv.RData"))
save(data, file =  file.path(PATH_RESULTS, "data_used_for_pb.process.RData"))
save(data.ext2, file =  file.path(PATH_RESULTS, "data_used_for_analysis.RData"))
save(bin.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_bin.RData"))
save(cont.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_cont.RData"))
save(surv.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_surv.RData"))
save(cor.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_cor.RData"))
save(zscore.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_zscore.RData"))
save(zscore.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_cohensd.RData"))
save(cum.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_cum.RData"))
