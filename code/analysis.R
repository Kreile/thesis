######################################################################
######################################################################
#META-ANALYSES AND PUBLICATION BIAS TESTS AND ADJUSTMENTS:############
######################################################################
######################################################################

#Load data:
rm(list = ls())
PATH_HOME = path.expand("~") # user home
PATH = file.path(PATH_HOME, 'Data/PubBias')
PATH2 = file.path(PATH_HOME, 'PubBias')
FILE = 'cochrane_2019-07-04.csv'
PATH_DATA = file.path(PATH, 'data')
PATH_CODE = file.path(PATH2, 'code')
PATH_RESULTS = file.path(PATH2, 'results_new')
# PATH_FIGURES = file.path(PATH_RESULTS, 'figures')

source(file.path(PATH_CODE, 'PubBias_functions.R'))

load(file.path(PATH_DATA, "PubBias_2019-07-19.RData"))
load(file.path(PATH_RESULTS, "data_used_for_analysis.RData"))
# data.ext2 <- pb.process3(data)


#--------------------------------------------------------------------------------------------------------------------#
# PRE-ANALYSIS PART - EXCLUDE UNSUITABLE META-ANALYSES:
#--------------------------------------------------------------------------------------------------------------------#

#To skip pre-analysis:
load(file.path(PATH_RESULTS, "meta_id_vector.RData"))
#--------------------------------------------------------------------------------------------------------------------#

# # ---------- Uncomment to run pre-analysis --------------
# # Get exclusion critetia (I-squared, n. sig. findings, s.e. ratio, dupl.index):
# meta.info.bin <- data.ext2 %>% filter(outcome.flag == "DICH") %>%
#   group_by(meta.id) %>%
#   mutate(n.pre = n()) %>% filter(n.pre >= 10) %>%
#   filter(!all(events1 == 0) & !all(events2 == 0)) %>% #No events
#   mutate(se.lrr = sqrt(var.lrr)) %>%
#   summarize(dupl.remove = unique(dupl.remove),
#             id = unique(id),
#             outcome.desired = unique(outcome.desired),
#             n.sig.single = sum(sig.single, na.rm = T),
#             se.min = min(se.lrr, na.rm = T),
#             se.max = max(se.lrr, na.rm = T))
# 
# meta.info.cont <- data.ext2 %>% filter(outcome.flag == "CONT") %>%
#   group_by(meta.id) %>%
#   mutate(n.pre = n()) %>% filter(n.pre >= 10) %>%
#   filter(!all(mean1 == 0) & !all(mean2 == 0)) %>% #No means
#   filter(!all(sd1 == 0) | !all(sd2 == 0)) %>% #No sd's
#   summarize(dupl.remove = unique(dupl.remove),
#             id = unique(id),
#             outcome.desired = unique(outcome.desired),
#             n.sig.single = sum(sig.single, na.rm = T),
#             se.min = min(se, na.rm = T),
#             se.max = max(se, na.rm = T))
# 
# meta.info.iv <- data.ext2 %>%
#   group_by(meta.id) %>%
#   mutate(n.pre = n()) %>% filter(n.pre >= 10) %>%
#   filter(outcome.flag == "IV") %>%
#   summarize(dupl.remove = unique(dupl.remove),
#             id = unique(id),
#             outcome.desired = unique(outcome.desired),
#             n.sig.single = sum(sig.single, na.rm = T),
#             se.min = min(se, na.rm = T),
#             se.max = max(se, na.rm = T))
# 
# meta.info.pre <- rbind(meta.info.bin, meta.info.cont, meta.info.iv) #To save to control how many are excluded
# 
# meta.info.pre <- meta.info.pre %>%
#   filter(id %in% data.review$id[which(data.review$withdrawn == FALSE)]) #No withdrawn meta-analyses
# 
# 
# meta.info <- meta.info.pre %>% filter(n.sig.single > 0) %>% #At least one significant result
#   filter((se.max^2)/(se.min^2) > 4) %>% #variance ratio > 4
#   filter(dupl.remove == 0) %>% #no duplicates
#   filter(outcome.desired == "efficacy") #Only efficacy outcomes
# #--------------------------------------------------------------------------------------------------------------------#
# 
# #Meta-analysis to get I2 (random effects meta-analysis):
# meta.id.vector.org.ranef <- meta.info$meta.id #Do meta-analysis for each of those
# 
# #Not for tobe.excluded meta-analyses: tau^2.max = 1000 (paulemandel method)
# meta.org.ranef.list <- list()
# counter <- 0
# for(u in meta.id.vector.org.ranef){
#   counter <- counter + 1
#   temp <- data.ext2 %>% filter(meta.id == u)
#   print(c(u, counter))
#   meta.org.ranef.list[[counter]] <- meta.fct.ranef.mom(temp)
# }
# 
# exclusion.estimates <- cbind(
#   meta.id = meta.id.vector.org.ranef,
#   n = unlist(lapply(meta.org.ranef.list, FUN = function(meta.analysis){meta.analysis$k})),
#   Q = unlist(lapply(meta.org.ranef.list, FUN = function(meta.analysis){meta.analysis$QE})))
# 
# 
# meta.info.pre2 <- merge(meta.info, exclusion.estimates, by = "meta.id")
# meta.info.I2 <- meta.info.pre2 %>% rowwise() %>% mutate(I2  = max(0, (Q - n + 1)/Q)) #calculate I2
# meta.info <- meta.info.I2 %>% filter(n > 9) %>% filter(I2 < 0.5)
# meta.id.vector <- meta.info$meta.id #Vector of meta-ids that match the criteria.
# 
# save(meta.id.vector, file =  file.path(PATH_RESULTS, "meta_id_vector.RData"))
# save(meta.info.I2, file =  file.path(PATH_RESULTS, "meta_id_I2.RData"))
#--------------------------------------------------------------------------------------------------------------------#

#Get more extensive information from the original dataset about the meta-analyses:
meta.info.extended <- data.ext2 %>% 
  filter(meta.id %in% meta.id.vector) %>% 
  group_by(meta.id) %>%
  summarize(id = unique(id),
            outcome.flag = unique(outcome.flag),
            comparison.nr = unique(comparison.nr),
            comparison.name = unique(comparison.name),
            outcome.name = unique(outcome.name),
            outcome.nr = unique(outcome.nr),
            subgroup.name = unique(subgroup.name),
            subgroup.nr = unique(subgroup.nr),
            outcome.measure.merged = unique(outcome.measure.merged),
            outcome.measure = unique(outcome.measure),
            outcome.desired = unique(outcome.desired),
            mean.samplesize = mean(total1 + total2, na.rm = T),
            total.samplesize = sum(total1 + total2),
            var.samplesize = var(total1 + total2, na.rm = T),
            min.samplesize = min(min(total1, na.rm = T), min(total2, na.rm = T)),
            total.events = sum(events1 + events2),
            mean.events = mean(events1 + events2),
            mean.publication.year = mean(study.year, na.rm = TRUE),
            first.publication.year = min(study.year, na.rm = T),
            side = bias.side.fct2(outcome = outcome.flag, outcome.measure.merged, lrr = lrr, var.lrr = var.lrr, 
                                  smd = cohensd, var.smd = var.cohensd, effect = effect, se = se),
            n.sig.single = sum(sig.single, na.rm = T),
            NA.sig.single = sum(is.na(sig.single)),
            se.min = min(se, na.rm = T),
            se.max = max(se, na.rm = T))

######################################################################################
######################################################################################
######################################################################################
# ANALYSIS PART: META-ANALYSES, TESTS AND ADJUSTMENTS: ###############################
######################################################################################
######################################################################################
######################################################################################

settings.meta(method.tau = "PM", method = "Inverse")

#Load previously constructed lists (to skip analysis):
load(file.path(PATH_RESULTS, "meta_complete_list_bin.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_cont.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_iv.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_zscore.RData"))
load(file.path(PATH_RESULTS, "meta_complete_list_cohensd.RData"))



#--------------------------------------------------------------------------------------------------------------------#
# BINARY DATA META-ANALYSES: ----------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#

# # -------- Uncomment to run analysis ----------
# meta.id.vector.bin <- meta.info.extended$meta.id[which(meta.info.extended$outcome.flag == "DICH")]
# meta.analyses <- list()
# meta.analyses.asd <- list()
# meta.analyses.reg <- list()
# meta.analyses.copas <- list()
# meta.tests.harbord <- list()
# meta.tests.peter <- list()
# meta.tests.rucker.mm <- list()
# meta.tests.rucker.linreg <- list()
# meta.tests.schwarzer <- list()
# meta.tests.excess <- list()
# counter <- 0
# for(u in meta.id.vector.bin){
# 	counter <- counter + 1
# 	print(c(u, counter))
# 	meta.analyses[[counter]] <- metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2,
# 																			studlab = study.name, sm = "RR",
# 																			method = "Inverse", data = data.ext2[data.ext2$meta.id == u,])
# 	meta.analyses.asd[[counter]] <- metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2,
# 																			studlab = study.name, sm = "ASD",
# 																			method = "Inverse", data = data.ext2[data.ext2$meta.id == u,])
# 
# 	meta.analyses.reg[[counter]] <- limitmeta(meta.analyses[[counter]])
# 	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
# 
# 	meta.tests.excess[[counter]] <- tes.fct2(data = data.ext2[data.ext2$meta.id == u,])
# 	meta.tests.harbord[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "score", k.min = 2)
# 	meta.tests.peter[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "peters", k.min = 2)
# 	meta.tests.schwarzer[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "count", k.min = 2)
# 	meta.tests.rucker.mm[[counter]] <- metabias(meta.analyses.asd[[counter]], method.bias = "mm", k.min = 2)
# 	meta.tests.rucker.linreg[[counter]] <- metabias(meta.analyses.asd[[counter]], method.bias = "linreg", k.min = 2)
# }
# 
# counter.unknown.bin.copas.error <- 0
# for(u in meta.id.vector.bin){
# 	if(length(meta.analyses.copas[[counter]]) < 3){
# 		meta.analyses.copas[[counter]] <- c(NA, NA, NA)
# 		counter.unknown.bin.copas.error <- counter.unknown.bin.copas.error + 1
# 	}
# }
# print(c("unknown.bin.copas.errors = ", counter.unknown.bin.copas.error))
# 
# bin.meta.list <- list(data.ext2, meta.id.vector.bin, meta.analyses,
# 											meta.analyses.asd,
# 											meta.analyses.reg,
# 											meta.analyses.copas,
# 											meta.tests.excess,
# 											meta.tests.harbord,
# 											meta.tests.peter,
# 											meta.tests.rucker.mm,
# 											meta.tests.schwarzer,
# 											meta.tests.excess,
# 											meta.tests.rucker.linreg)
# save(bin.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_bin.RData"))
#--------------------------------------------------------------------------------------------------------------------#

#Load analysis data:
meta.id.vector.bin <- bin.meta.list[[2]]
meta.analyses <- bin.meta.list[[3]]
meta.analyses.reg <- bin.meta.list[[5]]
meta.analyses.copas <- bin.meta.list[[6]]
meta.analyses.trimfill <- bin.meta.list[[7]]
meta.tests.harbord <- bin.meta.list[[8]]
meta.tests.peter <- bin.meta.list[[9]]
meta.tests.rucker.mm <- bin.meta.list[[10]]
meta.tests.schwarzer <- bin.meta.list[[11]]
meta.tests.excess <- bin.meta.list[[12]]
meta.tests.rucker.linreg <- bin.meta.list[[13]]


#Extract from lists:
bin.results <- data.frame(
  meta.id = meta.id.vector.bin,
  meta.es.measure = rep(times = length(meta.id.vector.bin), "log risk ratio"),
  
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
  method.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$method.adjust})),
  est.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
  se.est.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
  zval.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$zval.adjust})),
  pval.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
  alpha.r = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
  beta.r = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
  Q.small = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$Q.small})),
  Q.resid = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$Q.resid})),
  G.squared = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$G.squared})),
  
  est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
  se.est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})),
  missing.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[3]})),
  
  #Tests;
  pval.d.tes.org = unlist(lapply(meta.tests.excess, FUN = function(object){object[3]})),
  pval.d.tes = unlist(lapply(meta.tests.excess, FUN = function(object){object[4]})),
  stat.d.tes = unlist(lapply(meta.tests.excess, FUN = function(object){object[2]})),
  excess.d = unlist(lapply(meta.tests.excess, FUN = function(object){object[5]})) - 
    unlist(lapply(meta.tests.excess, FUN = function(object){object[6]})),
  pval.harbord = unlist(lapply(meta.tests.harbord, FUN = function(meta.test){meta.test$p.value})),
  stat.harbord = unlist(lapply(meta.tests.harbord, FUN = function(meta.test){meta.test$statistic})),
  pval.peter = unlist(lapply(meta.tests.peter, FUN = function(meta.test){meta.test$p.value})),
  stat.peter = unlist(lapply(meta.tests.peter, FUN = function(meta.test){meta.test$statistic})),
  pval.rucker = unlist(lapply(meta.tests.rucker.mm, FUN = function(meta.test){meta.test$p.value})),
  stat.rucker = unlist(lapply(meta.tests.rucker.mm, FUN = function(meta.test){meta.test$statistic})),
  pval.rucker.linreg = unlist(lapply(meta.tests.rucker.linreg, FUN = function(meta.test){meta.test$p.value})),
  stat.rucker.linreg = unlist(lapply(meta.tests.rucker.linreg, FUN = function(meta.test){meta.test$statistic})),
  pval.schwarzer = unlist(lapply(meta.tests.schwarzer, FUN = function(meta.test){meta.test$p.value})),
  stat.schwarzer = unlist(lapply(meta.tests.schwarzer, FUN = function(meta.test){meta.test$statistic}))
)



#--------------------------------------------------------------------------------------------------------------------#
# CONTINUOUS DATA META-ANALYSES: ------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#


# # -------- Uncomment to run analysis ----------
# meta.id.vector.cont <- meta.info.extended$meta.id[which(meta.info.extended$outcome.flag == "CONT")]
# meta.analyses <- list()
# meta.analyses.reg <- list()
# meta.analyses.copas <- list()
# meta.analyses.trimfill <- list()
# meta.tests.excess <- list()
# meta.tests.linreg <- list()
# meta.tests.begg <- list()
# meta.tests.mm <- list()
# counter <- 0
# for(u in meta.id.vector.cont){
# 	counter <- counter + 1
# 	print(c(u, counter))
# 	outcome.measure.merged <- as.character(unique(data.ext2[data.ext2$meta.id == u,"outcome.measure.merged"])$outcome.measure.merged)
# 	meta.analyses[[counter]] <- metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2,
# 																			 mean.c = mean2, sd.c = sd2, sm = outcome.measure.merged, studlab = study.name,
# 																			 data = data.ext2[data.ext2$meta.id == u,])
# 	meta.analyses.reg[[counter]] <- limitmeta(meta.analyses[[counter]])
# 	meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
# 	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
# 
# 	meta.tests.excess[[counter]] <- tes.fct2(data = data.ext2[data.ext2$meta.id == u,])
# 	meta.tests.linreg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "linreg", k.min = 2)
# 	meta.tests.begg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "rank", k.min = 2)
# 	meta.tests.mm[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "mm", k.min = 2)
# }
# cont.meta.list <- list(data.ext2, meta.id.vector.cont, meta.analyses,
# 											 meta.analyses.reg,
# 											 meta.analyses.copas,
# 											 meta.analyses.trimfill,
# 											 meta.tests.linreg,
# 											 meta.tests.begg,
# 											 meta.tests.mm,
# 											 meta.tests.excess)
# save(cont.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_cont.RData"))
#--------------------------------------------------------------------------------------------------------------------#

#Load analysis data:
meta.id.vector.cont <- cont.meta.list[[2]]
meta.analyses <- cont.meta.list[[3]]
meta.analyses.reg <- cont.meta.list[[4]]
meta.analyses.copas <- cont.meta.list[[5]]
meta.analyses.trimfill <- cont.meta.list[[6]]
meta.tests.linreg <- cont.meta.list[[7]]
meta.tests.begg <- cont.meta.list[[8]]
meta.tests.mm <- cont.meta.list[[9]]
meta.tests.excess <- cont.meta.list[[10]]

#Extract from lists:
cont.results <- data.frame(
  meta.id = meta.id.vector.cont,
  meta.es.measure = rep("std. mean difference", times = length(meta.id.vector.cont)),
  
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
  method.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$method.adjust})),
  est.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
  se.est.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
  zval.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$zval.adjust})),
  pval.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
  alpha.r = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
  beta.r = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
  Q.small = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$Q.small})),
  Q.resid = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$Q.resid})),
  G.squared = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$G.squared})),
  
  est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
  se.est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})),
  missing.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[3]})),
  
  #Tests:
  pval.d.tes.org = unlist(lapply(meta.tests.excess, FUN = function(object){object[3]})),
  pval.d.tes = unlist(lapply(meta.tests.excess, FUN = function(object){object[4]})),
  stat.d.tes = unlist(lapply(meta.tests.excess, FUN = function(object){object[2]})),
  excess.d = unlist(lapply(meta.tests.excess, FUN = function(object){object[5]})) - 
    unlist(lapply(meta.tests.excess, FUN = function(object){object[6]})),
  pval.egger = unlist(lapply(meta.tests.linreg, FUN = function(meta.test){meta.test$p.value})),
  stat.egger = unlist(lapply(meta.tests.linreg, FUN = function(meta.test){meta.test$statistic})),
  pval.begg = unlist(lapply(meta.tests.begg, FUN = function(meta.test){meta.test$p.value})),
  stat.begg = unlist(lapply(meta.tests.begg, FUN = function(meta.test){meta.test$statistic})),
  pval.thompson = unlist(lapply(meta.tests.mm, FUN = function(meta.test){meta.test$p.value})),
  stat.thompson = unlist(lapply(meta.tests.mm, FUN = function(meta.test){meta.test$statistic}))
)

#--------------------------------------------------------------------------------------------------------------------#
# IV OUTCOME META-ANALYSIS: -----------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#


# # -------- Uncomment to run analysis ----------
# meta.id.vector.iv <- meta.info.extended$meta.id[which(meta.info.extended$outcome.flag == "IV")]
# meta.analyses <- list()
# meta.analyses.reg <- list()
# meta.analyses.copas <- list()
# meta.analyses.trimfill <- list()
# meta.tests.linreg <- list()
# meta.tests.begg <- list()
# meta.tests.mm <- list()
# meta.tests.excess <- list()
# counter <- 0
# for(u in meta.id.vector.iv){
# 	counter <- counter + 1
# 
# 	meta.analyses[[counter]] <- metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR",
#                            data = data.ext2[data.ext2$meta.id == u,])
# 	print(c(u, counter))
# 	# print(c(u, counter, unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.fixed}))[counter]))
# 	meta.analyses.reg[[counter]] <- limitmeta(meta.analyses[[counter]])
# 	meta.analyses.trimfill[[counter]] <- trimfill(meta.analyses[[counter]])
# 	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
# 
# 	meta.tests.excess[[counter]] <- tes.fct2(data = data.ext2[data.ext2$meta.id == u,])
# 	meta.tests.linreg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "linreg", k.min = 2)
# 	meta.tests.begg[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "rank", k.min = 2)
# 	meta.tests.mm[[counter]] <- metabias(meta.analyses[[counter]], method.bias = "mm", k.min = 2)
# }
# iv.meta.list <- list(data.ext2, meta.id.vector.iv, meta.analyses,
# 											 meta.analyses.reg,
# 											 meta.analyses.copas,
# 											 meta.analyses.trimfill,
# 											 meta.tests.linreg,
# 											 meta.tests.begg,
# 											 meta.tests.mm,
# 											 meta.tests.excess)
# save(iv.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_iv.RData"))
#--------------------------------------------------------------------------------------------------------------------#

#Load analysis data:
meta.id.vector.iv <- iv.meta.list[[2]]
meta.analyses <- iv.meta.list[[3]]
meta.analyses.reg <- iv.meta.list[[4]]
meta.analyses.copas <- iv.meta.list[[5]]
meta.analyses.trimfill <- iv.meta.list[[6]]
meta.tests.linreg <- iv.meta.list[[7]]
meta.tests.begg <- iv.meta.list[[8]]
meta.tests.mm <- iv.meta.list[[9]]
meta.tests.excess <- iv.meta.list[[10]]

#Extract from lists:
iv.results <- data.frame(
  meta.id = meta.id.vector.iv,
  meta.es.measure = rep(times = length(meta.id.vector.iv), "unknown.effect(IV)"),
  
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
  method.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$method.adjust})),
  est.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
  se.est.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
  zval.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$zval.adjust})),
  pval.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
  alpha.r = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
  beta.r = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$pval.adjust})),
  Q.small = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$Q.small})),
  Q.resid = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$Q.resid})),
  G.squared = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$G.squared})),
  
  est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
  se.est.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})),
  missing.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[3]})),
  
  #Tests:
  pval.d.tes.org = unlist(lapply(meta.tests.excess, FUN = function(object){object[3]})),
  pval.d.tes = unlist(lapply(meta.tests.excess, FUN = function(object){object[4]})),
  stat.d.tes = unlist(lapply(meta.tests.excess, FUN = function(object){object[2]})),
  excess.d = unlist(lapply(meta.tests.excess, FUN = function(object){object[5]})) - 
    unlist(lapply(meta.tests.excess, FUN = function(object){object[6]})),
  pval.egger = unlist(lapply(meta.tests.linreg, FUN = function(meta.test){meta.test$p.value})),
  stat.egger = unlist(lapply(meta.tests.linreg, FUN = function(meta.test){meta.test$statistic})),
  pval.begg = unlist(lapply(meta.tests.begg, FUN = function(meta.test){meta.test$p.value})),
  stat.begg = unlist(lapply(meta.tests.begg, FUN = function(meta.test){meta.test$statistic})),
  pval.thompson = unlist(lapply(meta.tests.mm, FUN = function(meta.test){meta.test$p.value})),
  stat.thompson = unlist(lapply(meta.tests.mm, FUN = function(meta.test){meta.test$statistic}))
)


#--------------------------------------------------------------------------------------------------------------------#
# Z-SCORE BASED META-ANALYSES: --------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#

# # -------- Uncomment to run analysis ----------
# meta.id.vector.noiv <- meta.info.extended$meta.id[-which(meta.info.extended$outcome.flag == "IV")]
# meta.id.vector.smd <- meta.info.extended$meta.id[which(meta.info.extended$outcome.measure.merged == "SMD")]
# meta.id.vector.z <- union(meta.id.vector.smd, meta.id.vector.noiv)
# meta.id.vector.z.c <- meta.id.vector.z[-c(which(meta.id.vector.z == 104613),
#                                           which(meta.id.vector.z == 104616),
#                                           which(meta.id.vector.z == 153452),
#                                           which(meta.id.vector.z == 163246),
#                                           which(meta.id.vector.z == 182958))] #Have zero total1..
# tmp.z <- data.ext2 %>% group_by(meta.id) %>% mutate(n = n()) %>% filter(n > 9) %>%
#   filter(total1 + total2 > 3) %>% filter(total1 != 0 & total2 != 0)
# meta.analyses <- list()
# meta.analyses.reg <- list()
# meta.analyses.copas <- list()
# counter <- 0
# for(u in meta.id.vector.z.c){
# 	counter <- counter + 1
# 	print(c(u,counter))
# 	meta.analyses[[counter]] <- metacor(cor = z, n = total1 + total2, studlab = study.name, tmp.z[tmp.z$meta.id == u,])
# 	meta.analyses.reg[[counter]] <- limitmeta(meta.analyses[[counter]])
# 	meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
# }
# zscore.meta.list <- list(tmp.z, meta.id.vector.z.c, meta.analyses, meta.analyses.reg, meta.analyses.copas)
# save(zscore.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_zscore.RData"))
#--------------------------------------------------------------------------------------------------------------------#

#Load analysis data:
meta.id.vector.z.c <- zscore.meta.list[[2]]
meta.analyses <- zscore.meta.list[[3]]
meta.analyses.reg <- zscore.meta.list[[4]]
meta.analyses.copas <- zscore.meta.list[[5]]

#List extraction:
zscore.meta.analysis.estimates <- cbind(
  meta.id = meta.id.vector.z.c,
  est.z.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.fixed})),
  est.z.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.random})),
  se.est.z.fixef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.fixed})),
  se.est.z.ranef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.random})),
  
  est.z.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
  se.est.z.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
  
  est.z.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
  se.est.z.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})))


#--------------------------------------------------------------------------------------------------------------------#
# HEDGES G BASED META-ANALYSES: -------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#

# # -------- Uncomment to run analysis ----------
# meta.id.vector.d <- meta.id.vector.z[-c(which(meta.id.vector.z == 25407))]
# meta.analyses <- list()
# meta.analyses.reg <- list()
# meta.analyses.copas <- list()
# counter <- 0
# for(u in meta.id.vector.d){
#   counter <- counter + 1
#   print(c(u,counter))
#   meta.analyses[[counter]] <- metagen.bincont2(data.ext2[data.ext2$meta.id == u,])
#   meta.analyses.reg[[counter]] <- limitmeta(meta.analyses[[counter]])
#   meta.analyses.copas[[counter]] <- auto.copas(meta.analyses[[counter]], sig.level = 0.1)
# }
# cohensd.meta.list <- list(data.ext2, meta.id.vector.d, meta.analyses, meta.analyses.reg, meta.analyses.copas)
# save(cohensd.meta.list, file =  file.path(PATH_RESULTS, "meta_complete_list_cohensd.RData"))
#--------------------------------------------------------------------------------------------------------------------#

#Load analysis data:
meta.id.vector.d <- cohensd.meta.list[[2]]
meta.analyses <- cohensd.meta.list[[3]]
meta.analyses.reg <- cohensd.meta.list[[4]]
meta.analyses.copas <- cohensd.meta.list[[5]]

#List extraction:
cohensd.meta.analysis.estimates <- cbind(
  meta.id = meta.id.vector.d,
  est.d.fixef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.fixed})),
  est.d.ranef = unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$TE.random})),
  se.est.d.fixef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.fixed})),
  se.est.d.ranef =  unlist(lapply(meta.analyses, FUN = function(meta.analysis){meta.analysis$seTE.random})),
  
  est.d.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$TE.adjust})),
  se.est.d.reg = unlist(lapply(meta.analyses.reg, FUN = function(meta.adjust){meta.adjust$seTE.adjust})),
  
  est.d.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[1]})),
  se.est.d.copas = unlist(lapply(meta.analyses.copas, FUN = function(meta.adjust){meta.adjust[2]})))


#--------------------------------------------------------------------------------------------------------------------#
# BUILD FINAL DATASET OF RESULTS
#--------------------------------------------------------------------------------------------------------------------#

meta.tp5 <- bind_rows(bin.results, cont.results, iv.results)
meta.tp1 <- merge(meta.info.extended, meta.tp5, by = c("meta.id")) #Temporary versions
meta.tp4 <- merge(by = "meta.id", x = meta.tp1, y = zscore.meta.analysis.estimates, all.x = T)
meta.f <- merge(by = "meta.id", x = meta.tp4, y = cohensd.meta.analysis.estimates, all.x = T)

#--------------------------------------------------------------------------------------------------------------------#
# ROUND-UP: DETECT AND REPLACE MISSING VALUES
#--------------------------------------------------------------------------------------------------------------------#

meta.f$se.est.copas.na <- ifelse(is.na(meta.f$se.est.copas), 1, 0) #se.est.copas.na = 1 -> is missing
meta.f$se.est.z.copas.na <- ifelse(is.na(meta.f$se.est.z.copas), 1, 0)
meta.f$se.est.d.copas.na <- ifelse(is.na(meta.f$se.est.d.copas), 1, 0)

meta.f$se.est.reg.na <- ifelse(is.na(meta.f$se.est.reg), 1, 0)
meta.f$se.est.z.reg.na <- ifelse(is.na(meta.f$se.est.z.reg), 1, 0)
meta.f$se.est.d.reg.na <- ifelse(is.na(meta.f$se.est.d.reg), 1, 0)

meta.f$est.copas.na <- ifelse(is.na(meta.f$est.copas), 1, 0)
meta.f$est.z.copas.na <- ifelse(is.na(meta.f$est.z.copas), 1, 0)
meta.f$est.d.copas.na <- ifelse(is.na(meta.f$est.d.copas), 1, 0)

meta.f$est.reg.na <- ifelse(is.na(meta.f$est.reg), 1, 0)
meta.f$est.z.reg.na <- ifelse(is.na(meta.f$est.z.reg), 1, 0)
meta.f$est.d.reg.na <- ifelse(is.na(meta.f$est.d.reg), 1, 0)
#--------------------------------------------------------------------------------------------------------------------#

#Inpute random effect estimate for missing copas:
copas.names <- c("est.copas", "se.est.copas", "est.z.copas", "se.est.z.copas", "est.d.copas", "se.est.d.copas")
ranef.names <- c("est.ranef", "se.est.ranef", "est.z.ranef", "se.est.z.ranef", "est.d.ranef", "se.est.d.ranef")
missing.names <- paste(copas.names, ".missing", sep = "")
missing <- c()

for(u in 1:length(copas.names)){
  missing.count <- 0
  for(k in 1:dim(meta.f)[1]){
    if(is.na(meta.f[k, copas.names[u]])){
      missing.count <- missing.count + 1
      meta.f[k, copas.names[u]] <- meta.f[k, ranef.names[u]]
    }
  }
  missing[u] <- missing.count
  meta.f[, missing.names[u]] <- missing.count
}
#--------------------------------------------------------------------------------------------------------------------#

#Some useful variables:
meta.f <- meta.f %>% rowwise() %>% mutate(pval1.egger = onesided.p(stat = stat.egger, side = side, n = k, test.type = "reg"),
                                          pval1.thompson = onesided.p(stat = stat.thompson, side = side, n = k, test.type = "reg"),
                                          pval1.begg = onesided.p(stat = stat.begg, side = side, n = k, test.type = "rank"),
                                          
                                          pval1.harbord = onesided.p(stat = stat.harbord, side = side, n = k, test.type = "reg"),
                                          pval1.rucker = onesided.p(stat = stat.rucker, side = side, n = k, test.type = "reg"),
                                          pval1.rucker.linreg = onesided.p(stat = stat.rucker.linreg, side = side, n = k, test.type = "reg"),
                                          pval1.peter = onesided.p(stat = stat.peter, side = side, n = k, test.type = "reg"),
                                          pval1.schwarzer = onesided.p(stat = stat.schwarzer, side = side, n = k, test.type = "rank"),
                                          pval1.d.tes = pval.d.tes)

meta.f <- meta.f %>% rowwise() %>% mutate(I2  = max(0, (Q - k + 1)/Q),
                            var.ratio = (se.max^2)/(se.min^2),
                            pval.se = case_when(outcome.flag == "DICH" ~ pval1.rucker.linreg, 
                                                TRUE ~ pval1.egger),
                            pval.se.het = case_when(outcome.flag == "DICH" ~ pval1.rucker, 
                                                    TRUE ~ pval1.thompson))

sig.level <- 0.1
meta.f <- meta.f %>% rowwise() %>% 
  mutate(egger.test = ifelse(pval1.egger < sig.level, 1, 0),
         thompson.test = ifelse(pval1.thompson < sig.level, 1, 0),
         begg.test = ifelse(pval1.begg < sig.level, 1, 0),
         
         tes.d.test = ifelse(pval1.d.tes < sig.level, 1, 0),
         
         schwarzer.test = ifelse(pval1.schwarzer < sig.level, 1, 0),
         rucker.test = ifelse(pval1.rucker < sig.level, 1, 0),
         rucker.test.linreg = ifelse(pval1.rucker.linreg < sig.level, 1, 0),
         harbord.test = ifelse(pval1.harbord < sig.level, 1, 0),
         peter.test = ifelse(pval1.peter < sig.level, 1, 0))
#--------------------------------------------------------------------------------------------------------------------#

#Separate outcome.types:
meta.bin <- meta.f %>% filter(outcome.flag == "DICH")
meta.cont <- meta.f %>% filter(outcome.flag == "CONT")
meta.iv <- meta.f %>% filter(outcome.flag == "IV")
#--------------------------------------------------------------------------------------------------------------------#

#Save the result datasets:
save(meta.f, file =  file.path(PATH_RESULTS, "meta_analyses_summary_complete.RData"))
save(meta.bin, file =  file.path(PATH_RESULTS, "meta.bin.RData"))
save(meta.cont, file =  file.path(PATH_RESULTS, "meta.cont.RData"))
save(meta.iv, file =  file.path(PATH_RESULTS, "meta.iv.RData"))
save(meta.id.vector, file =  file.path(PATH_RESULTS, "meta_id_vector.RData"))
save(data.ext2, file =  file.path(PATH_RESULTS, "data_used_for_analysis.RData")) #Save data used for analysis




# x <- foreach(i = 1:length(meta.id.vector.bin)) %dopar% 
#   # u <- meta.id.vector.bin[i]
#   # counter <- counter + 1
#   metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2,
#           studlab = study.name, sm = "RR",
#           method = "Inverse", data = data.ext2[data.ext2$meta.id == meta.id.vector.bin[i],])
# 
# cl<-makeCluster(7)
# registerDoParallel(cl)
# 
# foreach(i = 1:length(meta.id.vector.bin)) %dopar%
#   sqrt(i)
