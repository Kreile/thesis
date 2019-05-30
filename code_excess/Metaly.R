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

load(file.path(PATH_RESULTS, "meta_analyses_summary_complete.RData"))
load(file.path(PATH_RESULTS, "meta.bin.RData"))
load(file.path(PATH_RESULTS, "meta.cont.RData"))
load(file.path(PATH_RESULTS, "meta.surv.RData"))
load(file.path(PATH_RESULTS, "data_used_for_analysis.RData"))



#Explanation of coding: 0 = (unsig, unsig), 1 = (unsig, sig), 2 = (sig, sig), 3 = (sig, unsig)
sig.level <- 0.05
meta.f <- meta.f %>% 
	mutate(
		pval.copas = 2*(1-pnorm(abs(est.copas/se.est.copas))),
		
		sig.fixef = ifelse(pval.fixef > sig.level, 0, 1),
		sig.ranef = ifelse(pval.ranef > sig.level, 0, 1),
		sig.reg.ranef = ifelse(pval.reg > sig.level, 0, 1),
		sig.copas = ifelse(pval.copas > sig.level, 0, 1),
		
		sig.change.ranef.reg = sig.change(sig.before = sig.ranef, sig.after =  sig.reg), #Function to register if significance of estimate has changed after correction
		sig.change.ranef.copas = sig.change(sig.before = sig.ranef, sig.after =  sig.copas),
		sig.change.fixef.reg = sig.change(sig.before = sig.fixef, sig.after =  sig.reg),
		sig.change.fixef.copas = sig.change(sig.before = sig.fixef, sig.after =  sig.copas))


sig.level <- 0.1
meta.f <- meta.f %>% 
	mutate(egger.test = ifelse(pval.peter < sig.level, 1, 0),
				 thompson.test = ifelse(pval.thompson < sig.level, 1, 0),
				 begg.test = ifelse(pval.begg < sig.level, 1, 0),
				 
				 schwarzer.test = ifelse(pval.schwarzer < sig.level, 1, 0),
				 rucker.test = ifelse(pval.rucker < sig.level, 1, 0),
				 harbord.test = ifelse(pval.harbord < sig.level, 1, 0),
				 peter.test = ifelse(pval.peter < sig.level, 1, 0))

meta.bin <- meta.f %>% filter(outcome.type == "bin")
meta.cont <- meta.f %>% filter(outcome.type == "cont")
meta.surv <- meta.f %>% filter(outcome.type == "surv")

#Comparison of "original effect sizes":
effect.diff <- meta.f %>% mutate(
	est.fixef = ifelse(outcome.type == "bin", exp(est.fixef), est.fixef),
	est.ranef = ifelse(outcome.type == "bin", exp(est.ranef), est.ranef),
	est.copas = ifelse(outcome.type == "bin", exp(est.copas), est.copas),
	est.reg = ifelse(outcome.type == "bin", exp(est.reg.ranef), est.reg.ranef)
)

effect.diff <- effect.diff %>% mutate(
	est.fixef = ifelse(outcome.type == "bin", est.fixef - 1, est.fixef),
	est.ranef = ifelse(outcome.type == "bin", est.ranef - 1, est.ranef),
	est.copas = ifelse(outcome.type == "bin", est.copas - 1, est.copas),
	est.reg = ifelse(outcome.type == "bin", est.reg.ranef - 1, est.reg.ranef),
	fixef.ranef = ifelse(abs(est.ranef) > abs(est.fixef), "Reduction", "Amplification"),
	fixef.copas = ifelse(abs(est.fixef) > abs(est.copas), "Reduction", "Amplification"),
	fixef.reg = ifelse(abs(est.fixef) > abs(est.reg), "Reduction", "Amplification"),
	fixef.copas.s = ifelse(sign(est.fixef) == sign(est.copas), "Unchanged", "Reversed"),
	fixef.reg.s = ifelse(sign(est.fixef) == sign(est.reg), "Unchanged", "Reversed")) 






mly.bin = function(data, sig.level, min.study.number) {
  meta.bin <- data %>% filter(outcome.measure.new == "Risk Ratio" | outcome.measure.new == "Odds Ratio" | outcome.measure.new == "Risk Difference"
                              | outcome.measure.new == "Peto Odds Ratio") %>%  
    filter(file.nr != 3014) %>% #file.nr 3014 gibt eine seltsame Fehlermeldung, muss ausgeschlossen werden
    filter(events1 > 0 | events2 > 0) %>% #comparisons mit events  = 0 ausschliessen
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>% #comparisons mit events = total ausschliessen - vertraegt der alg. irgendwie nicht.
    filter(total1 > 12 & total2 > 12) %>% 
    group_by(meta.id) %>% #Die nachfolgenden Rechnungen werden jeweils fuer diese Gruppen gemacht.
    mutate(n = n()) %>% filter(n >= min.study.number) %>% #Nur meta-analysen mit mehr als 9 comparisons behalten
    mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
    mutate(study.year = ifelse(study.year > 1920, study.year, NA)) %>% 
    summarize(doi = unique(doi),
              n = n(),
              #outcome.type = unique(outcome.type),
              comparison.name = unique(comparison.name),
              outcome.name = unique(outcome.name),
              sungroup.name = unique(subgroup.name),
              outcome.measure.new = unique(outcome.measure.new),
              mean.samplesize = mean(total1 + total2, na.rm = T),
              total.samplesize = sum(total1 + total2),
              var.samplesize = var(total1 + total2, na.rm = T),
              total.events = sum(events1 + events2),
              mean.events = mean(events1 + events2),
              mean.publication.year = mean(study.year, na.rm = TRUE),
              first.publication.year = min(study.year, na.rm = T),
              n.sig.type.bin = sum(sig.type, na.rm = T),
              n.sig.fishersz.bin = sum(sig.fishersz, na.rm = T),
              mean.sig.type = mean(sig.type, na.rm = T),
              NA.sig.type = sum(is.na(sig.type)),
              mean.sig.fishersz = mean(sig.fishersz, na.rm = T),
              NA.sig.fishersz = sum(is.na(sig.fishersz)),
              lor.fixef.bin = #Pooled odds ratio estimate of fixed effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.fixed),
              se.lor.fixef.bin = #Pooled odds ratio estimate of fixed effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$seTE.fixed),
              pval.fixef.bin = #P-value
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.fixed,
              sig.fixef.bin = ifelse(pval.fixef.bin > sig.level, 0, 1), #If significant
              lor.ranef.bin = #Pooled odds ratio estimate of random effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.random),
              se.lor.ranef.bin = #Pooled odds ratio estimate of random effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.random),
              pval.ranef.bin = 
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.random,
              sig.ranef.bin =  ifelse(pval.ranef.bin > sig.level, 0, 1),
              
              #Hartung Knapp Modification
              lor.ranef.hkn.bin = #Pooled odds ratio estimate of random effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", hakn = T)$TE.random),
              se.lor.ranef.hkn.bin = #Pooled odds ratio estimate of random effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", hakn = T)$TE.random),
              pval.ranef.hkn.bin = 
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", hakn = T)$pval.random,
              sig.ranef.hkn.bin = ifelse(pval.ranef.hkn.bin > sig.level, 0, 1),
              
              Q.bin = #Q between studies heterogeneity statistic
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$Q,
              pval.Q.bin = #P-value for between study heterogeneity statistic = 0
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$pval.Q,
              sig.Q.bin = ifelse(pval.Q.bin < sig.level, 1, 0),
              I2  = max(0, (Q.bin - n + 1)/Q.bin), #Proportion of variance of study estimates that is due to real variance
              sig.change.ranef.reg = sig.change(sig.fixef.bin, sig.after =  sig.ranef.bin)) #Function to register if significance of estimate has changed after correction
  
  return(meta.bin)
}


file.bin <- "mly.bin.RData"
if (file.exists(file.path(PATH_RESULTS, file.bin))) {
  load(file.path(PATH_RESULTS, file.bin))
} else { 
  data.bin <- mly.bin(data.ext, 0.05, min.study.number = 2)
  save(data.bin, file =  file.path(PATH_RESULTS, file.bin))
}


mly.cont = function(data, sig.level, min.study.number) {
  meta.cont <- data %>% filter(outcome.measure.new == "Mean Difference" | outcome.measure.new == "Std. Mean Difference" | 
                                 outcome.measure.new == "Hedges' G") %>%
    filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
    filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
    filter(total1 > 12 & total2 > 12) %>% 
    group_by(meta.id) %>% #Die nachfolgenden Rechnungen werden jeweils fuer diese Gruppen gemacht.
    mutate(n = n()) %>% filter(n >= min.study.number) %>% #Nur meta-analysen mit mehr als 9 comparisons behalten
    mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
    mutate(study.year = ifelse(study.year > 1920, study.year, NA)) %>% 
    summarize(doi = unique(doi),
              #outcome.type = unique(outcome.type),
              comparison.name = unique(comparison.name),
              outcome.name = unique(outcome.name),
              sungroup.name = unique(subgroup.name),
              outcome.measure.new = unique(outcome.measure.new),
              n = n(), #Number of studies in meta-analysis
              mean.samplesize = mean(total1 + total2, na.rm = T),
              total.samplesize = sum(total1 + total2),
              var.samplesize = var(total1 + total2, na.rm = T),
              total.events = sum(events1 + events2),
              mean.events = mean(events1 + events2),
              mean.publication.year = mean(study.year, na.rm = TRUE),
              mean.sig.type = mean(sig.type, na.rm = T),
              mean.sig.fishersz = mean(sig.fishersz, na.rm = T),
              n.sig.type.cont = sum(sig.type, na.rm = T),
              NA.sig.type = sum(is.na(sig.type)),
              n.sig.fishersz.cont = sum(sig.fishersz, na.rm = T),
              NA.sig.fishersz = sum(is.na(sig.fishersz)),
              first.publication.year = min(study.year, na.rm = T),
              
              est.fixef.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$TE.fixed,
              se.fixef.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$seTE.fixed,
              pval.fixef.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$pval.fixed,
              sig.fixef.cont = ifelse(pval.fixef.cont > sig.level, 0, 1),
              est.ranef.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$TE.random,
              se.ranef.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$seTE.random,
              pval.ranef.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$pval.random,
              sig.ranef.cont = ifelse(pval.ranef.cont > sig.level, 0, 1),
              
              #Hartung Knapp Method
              est.ranef.hkn.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, hakn = T)$TE.random,
              se.ranef.hkn.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, hakn = T)$seTE.random,
              pval.ranef.hkn.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, hakn = T)$pval.random,
              sig.ranef.hkn.cont = ifelse(pval.ranef.hkn.cont > sig.level, 0, 1),
              
              Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$Q,
              pval.Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$pval.Q,
              sig.Q.cont = ifelse(pval.Q.cont < sig.level, 1, 0),
              I2  = max(0, (Q.cont - n + 1)/Q.cont))
  return(meta.cont)
}

file.cont <- "mly.cont.RData"
if (file.exists(file.path(PATH_RESULTS, file.cont))) {
  load(file.path(PATH_RESULTS, file.cont))
} else { 
  data.cont <- mly.cont(data.ext, 0.05, min.study.number = 2)
  save(data.cont, file =  file.path(PATH_RESULTS, file.cont))
}


mly.merge <- function(mly.result.bin, mly.result.cont){
  mly.result <- bind_rows(mly.result.bin, mly.result.cont)
  mly.result <- mly.result %>% mutate(est.fixef = ifelse(!is.na(lor.fixef.bin), lor.fixef.bin, est.fixef.cont),
                                      est.ranef = ifelse(!is.na(lor.ranef.bin), lor.ranef.bin, est.ranef.cont),
                                      pval.fixef = ifelse(!is.na(pval.fixef.bin), pval.fixef.bin, pval.fixef.cont),
                                      pval.ranef = ifelse(!is.na(pval.ranef.bin), pval.ranef.bin, pval.ranef.cont),
                                      sig.fixef = ifelse(!is.na(sig.fixef.bin), sig.fixef.bin, sig.fixef.cont),
                                      sig.ranef = ifelse(!is.na(sig.ranef.bin), sig.ranef.bin, sig.fixef.cont),
                                      est.ranef.hkn = ifelse(!is.na(lor.ranef.hkn.bin), lor.ranef.hkn.bin, est.ranef.hkn.cont),
                                      pval.ranef.hkn = ifelse(!is.na(pval.ranef.hkn.bin), pval.ranef.hkn.bin, pval.ranef.hkn.cont),
                                      sig.ranef.hkn = ifelse(!is.na(sig.ranef.hkn.bin), sig.ranef.hkn.bin, sig.ranef.hkn.cont),
                                      Q = ifelse(!is.na(Q.bin), Q.bin, Q.cont),
                                      pval.Q = ifelse(!is.na(pval.Q.bin), pval.Q.bin, pval.Q.cont),
                                      sig.Q = ifelse(!is.na(sig.Q.bin), sig.Q.bin, sig.Q.cont),
                                      n.sig.type = ifelse(!is.na(n.sig.type.bin), n.sig.type.bin, n.sig.type.cont),
                                      n.sig.fishersz = ifelse(!is.na(n.sig.fishersz.bin), n.sig.fishersz.bin, n.sig.fishersz.cont)
                                      #= ifelse(!is.na(.bin), .bin, .cont),
                                      )
  return(mly.result)
}

mly <- mly.merge(data.bin, data.cont)
save(mly, file =  file.path(PATH_RESULTS, "mly.RData"))






print(data.ext %>% group_by(file.nr, comparison.nr, outcome.name, subgroup.nr) %>% distinct(outcome.type), n = 400)
print(data.ext %>% filter(file.nr == 2 & outcome.name == "Resuturing of wound - up to 3 months" & subgroup.nr == 1) %>% select(outcome.measure, events1, events2, total1, total2))








#Tests
tmp <- data.ext %>% filter(outcome.measure.new == "Odds Ratio") %>% select(study.name, events1, total1, events2, total2, log.OR, se.log.OR, pval.type, effect, se)
dd <- tmp[40:45,]

mt <- metabin(event.e = events1, event.c = events2, n.e = total1, n.c = total2, studlab = study.name, data = dd, sm = "OR")

print(mt)
data.frame(mt$TE, mt$seTE)
dd %>% select(study.name, std.m.d, se.std.m.d, effect, se)

log(dd$effect)



tmp <- data.ext %>% filter(outcome.measure.new == "Mean Difference") %>% select(study.name, mean1, sd1, total1, mean2, sd2, total2, std.m.d, se.std.m.d, pval.type, effect, se)
dd <- tmp[40:45,]

mt <- metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, data = dd, sm = "SMD", method.smd = "Cohen")

print(mt)
data.frame(mt$TE, mt$seTE, mt$pval)
dd %>% select(study.name, std.m.d, se.std.m.d, effect, se, pval.type)
?metacont























data.ext %>% ggplot(aes(x = study.year)) + geom_histogram(binwidth = 1) + scale_x_continuous(limits = c(1980, 2020))

data.ext %>% filter(!is.na(sig.type)) %>%  ggplot(aes(x = study.year, fill = factor(sig.type), stat(count))) + geom_density(na.rm = T, position = "fill") + 
  labs(fill = "Significance of treatment effect") + scale_fill_discrete(labels= c("No", "Yes"))+ scale_x_continuous(limits = c(1970, 2017))

###############################################################################################################################################
###############################################################################################################################################
#Merge comparisons with meta-analysis result:
mly.tomerge <- mly %>% select(meta.id, sig.Q, n.sig.type, n.sig.fishersz, sig.fixef, sig.ranef, sig.ranef.hkn)

data.mly <- merge(x = data.ext, y = mly.tomerge, by = c("meta.id")) #ideally, 278626 comparisons ( = sum(mly$n)), 
# but 300831 (18261 have smaller sample size than mly, but same meta.id)

data.mly <- data.mly %>% filter(total1 > 11 & total2 > 11) #Now 282550



sig.meta.c <- data.mly %>% select(sig.fixef, sig.ranef, sig.ranef.hkn, sig.type) %>% 
  gather(key = "meta.analysis.method", value = "sig.meta.analysis") %>% mutate(sig.meta.analysis = factor(sig.meta.analysis))

sig.meta.count.c <- sig.meta.c %>% 
  group_by(sig.meta.analysis, meta.analysis.method) %>% count()

sig.meta.c %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("meta.nonsig", "meta.sig")) +
  annotate("text", x = sig.meta.count.c[sig.meta.count.c$sig.meta.analysis == 0,]$meta.analysis.method, y = 100000, 
           label = paste(sig.meta.count.c[sig.meta.count.c$sig.meta.analysis == 0,]$n, "nonsig."), 
           color = "white")


#Change in significance after meta-analysis, separated by significant effect size estimate
sig.meta <- data.mly %>% 
  select(meta.id, study.id, sig.fixef, sig.ranef, sig.ranef.hkn) %>% 
  gather(key = "meta.analysis.method", value = "sig.meta.analysis", sig.fixef:sig.ranef.hkn) %>% mutate(sig.meta.analysis = factor(sig.meta.analysis))

# duplicate.studies <- data.mly %>% group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr, study.name) %>% 
#   count() %>% filter(n > 1) #Merging is not possible for these since studies bear multiple results, but are not distinguishable and are multiplied while merging.

sig.meta.plot <- merge(y = select(data.ext, meta.id, study.id, sig.type), 
                       x = sig.meta, by = c("meta.id", "study.id")) %>% 
  mutate(sig.type = factor(sig.type))

sig.meta.count <- sig.meta.plot %>% 
  group_by(sig.meta.analysis, meta.analysis.method, sig.type) %>% count()

sig.meta.plot %>% filter(sig.type == 1) %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() + ggtitle("Significant studies") +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("meta.nonsig", "meta.sig")) +
  annotate("text", x = sig.meta.count[which(sig.meta.count$sig.type == 1)[1:3],]$meta.analysis.method, y = 20000, 
           label = paste(sig.meta.count[which(sig.meta.count$sig.type == 1)[1:3],]$n, "changed to nonsig."), 
           color = "white")

sig.meta.plot %>% filter(sig.type == 0) %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() + ggtitle("Non-significant studies") +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("meta.nonsig", "meta.sig")) +
  annotate("text", x = sig.meta.count[which(sig.meta.count$sig.type == 0)[4:6],]$meta.analysis.method, y = 20000, 
           label = paste(sig.meta.count[which(sig.meta.count$sig.type == 0)[4:6],]$n, "changed to sig."), 
           color = "white")


#Change in significance after meta-analysis, separated by significant heterogeneity
sig.meta.Q <- data.mly %>% 
  select(meta.id, sig.Q, sig.fixef, sig.ranef, sig.ranef.hkn) %>% 
  gather(key = "meta.analysis.method", value = "sig.meta.analysis", sig.fixef:sig.ranef.hkn) %>% 
  mutate(sig.meta.analysis = factor(sig.meta.analysis))


sig.meta.count.Q <- sig.meta.Q %>% 
  group_by(sig.meta.analysis, meta.analysis.method, sig.Q) %>% count()

sig.meta.Q %>% filter(sig.Q == 1) %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() + ggtitle("Significant heterogeneity Q studies") +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("meta.nonsig", "meta.sig")) +
  annotate("text", x = sig.meta.count.Q[which(sig.meta.count.Q$sig.Q == 1)[1:3],]$meta.analysis.method, y = 20000, 
           label = paste(sig.meta.count.Q[which(sig.meta.count.Q$sig.Q == 1)[1:3],]$n, "changed to nonsig."), 
           color = "white")

sig.meta.Q %>% filter(sig.Q == 0) %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() + ggtitle("Non-significant heterogeneity Q studies") +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("meta.nonsig", "meta.sig")) +
  annotate("text", x = sig.meta.count.Q[which(sig.meta.count.Q$sig.Q == 0)[1:3],]$meta.analysis.method, y = 100000, 
           label = paste(sig.meta.count.Q[which(sig.meta.count.Q$sig.Q == 0)[1:3],]$n, "changed to nonsig."), 
           color = "white")



######Report Backup:
#Merge comparisons with meta-analysis result:
mly.tomerge <- mly %>% select(meta.id, sig.Q, n.sig.type, n.sig.fishersz, sig.fixef, sig.ranef, sig.ranef.hkn)

data.mly <- merge(x = data.ext, y = mly.tomerge, by = c("meta.id")) #ideally, 278626 comparisons ( = sum(mly$n)), 
# but 300831 (18261 have smaller sample size than mly, but same meta.id)

data.mly <- data.mly %>% filter(total1 > 11 & total2 > 11) #Now 282550



sig.meta.c <- data.mly %>% select(sig.type, sig.fixef, sig.ranef, sig.ranef.hkn) %>% 
  gather(key = "meta.analysis.method", value = "sig.meta.analysis") %>% mutate(sig.meta.analysis = factor(sig.meta.analysis))

sig.meta.count.c <- sig.meta.c %>% 
  group_by(sig.meta.analysis, meta.analysis.method) %>% count()

sig.meta.c <- sig.meta.c %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("meta.nonsig", "meta.sig")) +
  annotate("text", x = sig.meta.count.c[sig.meta.count.c$sig.meta.analysis == 0,]$meta.analysis.method, y = 100000, 
           label = paste(sig.meta.count.c[sig.meta.count.c$sig.meta.analysis == 0,]$n, "nonsig."), 
           color = "white")


#Change in significance after meta-analysis, separated by significant effect size estimate
sig.meta <- data.mly %>% 
  select(meta.id, study.id, sig.fixef, sig.ranef, sig.ranef.hkn) %>% 
  gather(key = "meta.analysis.method", value = "sig.meta.analysis", sig.fixef:sig.ranef.hkn) %>% mutate(sig.meta.analysis = factor(sig.meta.analysis))

# duplicate.studies <- data.mly %>% group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr, study.name) %>% 
#   count() %>% filter(n > 1) #Merging is not possible for these since studies bear multiple results, but are not distinguishable and are multiplied while merging.

sig.meta.plot <- merge(y = select(data.ext, meta.id, study.id, sig.type), 
                       x = sig.meta, by = c("meta.id", "study.id")) %>% 
  mutate(sig.type = factor(sig.type))

sig.meta.count <- sig.meta.plot %>% 
  group_by(sig.meta.analysis, meta.analysis.method, sig.type) %>% count()

sig.meta.plot1 <- sig.meta.plot %>% filter(sig.type == 1) %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() + ggtitle("Significant studies") +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("non-significant", "significant")) +
  annotate("text", x = sig.meta.count[which(sig.meta.count$sig.type == 1)[1:3],]$meta.analysis.method, y = 20000, 
           label = paste(sig.meta.count[which(sig.meta.count$sig.type == 1)[1:3],]$n, "changed to nonsig."), 
           color = "white")

nonsig.meta.plot <- sig.meta.plot %>% filter(sig.type == 0) %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() + ggtitle("Non-significant studies") +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("non-significant", "significant")) +
  annotate("text", x = sig.meta.count[which(sig.meta.count$sig.type == 0)[4:6],]$meta.analysis.method, y = 100000, 
           label = paste(sig.meta.count[which(sig.meta.count$sig.type == 0)[4:6],]$n, "changed to sig."), 
           color = "white")


#Change in significance after meta-analysis, separated by significant heterogeneity
sig.meta.Q <- data.mly %>% 
  select(meta.id, sig.Q, sig.fixef, sig.ranef, sig.ranef.hkn) %>% 
  gather(key = "meta.analysis.method", value = "sig.meta.analysis", sig.fixef:sig.ranef.hkn) %>% 
  mutate(sig.meta.analysis = factor(sig.meta.analysis))


sig.meta.count.Q <- sig.meta.Q %>% 
  group_by(sig.meta.analysis, meta.analysis.method, sig.Q) %>% count()

sig.meta.Q1 <- sig.meta.Q %>% filter(sig.Q == 1) %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() + ggtitle("Significant heterogeneity Q studies") +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("non-significant", "significant")) +
  annotate("text", x = sig.meta.count.Q[which(sig.meta.count.Q$sig.Q == 1)[1:3],]$meta.analysis.method, y = 20000, 
           label = paste(sig.meta.count.Q[which(sig.meta.count.Q$sig.Q == 1)[1:3],]$n, "nonsig."), 
           color = "white")

nonsig.meta.Q <- sig.meta.Q %>% filter(sig.Q == 0) %>%  ggplot(aes(x = meta.analysis.method, fill = sig.meta.analysis)) + 
  geom_bar() + coord_flip() + ggtitle("Non-significant heterogeneity Q studies") +
  theme_bw() + xlab(label = NULL) + scale_fill_discrete(labels= c("non-significant", "significant")) +
  annotate("text", x = sig.meta.count.Q[which(sig.meta.count.Q$sig.Q == 0)[1:3],]$meta.analysis.method, y = 100000, 
           label = paste(sig.meta.count.Q[which(sig.meta.count.Q$sig.Q == 0)[1:3],]$n, "nonsig."), 
           color = "white")

p.secondary.over.meansig <- mly %>% filter(NA.sig.type != 0) %>% ggplot(aes(x = mean.sig.type, fill = factor(sig.ranef), stat(count))) + 
  geom_density(na.rm = T, position = "fill") + 
  labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))

