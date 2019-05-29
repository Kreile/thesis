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





pb.meta.binary = function(data, sig.level, min.study.number) {
  file_fetch = 'pb.bin.RData'
  
  # do not run if file already exists, takes ~60 min 
  if (file.exists(file.path(PATH_RESULTS, file_fetch))) {
    load(file.path(PATH_RESULTS, file_fetch))
  } else {
    meta.bin <- data %>% filter(outcome.measure.new == "Risk Ratio" | outcome.measure.new == "Odds Ratio" | outcome.measure.new == "Risk Difference") %>%  
      filter(file.nr != 3014) %>% #file.nr 3014 gibt eine seltsame Fehlermeldung, muss ausgeschlossen werden
      filter(events1 > 0 | events2 > 0) %>% #comparisons mit events  = 0 ausschliessen
      filter(total1 - events1 > 0 | total2 - events2 > 0) %>% #comparisons mit events = total ausschliessen - vertraegt der alg. irgendwie nicht.
      group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% #Die nachfolgenden Rechnungen werden jeweils fuer diese Gruppen gemacht.
      mutate(n = n()) %>% filter(n > 9) %>% #Nur meta-analysen mit mehr als 9 comparisons behalten
      mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
      mutate(study.year = ifelse(study.year > 1950, study.year, NA)) %>% 
      summarize(doi = unique(doi),
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
                first.publication.year = min(study.year, na.rm = T),
                OR.fixef.bin = #Pooled odds ratio estimate of fixed effects meta-analysis model
                  exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.fixed),
                pval.fixef.bin = #P-value
                  metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.fixed,
                sig.fixef.bin = ifelse(pval.fixef.bin > sig.level, 0, 1), #If significant
                OR.ranef.bin = #Pooled odds ratio estimate of random effects meta-analysis model
                  exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.random),
                pval.ranef.bin = 
                  metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.random,
                sig.ranef.bin = ifelse(pval.ranef.bin > sig.level, 0, 1),
                Q.bin = #Q between studies heterogeneity statistic
                  metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$Q,
                pval.Q.bin = #P-value for between study heterogeneity statistic = 0
                  metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$pval.Q,
                sig.Q.bin = ifelse(pval.Q.bin < sig.level, 1, 0),
                I2  = max(0, (Q.bin - n + 1)/Q.bin), #Proportion of variance of study estimates that is due to real variance
                
                #Publication bias correction:
                lor.trimfill.fixef.bin = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
                  trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE.fixed,
                lor.trimfill.ranef.bin = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
                  trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE.random,
                se.lor.trimfill.ranef.bin = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
                  trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$seTE.random,
                missing.trim.bin = #Estimated number of missing studies in trimfill model
                  trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name))$k0 / n(),
                pval.trimfill.bin = 2*(1-pnorm(lor.trimfill.ranef.bin/se.lor.trimfill.ranef.bin)),
                sig.trimfill.bin = ifelse(pval.trimfill.bin > sig.level, 0, 1),
                lor.reg.ranef.bin = #Adjusted log odds ratio estimate of regression model
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE.adjust,
                se.lor.reg.ranef.bin = #Standard error of adjusted log odds ratio of regression model
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$seTE.adjust,
                pval.reg.ranef.bin = #Pvalue of adjusted log odds ratio estimate of regression model
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$pval.adjust,
                Q.resid = #Residual between study heterogeneity (after adjustment with regression model)
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$Q.resid,
                G.squared = 
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$G.squared,
                sig.reg.ranef.bin = ifelse(pval.reg.ranef.bin > sig.level, 0, 1),
                lor.copas.bin = #Adjusted log odds ratio estimate of copas model
                  auto.copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"), sig.level = sig.level)[1],
                se.lor.copas.bin = #Standard error of adjusted log odds ratio of copas model
                  auto.copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"), sig.level = sig.level)[2],
                missing.copas.bin = #Estimated number of missing studies in copas model
                  auto.copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"), sig.level = sig.level)[3]/n(),
                pval.copas.bin = 2*(1 - pnorm(abs(lor.copas.bin/se.lor.copas.bin))), #Pvalue of adjusted log odds ratio estimate of regression model
                sig.copas.bin = ifelse(pval.trimfill.bin > sig.level, 0, 1),
                sig.change.ranef.reg = sig.change(sig.ranef.bin, sig.after =  sig.reg.ranef.bin), #Function to register if significance of estimate has changed after correction
                sig.change.ranef.copas = sig.change(sig.copas.bin, sig.after =  sig.reg.ranef.bin),
                
                #Publication bias tests:
                pval.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "linreg")$p.val, #P-value (of eggers test)
                egger.test = ifelse(pval.egger.bin < sig.level, 1, 0),
                pval.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                              method = "count")$p.val, #Significance of p-value
                schwarzer.test = ifelse(pval.schwarzer.bin < sig.level, 1, 0),
                pval.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "peter")$p.val,
                peter.test = ifelse(pval.peter.bin < sig.level, 1, 0),
                pval.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                            method = "score")$p.val,
                harbord.test = ifelse(pval.harbord.bin < sig.level, 1, 0),
                pval.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), 
                                           method = "mm")$p.val,
                rucker.test = ifelse(pval.rucker.bin < sig.level, 1, 0),
                stat.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "linreg")$statistic, #Test statistic (of eggers test)
                stat.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                              method = "count")$statistic,
                stat.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "peter")$statistic,
                stat.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                            method = "score")$statistic,
                stat.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), 
                                           method = "mm")$statistic)
    save(meta.bin, file = file.path(PATH_RESULTS, file_fetch))
    load(file.path(PATH_RESULTS, file_fetch))
  }
}



pb.meta.continuous = function(data, sig.level, min.study.number) {
  file_fetch <- 'pb.cont.RData'
  
  # do not run if file already exists, takes ~20 min 
  if (file.exists(file.path(PATH_RESULTS, file_fetch))) {
    load(file.path(PATH_RESULTS, file_fetch))
  } else {
    meta.cont <- data %>% filter(outcome.measure.new == "Mean Difference" | outcome.measure.new == "Std. Mean Difference" | outcome.measure.new == "Hedges' G") %>%
      filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
      filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
      mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
      mutate(study.year = ifelse(study.year > 1950, study.year, NA)) %>% 
      group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% 
      mutate(n = n()) %>% filter(n >= min.study.number) %>% 
      summarize(doi = unique(doi),
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
                first.publication.year = min(study.year, na.rm = T),
                est.fixef.cont = 
                  exp(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$TE.fixed),
                pval.fixef.cont = 
                  metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$pval.fixed,
                sig.fixef.cont = ifelse(pval.fixef.cont > sig.level, 0, 1),
                est.ranef.cont = 
                  exp(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$TE.random),
                pval.ranef.cont = 
                  metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$pval.random,
                sig.ranef.cont = ifelse(pval.ranef.cont > sig.level, 0, 1),
                Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$Q,
                pval.Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$pval.Q,
                sig.Q.cont = ifelse(pval.Q.cont < sig.level, 1, 0),
                I2  = max(0, (Q.cont - n + 1)/Q.cont),
                
                #Pubbias correction:
                est.trimfill.fixef.cont = #Trimfill
                  trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$TE.fixed,
                est.trimfill.ranef.cont = 
                  trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$TE.random,
                se.est.trimfill.ranef.cont = 
                  trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$seTE.random,
                pval.trimfill.cont = 2*(1-pnorm(abs(est.trimfill.ranef.cont/se.est.trimfill.ranef.cont))),
                sig.trimfill.cont = ifelse(pval.trimfill.cont > sig.level, 0, 1),
                missing.trim.cont = 
                  trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$k0 / n(),
                est.reg.ranef.cont = #Regression
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$TE.adjust,
                se.est.reg.ranef.cont = 
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$seTE.adjust,
                pval.reg.ranef.cont = 
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$pval.adjust,
                Q.resid = 
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$Q.resid,
                G.squared = 
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$G.squared,
                sig.reg.ranef.cont = ifelse(pval.reg.ranef.cont > sig.level, 0, 1),
                est.copas.cont = #Copas
                  auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), sig.level = sig.level)[1],
                se.est.copas.cont = 
                  auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), sig.level = sig.level)[2],
                missing.copas.cont = 
                  auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), sig.level = sig.level)[3]/n(),
                pval.copas.cont = 2*(1 - pnorm(abs(est.copas.cont/se.est.copas.cont))),
                sig.copas.cont = ifelse(pval.copas.cont > sig.level, 0, 1),
                sig.change.ranef.reg = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.reg.ranef.cont),
                sig.change.ranef.copas = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.copas.cont),
                
                #Pubbias tests
                pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                           method = "linreg")$p.val,
                egger.test = ifelse(pval.egger.cont < sig.level, 1, 0),
                pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                             method = "mm")$p.val,
                thomson.test = ifelse(pval.thomson.cont < sig.level, 1, 0),
                pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                          method = "rank")$p.val,
                begg.test = ifelse(pval.begg.cont < sig.level, 1, 0),
                stat.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                           method = "linreg")$statistic,
                stat.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                             method = "mm")$statistic,
                stat.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                          method = "rank")$statistic) 
    save(meta.cont, file = file.path(PATH_RESULTS, file_fetch))
    return(meta.cont)
  }
  
  
}

pb.meta.binary = function(data, sig.level, min.study.number) {
  file_fetch = 'pb.bin.RData'
  
  # do not run if file already exists, takes ~60 min 
  if (file.exists(file.path(PATH_RESULTS, file_fetch))) {
    load(file.path(PATH_RESULTS, file_fetch))
  } else {
    meta.bin <- data %>% filter(outcome.measure.new == "Risk Ratio" | outcome.measure.new == "Odds Ratio" | outcome.measure.new == "Risk Difference") %>%  
      filter(file.nr != 3014) %>% #file.nr 3014 gibt eine seltsame Fehlermeldung, muss ausgeschlossen werden
      filter(events1 > 0 | events2 > 0) %>% #comparisons mit events  = 0 ausschliessen
      filter(total1 - events1 > 0 | total2 - events2 > 0) %>% #comparisons mit events = total ausschliessen - vertraegt der alg. irgendwie nicht.
      group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% #Die nachfolgenden Rechnungen werden jeweils fuer diese Gruppen gemacht.
      mutate(n = n()) %>% filter(n > 9) %>% #Nur meta-analysen mit mehr als 9 comparisons behalten
      mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
      mutate(study.year = ifelse(study.year > 1950, study.year, NA)) %>% 
      summarize(doi = unique(doi),
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
                first.publication.year = min(study.year, na.rm = T),
                OR.fixef.bin = #Pooled odds ratio estimate of fixed effects meta-analysis model
                  exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.fixed),
                pval.fixef.bin = #P-value
                  metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.fixed,
                sig.fixef.bin = ifelse(pval.fixef.bin > sig.level, 0, 1), #If significant
                OR.ranef.bin = #Pooled odds ratio estimate of random effects meta-analysis model
                  exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.random),
                pval.ranef.bin = 
                  metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.random,
                sig.ranef.bin = ifelse(pval.ranef.bin > sig.level, 0, 1),
                Q.bin = #Q between studies heterogeneity statistic
                  metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$Q,
                pval.Q.bin = #P-value for between study heterogeneity statistic = 0
                  metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$pval.Q,
                sig.Q.bin = ifelse(pval.Q.bin < sig.level, 1, 0),
                I2  = max(0, (Q.bin - n + 1)/Q.bin), #Proportion of variance of study estimates that is due to real variance
                
                #Publication bias correction:
                lor.trimfill.fixef.bin = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
                  trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE.fixed,
                lor.trimfill.ranef.bin = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
                  trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE.random,
                se.lor.trimfill.ranef.bin = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
                  trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$seTE.random,
                missing.trim.bin = #Estimated number of missing studies in trimfill model
                  trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name))$k0 / n(),
                pval.trimfill.bin = 2*(1-pnorm(lor.trimfill.ranef.bin/se.lor.trimfill.ranef.bin)),
                sig.trimfill.bin = ifelse(pval.trimfill.bin > sig.level, 0, 1),
                lor.reg.ranef.bin = #Adjusted log odds ratio estimate of regression model
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE.adjust,
                se.lor.reg.ranef.bin = #Standard error of adjusted log odds ratio of regression model
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$seTE.adjust,
                pval.reg.ranef.bin = #Pvalue of adjusted log odds ratio estimate of regression model
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$pval.adjust,
                Q.resid = #Residual between study heterogeneity (after adjustment with regression model)
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$Q.resid,
                G.squared = 
                  limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$G.squared,
                sig.reg.ranef.bin = ifelse(pval.reg.ranef.bin > sig.level, 0, 1),
                lor.copas.bin = #Adjusted log odds ratio estimate of copas model
                  auto.copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"), sig.level = sig.level)[1],
                se.lor.copas.bin = #Standard error of adjusted log odds ratio of copas model
                  auto.copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"), sig.level = sig.level)[2],
                missing.copas.bin = #Estimated number of missing studies in copas model
                  auto.copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"), sig.level = sig.level)[3]/n(),
                pval.copas.bin = 2*(1 - pnorm(abs(lor.copas.bin/se.lor.copas.bin))), #Pvalue of adjusted log odds ratio estimate of regression model
                sig.copas.bin = ifelse(pval.trimfill.bin > sig.level, 0, 1),
                sig.change.ranef.reg = sig.change(sig.ranef.bin, sig.after =  sig.reg.ranef.bin), #Function to register if significance of estimate has changed after correction
                sig.change.ranef.copas = sig.change(sig.copas.bin, sig.after =  sig.reg.ranef.bin),
                
                #Publication bias tests:
                pval.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "linreg")$p.val, #P-value (of eggers test)
                egger.test = ifelse(pval.egger.bin < sig.level, 1, 0),
                pval.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                              method = "count")$p.val, #Significance of p-value
                schwarzer.test = ifelse(pval.schwarzer.bin < sig.level, 1, 0),
                pval.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "peter")$p.val,
                peter.test = ifelse(pval.peter.bin < sig.level, 1, 0),
                pval.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                            method = "score")$p.val,
                harbord.test = ifelse(pval.harbord.bin < sig.level, 1, 0),
                pval.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), 
                                           method = "mm")$p.val,
                rucker.test = ifelse(pval.rucker.bin < sig.level, 1, 0),
                stat.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "linreg")$statistic, #Test statistic (of eggers test)
                stat.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                              method = "count")$statistic,
                stat.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "peter")$statistic,
                stat.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                            method = "score")$statistic,
                stat.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), 
                                           method = "mm")$statistic)
    save(meta.bin, file = file.path(PATH_RESULTS, file_fetch))
    load(file.path(PATH_RESULTS, file_fetch))
  }
}



pb.meta.continuous = function(data, sig.level, min.study.number) {
  file_fetch <- 'pb.cont.RData'
  
  # do not run if file already exists, takes ~20 min 
  if (file.exists(file.path(PATH_RESULTS, file_fetch))) {
    load(file.path(PATH_RESULTS, file_fetch))
  } else {
    meta.cont <- data %>% filter(outcome.measure.new == "Mean Difference" | outcome.measure.new == "Std. Mean Difference" | outcome.measure.new == "Hedges' G") %>%
      filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
      filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
      mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
      mutate(study.year = ifelse(study.year > 1950, study.year, NA)) %>% 
      group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% 
      mutate(n = n()) %>% filter(n >= min.study.number) %>% 
      summarize(doi = unique(doi),
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
                first.publication.year = min(study.year, na.rm = T),
                est.fixef.cont = 
                  exp(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$TE.fixed),
                pval.fixef.cont = 
                  metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$pval.fixed,
                sig.fixef.cont = ifelse(pval.fixef.cont > sig.level, 0, 1),
                est.ranef.cont = 
                  exp(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$TE.random),
                pval.ranef.cont = 
                  metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$pval.random,
                sig.ranef.cont = ifelse(pval.ranef.cont > sig.level, 0, 1),
                Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$Q,
                pval.Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$pval.Q,
                sig.Q.cont = ifelse(pval.Q.cont < sig.level, 1, 0),
                I2  = max(0, (Q.cont - n + 1)/Q.cont),
                
                #Pubbias correction:
                est.trimfill.fixef.cont = #Trimfill
                  trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$TE.fixed,
                est.trimfill.ranef.cont = 
                  trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$TE.random,
                se.est.trimfill.ranef.cont = 
                  trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$seTE.random,
                pval.trimfill.cont = 2*(1-pnorm(abs(est.trimfill.ranef.cont/se.est.trimfill.ranef.cont))),
                sig.trimfill.cont = ifelse(pval.trimfill.cont > sig.level, 0, 1),
                missing.trim.cont = 
                  trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$k0 / n(),
                est.reg.ranef.cont = #Regression
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$TE.adjust,
                se.est.reg.ranef.cont = 
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$seTE.adjust,
                pval.reg.ranef.cont = 
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$pval.adjust,
                Q.resid = 
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$Q.resid,
                G.squared = 
                  limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$G.squared,
                sig.reg.ranef.cont = ifelse(pval.reg.ranef.cont > sig.level, 0, 1),
                est.copas.cont = #Copas
                  auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), sig.level = sig.level)[1],
                se.est.copas.cont = 
                  auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), sig.level = sig.level)[2],
                missing.copas.cont = 
                  auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), sig.level = sig.level)[3]/n(),
                pval.copas.cont = 2*(1 - pnorm(abs(est.copas.cont/se.est.copas.cont))),
                sig.copas.cont = ifelse(pval.copas.cont > sig.level, 0, 1),
                sig.change.ranef.reg = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.reg.ranef.cont),
                sig.change.ranef.copas = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.copas.cont),
                
                #Pubbias tests
                pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                           method = "linreg")$p.val,
                egger.test = ifelse(pval.egger.cont < sig.level, 1, 0),
                pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                             method = "mm")$p.val,
                thomson.test = ifelse(pval.thomson.cont < sig.level, 1, 0),
                pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                          method = "rank")$p.val,
                begg.test = ifelse(pval.begg.cont < sig.level, 1, 0),
                stat.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                           method = "linreg")$statistic,
                stat.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                             method = "mm")$statistic,
                stat.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                          method = "rank")$statistic) 
    save(meta.cont, file = file.path(PATH_RESULTS, file_fetch))
    return(meta.cont)
  }
  
  
}







# sig.change <- function(sig.before, sig.after){
#   if(!is.na(sig.before & !is.na(sig.after))){
#     if(sig.before == 0){
#       results <- ifelse(sig.after == 0, "unchanged.nonsig", "change.to.sig")
#     } else{
#       results <- ifelse(sig.after == 1, "unchanged.sig", "change.to.nonsig")}
#   } else{
#     results <- NA
#   }
#   return(results)
# }


meta.cont %>% ggplot(aes(y = pval.egger.cont, x = pval.thomson.cont)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta.cont %>% ggplot(aes(x = pval.thomson.cont, y = pval.begg.cont)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta.cont %>% filter(n < 100) %>% ggplot(aes(y = pval.thomson.cont, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess", se = F) 
meta.cont %>% filter(n < 100) %>% ggplot(aes(y = pval.egger.cont, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess", se = F) 
meta.cont %>% ggplot(aes(y = stat.egger.cont, x = stat.thomson.cont)) + geom_point(alpha = 0.3) + theme_bw() #+ geom_smooth(method = "lm")  
meta.cont %>% ggplot(aes(x = pval.thomson.cont, y = missing.trim.cont)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 

meta.bin %>% ggplot(aes(x = pval.harbord.bin, y = pval.peter.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta.bin %>% filter(n < 150)%>% ggplot(aes(y = pval.peter.bin, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta.bin %>% filter(n < 150)%>% ggplot(aes(y = pval.harbord.bin, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(y = pval.egger.bin, x = pval.peter.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(y = pval.egger.bin, x = pval.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(x = pval.schwarzer.bin, y = pval.peter.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(x = pval.schwarzer.bin, y = pval.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(y = pval.rucker.bin, y = pval.egger.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(y = pval.rucker.bin, y = pval.peter.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(y = pval.rucker.bin, y = pval.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 

meta.bin %>% ggplot(aes(x = stat.peter.bin, y = stat.schwarzer.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(x = stat.peter.bin, y = stat.egger.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(x = stat.peter.bin, y = stat.rucker.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(x = stat.peter.bin, y = stat.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess")

meta.bin %>% ggplot(aes(y = stat.egger.bin, x = stat.schwarzer.bin)) + geom_point(alpha = 0.5) + theme_bw()# + geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(y = stat.rucker.bin, x = stat.schwarzer.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(y = stat.harbord.bin, x = stat.schwarzer.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 

meta.bin %>% ggplot(aes(x = stat.rucker.bin, y = stat.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta.bin %>% ggplot(aes(x = stat.rucker.bin, y = stat.egger.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 

meta.bin %>% ggplot(aes(x = stat.harbord.bin, y = stat.egger.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess")
















# meta.bin.plot <- meta.bin %>% ungroup() %>% select(pval.peter.bin, pval.harbord.bin, trim.bin, Q, I2)
# meta.cont.plot <- meta.cont %>% ungroup() %>% select(pval.egger.cont, pval.thomson.cont, pval.begg.cont, missing.missing.trim.cont, Q, I2)
# scatterplotMatrix(meta.bin.plot)
# scatterplotMatrix(meta.cont.plot)
# pairs(meta.bin.plot)

# data %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% #filter(file.nr < 503) %>% 
#   filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
#   filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
#   group_by(file.nr, outcome.nr, subgroup.nr) %>% 
#   mutate(n = n()) %>% filter(n > 9) %>% 
#   summarize(pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
#                                          method = "mm")$p.val,
#             Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$Q, n = n(),
#             I2  = max(0, (Q - n + 1)/Q)) %>% 
#   ggplot(aes(x = I2, y = pval.thomson.cont)) + geom_point(alpha = 0.7) + theme_bw() + geom_smooth(method = "lm")
# 
# 
# 
# data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% filter(file.nr != 3014) %>% 
#   filter(events1 > 0 | events2 > 0) %>% 
#   filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
#   group_by(file.nr, outcome.nr, subgroup.nr) %>% 
#   mutate(n = n()) %>% filter(n > 9) %>% 
#   summarize(pval.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), method = "peter")$p.val,
#             Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$Q, n = n(),
#             I2  = max(0, (Q - n + 1)/Q)) %>%
#   ggplot(aes(x = I2, y = pval.peter.bin)) + geom_point(alpha = 0.7) + theme_bw() + geom_smooth(method = "lm")

# #Proportion of pbias according to 5% significance level:
# rejection.bin <- meta.bin %>% summarize(egger.rejection = mean(egger.test),
#                                             schwarzer.rejection = mean(schwarzer.test),
#                                             rucker.rejection = mean(rucker.test),
#                                             harbord.rejection = mean(harbord.test),
#                                             peter.rejection = mean(peter.test))
# rejection.bin <- rejection.bin %>% gather(key = test.type)
# rejection.bin %>% ggplot(aes(y = value, x = test.type)) + geom_col()
# bin.spread <- meta.bin %>% select(egger.test, schwarzer.test, harbord.test, peter.test) %>% gather(key = test.type, value = "null.hypothesis")
# bin.spread$value <- ifelse(bin.spread$value == 1, "reject", "accept")
# bin.spread$value <- factor(bin.spread$value)
# bin.spread %>% ggplot(aes(x = test.type, fill = value)) + geom_bar() + coord_flip()
# 
# rejection.cont <- meta.cont %>% summarize(egger.rejection = mean(egger.test),
#                                               begg.rejection = mean(begg.test),
#                                               thomson.rejection = mean(thomson.test))
# 
# cont.spread <- meta.cont %>% select(egger.test, begg.test, thomson.test) %>% gather(key = test.type, value = "null.hypothesis")
# cont.spread$value <- ifelse(cont.spread$value == 1, "rejected", "accepted")
# cont.spread$value <- factor(cont.spread$value)
# cont.spread %>% ggplot(aes(x = test.type, fill = value)) + geom_bar() + coord_flip() + theme_bw()
# 









pb.bias.cont.extended <- function(data, min.study.number, sig.level){
  metadat <- data %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>%# filter(file.nr < 503) %>% 
    filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
    filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), n = n(),
              pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                         method = "linreg")$p.val,
              egger.test = ifelse(pval.egger.cont < sig.level, 1, 0),
              pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                           method = "mm")$p.val,
              thomson.test = ifelse(pval.thomson.cont < sig.level, 1, 0),
              pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                        method = "rank")$p.val,
              begg.test = ifelse(pval.begg.cont < sig.level, 1, 0),
              stat.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                         method = "linreg")$statistic,
              stat.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                           method = "mm")$statistic,
              stat.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                        method = "rank")$statistic,
              missing.trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$k0 / n(),
              Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() %>% select(-outcome.nr, -subgroup.nr)
  return(metadat)
}


pb.bias.bin.extended <- function(data, min.study.number, sig.level){
  metadat <- data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% filter(file.nr != 3014) %>% 
    filter(events1 > 0 | events2 > 0) %>% 
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), n = n(),
              pval.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), method = "linreg")$p.val,
              egger.test = ifelse(pval.egger.bin < sig.level, 1, 0),
              pval.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), method = "count")$p.val,
              schwarzer.test = ifelse(pval.schwarzer.bin < sig.level, 1, 0),
              pval.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), method = "peter")$p.val,
              peter.test = ifelse(pval.peter.bin < sig.level, 1, 0),
              pval.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), method = "score")$p.val,
              harbord.test = ifelse(pval.harbord.bin < sig.level, 1, 0),
              pval.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), method = "mm")$p.val,
              rucker.test = ifelse(pval.rucker.bin < sig.level, 1, 0),
              stat.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), method = "linreg")$statistic,
              stat.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), method = "count")$statistic,
              stat.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), method = "peter")$statistic,
              stat.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), method = "score")$statistic,
              stat.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), method = "mm")$statistic,
              trim.bin = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name))$k0 / n(),
              Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() %>% select(-outcome.nr, -subgroup.nr)
  return(metadat)
} 


pb.bias.cont.extended <- function(data, min.study.number, sig.level){
  metadat <- data %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>%# filter(file.nr < 503) %>% 
    filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
    filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), n = n(),
              pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                         method = "linreg")$p.val,
              egger.test = ifelse(pval.egger.cont < sig.level, 1, 0),
              pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                           method = "mm")$p.val,
              thomson.test = ifelse(pval.thomson.cont < sig.level, 1, 0),
              pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                        method = "rank")$p.val,
              begg.test = ifelse(pval.begg.cont < sig.level, 1, 0),
              stat.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                         method = "linreg")$statistic,
              stat.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                           method = "mm")$statistic,
              stat.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
                                        method = "rank")$statistic,
              missing.trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$k0 / n(),
              Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() %>% select(-outcome.nr, -subgroup.nr)
  return(metadat)
}



pb.corr.cont <- function(data, min.study.number){
  metadat <- data %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>%# filter(file.nr < 503) %>% 
    filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
    filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), n = n(),
              # pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
              #                            method = "linreg")$p.val,
              # pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
              #                              method = "mm")$p.val,
              # pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name), 
              #                           method = "rank")$p.val,
              # missing.trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$k0 / n(),
              Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() %>% select(-outcome.nr, -subgroup.nr)
  return(metadat)
}

pb.corr.bin <- function(data, min.study.number){
  metadat <- data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio"| outcome.measure == "Risk Difference") %>% filter(file.nr != 3014) %>% 
    filter(events1 > 0 | events2 > 0) %>% 
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), n = n(),
              LOR.fixef.bin = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.fixed,
              pval.fixef.bin = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.fixed,
              LOR.ranef.bin = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.random,
              pval.fixef.bin = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.random,
              trimfill.bin = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE,
              copas.which = which.min(abs(0.05 - copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$pval.rsb
                                          [copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$pval.rsb > 0.05])),
              copas.max.bin = copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE[copas.which],
              reg.bin = limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE,
              Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() 
  return(metadat)
}  



pb.corr.bin <- function(data, min.study.number){
  metadat <- data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio" | outcome.measure == "") %>% filter(file.nr != 3014) %>% 
    filter(events1 > 0 | events2 > 0) %>% 
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), n = n(),
              copas.which = which.min(abs(0.05 - copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$pval.rsb
                                          [copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$pval.rsb > 0.05])),
              copas.bin = copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE[copas.which]) %>% 
    ungroup() 
  return(metadat)
}  

