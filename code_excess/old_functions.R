
meta.bin.complete = function(data, sig.level, min.study.number, sm1) {
	meta.bin <- data %>% filter(outcome.measure.new == "Risk Ratio" | outcome.measure.new == "Odds Ratio" | outcome.measure.new == "Risk Difference") %>%  
		filter(file.nr != 3014 & file.nr != 208) %>% #file.nr 3014 gibt eine seltsame Fehlermeldung, muss ausgeschlossen werden
		filter(events1 > 0 | events2 > 0) %>% #comparisons mit events  = 0 ausschliessen
		filter(total1 - events1 > 0 | total2 - events2 > 0) %>% #comparisons mit events = total ausschliessen - vertraegt der alg. irgendwie nicht.
		group_by(meta.id) %>% #Die nachfolgenden Rechnungen werden jeweils fuer diese Gruppen gemacht.
		mutate(n = n()) %>% filter(n >= min.study.number) %>% #Nur meta-analysen mit mehr als 9 comparisons behalten
		mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
		mutate(study.year = ifelse(study.year > 1950, study.year, NA)) %>% 
		summarize(doi = unique(doi),
							file.nr = unique(file.nr),
							comparison.nr = unique(comparison.nr),
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
							n.sig.type.bin = sum(sig.type, na.rm = T),
							NA.sig.type = sum(is.na(sig.type)),
							se.min = min(se, na.rm = T),
							se.max = max(se, na.rm = T),
							
							#Meta-analysis:
							est.fixef.bin = #Pooled odds ratio estimate of fixed effects meta-analysis model
								exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1,
														method = "Inverse")$TE.fixed),
							se.est.fixef.bin = #Pooled odds ratio estimate of fixed effects meta-analysis model
								exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1,
														method = "Inverse")$seTE.fixed),
							zval.fixef.bin = #P-value
								metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1,
												method = "Inverse")$zval.fixed,
							pval.fixef.bin = #P-value
								metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1,
												method = "Inverse")$pval.fixed,
							sig.fixef.bin = ifelse(pval.fixef.bin > sig.level, 0, 1), #If significant
							se.est.ranef.bin = #Pooled odds ratio estimate of random effects meta-analysis model
								exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1)$seTE.random),
							est.ranef.bin = #Pooled odds ratio estimate of random effects meta-analysis model
								exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1)$TE.random),
							zval.ranef.bin = 
								metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1)$zval.random,
							pval.ranef.bin = 
								metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1)$pval.random,
							sig.ranef.bin = ifelse(pval.ranef.bin > sig.level, 0, 1),
							Q.bin = #Q between studies heterogeneity statistic
								metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1)$Q,
							pval.Q.bin = #P-value for between study heterogeneity statistic = 0
								metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1)$pval.Q,
							sig.Q.bin = ifelse(pval.Q.bin < sig.level, 1, 0),
							I2  = max(0, (Q.bin - n + 1)/Q.bin), #Proportion of variance of study estimates that is due to real variance
							
							#Publication bias correction:
							logest.trimfill.fixef.bin = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
								(trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1,
																	method = "Inverse"))$TE.fixed),
							logest.trimfill.ranef.bin = #Adjusted pooled log odds ratio estimate of fixed effects meta-analysis model of trimfill model
								(trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1,
																	method = "Inverse"))$TE.random),
							se.logest.trimfill.ranef.bin = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1,
																 method = "Inverse"))$seTE.random,
							se.logest.trimfill.fixef.bin = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1,
																 method = "Inverse"))$seTE.fixed,
							zval.trimfill.fixef = #Standard error of adjusted log odds ratio of random effects meta-analysis of trimfill model
								trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1,
																 method = "Inverse"))$zval.fixed,
							missing.trim.bin = #Estimated number of missing studies in trimfill model
								trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1,
																 method = "Inverse"))$k0 / n(),
							pval.trimfill.bin = 2*(1-pnorm(logest.trimfill.ranef.bin/se.logest.trimfill.ranef.bin)),
							sig.trimfill.bin = ifelse(pval.trimfill.bin > sig.level, 0, 1),
							logest.reg.ranef.bin = #Adjusted log odds ratio estimate of regression model
								limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1))$TE.adjust,
							se.logest.reg.ranef.bin = #Standard error of adjusted log odds ratio of regression model
								limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1))$seTE.adjust,
							pval.reg.ranef.bin = #Pvalue of adjusted log odds ratio estimate of regression model
								limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1))$pval.adjust,
							Q.resid = #Residual between study heterogeneity (after adjustment with regression model)
								limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1))$Q.resid,
							G.squared = 
								limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1))$G.squared,
							sig.reg.ranef.bin = ifelse(pval.reg.ranef.bin > sig.level, 0, 1),
							logest.copas.bin = #Adjusted log odds ratio estimate of copas model
								auto.copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1), sig.level = sig.level)[1],
							se.logest.copas.bin = #Standard error of adjusted log odds ratio of copas model
								auto.copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1), sig.level = sig.level)[2],
							missing.copas.bin = #Estimated number of missing studies in copas model
								auto.copas(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= sm1), sig.level = sig.level)[3]/n(),
							pval.copas.bin = 2*(1 - pnorm(abs(logest.copas.bin/se.logest.copas.bin))), #Pvalue of adjusted log odds ratio estimate of regression model
							sig.copas.bin = ifelse(pval.copas.bin > sig.level, 0, 1),
							
							#Significance change
							sig.change.ranef.reg = sig.change(sig.before = sig.ranef.bin, sig.after =  sig.reg.ranef.bin), #Function to register if significance of estimate has changed after correction
							sig.change.ranef.copas = sig.change(sig.before = sig.ranef.bin, sig.after =  sig.copas.bin),
							sig.change.ranef.trimfill = sig.change(sig.before = sig.ranef.bin, sig.after =  sig.trimfill.bin),
							sig.change.fixef.reg = sig.change(sig.before = sig.fixef.bin, sig.after =  sig.reg.ranef.bin),
							sig.change.fixef.copas = sig.change(sig.before = sig.fixef.bin, sig.after =  sig.copas.bin),
							sig.change.fixef.trimfill = sig.change(sig.before = sig.fixef.bin, sig.after =  sig.trimfill.bin),
							
							#Publication bias tests:
							pval.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1), 
																				method.bias ="linreg")$p.val, #P-value (of eggers test)
							egger.test = ifelse(pval.egger.bin < sig.level, 1, 0),
							pval.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1), 
																						method.bias ="count")$p.val, #Significance of p-value
							schwarzer.test = ifelse(pval.schwarzer.bin < sig.level, 1, 0),
							pval.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1), 
																				method.bias ="peter")$p.val,
							peter.test = ifelse(pval.peter.bin < sig.level, 1, 0),
							pval.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1), 
																					method.bias ="score")$p.val,
							harbord.test = ifelse(pval.harbord.bin < sig.level, 1, 0),
							pval.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), 
																				 method.bias ="mm")$p.val,
							rucker.test = ifelse(pval.rucker.bin < sig.level, 1, 0),
							stat.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1), 
																				method.bias ="linreg")$statistic, #Test statistic (of eggers test)
							stat.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1), 
																						method.bias ="count")$statistic,
							stat.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1), 
																				method.bias ="peter")$statistic,
							stat.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = sm1), 
																					method.bias ="score")$statistic,
							stat.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), 
																				 method.bias ="mm")$statistic)
}



meta.cont.complete = function(data, sig.level, min.study.number) {
	meta.cont <- data %>% filter(outcome.measure.new == "Mean Difference" | outcome.measure.new == "Std. Mean Difference") %>%
		filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
		filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
		mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
		mutate(study.year = ifelse(study.year > 1950, study.year, NA)) %>% 
		group_by(meta.id) %>% 
		mutate(n = n()) %>% filter(n >= min.study.number) %>% 
		summarize(doi = unique(doi),
							file.nr = unique(file.nr),
							comparison.nr = unique(comparison.nr),
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
							n.sig.type.cont = sum(sig.type, na.rm = T),
							NA.sig.type = sum(is.na(sig.type)),
							se.min = min(se, na.rm = T),
							se.max = max(se, na.rm = T),
							
							#Meta-analysis:
							est.fixef.cont = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$TE.fixed,
							se.est.fixef.cont = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$seTE.fixed,
							zval.fixef.cont = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$zval.fixed,
							pval.fixef.cont = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$pval.fixed,
							sig.fixef.cont = ifelse(pval.fixef.cont > sig.level, 0, 1),
							est.ranef.cont = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$TE.random,
							se.est.ranef.cont = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$seTE.random,
							zval.ranef.cont = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$zval.random,
							pval.ranef.cont = 
								metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$pval.random,
							sig.ranef.cont = ifelse(pval.ranef.cont > sig.level, 0, 1),
							Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$Q,
							pval.Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name)$pval.Q,
							sig.Q.cont = ifelse(pval.Q.cont < sig.level, 1, 0),
							I2  = max(0, (Q.cont - n + 1)/Q.cont),
							
							#Pubbias correction:
							est.trimfill.fixef.cont = #Trimfill
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$TE.fixed,
							est.trimfill.ranef.cont = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$TE.random,
							se.est.trimfill.ranef.cont = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$seTE.random,
							se.est.trimfill.fixef.cont = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$seTE.fixed,
							zval.trimfill.fixef = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$zval.fixed,
							pval.trimfill.cont = 2*(1-pnorm(abs(est.trimfill.ranef.cont/se.est.trimfill.ranef.cont))),
							sig.trimfill.cont = ifelse(pval.trimfill.cont > sig.level, 0, 1),
							missing.trim.cont = 
								trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$k0 / n(),
							est.reg.ranef.cont = #Regression
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$TE.adjust,
							se.est.reg.ranef.cont = 
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$seTE.adjust,
							pval.reg.ranef.cont = 
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$pval.adjust,
							Q.resid = 
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$Q.resid,
							G.squared = 
								limitmeta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name))$G.squared,
							sig.reg.ranef.cont = ifelse(pval.reg.ranef.cont > sig.level, 0, 1),
							est.copas.cont = #Copas
								auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), sig.level = sig.level)[1],
							se.est.copas.cont = 
								auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), sig.level = sig.level)[2],
							missing.copas.cont = 
								auto.copas(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), sig.level = sig.level)[3]/n(),
							pval.copas.cont = 2*(1 - pnorm(abs(est.copas.cont/se.est.copas.cont))),
							sig.copas.cont = ifelse(pval.copas.cont > sig.level, 0, 1),
							
							#Significance change
							sig.change.ranef.reg = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.reg.ranef.cont),
							sig.change.ranef.copas = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.copas.cont),
							sig.change.ranef.trimfill = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.trimfill.cont),
							sig.change.fixef.reg = sig.change(sig.before = sig.fixef.cont, sig.after =  sig.reg.ranef.cont),
							sig.change.fixef.copas = sig.change(sig.before = sig.fixef.cont, sig.after =  sig.copas.cont),
							sig.change.fixef.trimfill = sig.change(sig.before = sig.fixef.cont, sig.after =  sig.trimfill.cont),
							
							#Pubbias tests
							pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																				 method.bias ="linreg")$p.val,
							egger.test = ifelse(pval.egger.cont < sig.level, 1, 0),
							pval.thompson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																						method.bias ="mm")$p.val,
							thompson.test = ifelse(pval.thompson.cont < sig.level, 1, 0),
							pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																				method.bias ="rank")$p.val,
							begg.test = ifelse(pval.begg.cont < sig.level, 1, 0),
							stat.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																				 method.bias ="linreg")$statistic,
							stat.thompson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																						method.bias ="mm")$statistic,
							stat.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																				method.bias ="rank")$statistic) 
	
	return(meta.cont)
}

pb.meta.merge <- function(meta.binary, meta.continous){
	meta <- bind_rows(meta.bin, meta.cont) %>%
		mutate(n.sig.type = ifelse(!is.na(n.sig.type.bin), n.sig.type.bin, n.sig.type.cont),
					 pval.thompson = ifelse(!is.na(pval.rucker.bin), pval.rucker.bin, pval.thompson.cont),
					 thompson.test = ifelse(!is.na(thompson.test), thompson.test, rucker.test),
					 
					 est.fixef = ifelse(!is.na(est.fixef.bin), est.fixef.bin, est.fixef.cont),
					 zval.fixef = ifelse(!is.na(zval.fixef.bin), zval.fixef.bin, zval.fixef.cont),
					 
					 est.ranef = ifelse(!is.na(est.ranef.bin), est.ranef.bin, est.ranef.cont),
					 zval.ranef = ifelse(!is.na(zval.ranef.bin), zval.ranef.bin, zval.ranef.cont),
					 
					 est.trimfill.fixef = ifelse(!is.na(logest.trimfill.fixef.bin), exp(logest.trimfill.fixef.bin), est.trimfill.fixef.cont),
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


#Cochrane data process function, constructs addtitional variables:
#meta.id: unique ID for each (possible) meta-analysis.
#outcome.type: classifies all data into binary ("bin"), continuous ("cont") and survival ("surv") outcome types
#pval.type: p-value of a test against the null hypothesis of no treatment effect
#pval.fishersz: p-value of a test for the fishers z correlation to be equal to 0 (not working well)
#sig.: significance of p-value based on 0.05 threshold
#larger12.ss: sample size of both treatment arms is larger than 12
pb.process <- function(data){
	data <- data %>% mutate(meta.id = group_indices(., file.nr, comparison.nr, outcome.nr, subgroup.nr)) %>%
		group_by(meta.id) %>% mutate(study.id = row_number())
	
	data <- data %>% mutate(odds.ratio = NA,
													std.mean.d = NA,
													correlation = NA,
													fishersz = NA, 
													fishersz.variance = NA,
													log.OR = NA,
													se.log.OR = NA,
													std.m.d = NA,
													se.std.m.d = NA, 
													outcome.type = ifelse(outcome.measure.new == "Hazard Ratio", "surv", 
																								ifelse(outcome.measure.new == "Odds Ratio" | outcome.measure.new == "Risk Ratio" |
																											 	outcome.measure.new == "Peto Odds Ratio" | outcome.measure.new == "Risk Difference", "bin", "cont")),
													pval.type = NA,
													pval.fishersz = NA,
													sig.type = NA,
													sig.fishersz = NA,
													larger12.ss = ifelse(total1 > 12 & total2 > 12, 1, 0))
	
	binary.unprob <- which(data$events1 > 0 & data$events2 > 0 & data$total1 - data$events1 > 0 & data$total2 - data$events2 > 0 &
												 	data$total1 + data$total2 > 3 & data$larger12.ss == 1& data$outcome.type == "bin")
	
	data[binary.unprob,] <- data[binary.unprob,] %>%
		mutate(log.OR = log( (events1 * (total2 - events2))/ (events2 * (total1 - events1)) ), 
					 se.log.OR = sqrt( (1/events1) + 1/(total1 - events1) + (1/events2) + 1/(total2 - events2)),
					 std.mean.d = log.OR * (sqrt(3)/pi),
					 correlation = std.mean.d/sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1- pnorm(abs(log.OR/se.log.OR))))
	
	#Indices of comparisons with zero events or total events
	zero1 <- which(data$events2 == 0 & data$total1 - data$events1 > 0 & data$outcome.type == "bin" & data$larger12.ss == 1)
	zero2 <- which(data$events1 == 0 & data$total2 - data$events2 > 0 & data$outcome.type == "bin" & data$larger12.ss == 1)
	binary.zeros <- c(zero1, zero2)
	totals1 <- which(data$larger12.ss == 1 & data$total1 - data$events1 == 0 & data$events2 > 0 & data$outcome.type == "bin")
	totals2 <- which(data$larger12.ss == 1 & data$total2 - data$events2 == 0 & data$events1 > 0 & data$outcome.type == "bin")
	binary.totals <- c(totals1, totals2)
	
	data[binary.zeros,] <- data[binary.zeros,] %>%
		mutate(events1 = events1 + 0.5, events2 = events2 + 0.5,
					 log.OR = log( (events1 * (total2 - events2))/ (events2 * (total1 - events1)) ), 
					 se.log.OR = sqrt( (1/events1) + 1/(total1 - events1) + (1/events2) + 1/(total2 - events2)),
					 std.mean.d = log.OR * (sqrt(3)/pi),
					 correlation = std.mean.d/sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1- pnorm(abs(log.OR/se.log.OR))))
	
	data[binary.totals,] <- data[binary.totals,] %>%
		mutate(events1 = events1 - 0.5, events2 = events2 - 0.5,
					 log.OR = log( (events1 * (total2 - events2))/ (events2 * (total1 - events1)) ), 
					 se.log.OR = sqrt( (1/events1) + 1/(total1 - events1) + (1/events2) + 1/(total2 - events2)),
					 std.mean.d = log.OR * (sqrt(3)/pi),
					 correlation = std.mean.d/sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1- pnorm(abs(log.OR/se.log.OR))))
	
	#from 307311 outcomes with outcome type "bin", 280886 have been used (difference = 26425). 26370 have too low sample size
	
	data[!is.na(data$outcome.measure.new) & data$outcome.measure.new == "Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Mean Difference") %>%
		mutate(std.m.d = (mean1 - mean2)/sqrt((sd1^2+sd2^2)/2),
					 se.std.m.d = sqrt( (total1 + total2)/(total1 * total2) + (std.m.d^2)/(2*(total1 + total2)) ),
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	data[!is.na(data$outcome.measure.new) & data$outcome.measure.new == "Std. Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Std. Mean Difference") %>%
		mutate(std.m.d = (mean1 - mean2)/sqrt((sd1^2+sd2^2)/2), #No correction factor used..
					 se.std.m.d = sqrt( (total1 + total2)/(total1 * total2) + (std.m.d^2)/(2*(total1 + total2)) ),
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	#Calculate pvalues etc. when no mean or sd values are given
	data[is.na(data$mean1) & !is.na(data$outcome.measure.new) & data$outcome.measure.new == "Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Mean Difference", is.na(mean1)) %>%
		mutate(std.m.d = effect/sqrt((sd1^2+sd2^2)/2),
					 se.std.m.d = sqrt( (total1 + total2)/(total1 * total2) + (std.m.d^2)/(2*(total1 + total2)) ),
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	data[is.na(data$mean2) & !is.na(data$outcome.measure.new) & data$outcome.measure.new == "Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Mean Difference", is.na(mean2)) %>%
		mutate(std.m.d = effect/sqrt((sd1^2+sd2^2)/2),
					 se.std.m.d = sqrt( (total1 + total2)/(total1 * total2) + (std.m.d^2)/(2*(total1 + total2)) ),
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	data[is.na(data$sd1) & !is.na(data$outcome.measure.new) & data$outcome.measure.new == "Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Mean Difference", is.na(sd1)) %>%
		mutate(std.m.d = effect/se,
					 se.std.m.d = sqrt( (total1 + total2)/(total1 * total2) + (std.m.d^2)/(2*(total1 + total2)) ),
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	data[is.na(data$sd2) & !is.na(data$outcome.measure.new) & data$outcome.measure.new == "Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Mean Difference", is.na(sd2)) %>%
		mutate(std.m.d = effect/se,
					 se.std.m.d = sqrt( (total1 + total2)/(total1 * total2) + (std.m.d^2)/(2*(total1 + total2)) ),
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	#Calculate pvalues etc. when no mean or sd values are given
	data[is.na(data$mean1) & !is.na(data$outcome.measure.new) & data$outcome.measure.new == "Std. Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Std. Mean Difference", is.na(mean1)) %>%
		mutate(std.m.d = effect,
					 se.std.m.d = se,
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	data[is.na(data$mean2) & !is.na(data$outcome.measure.new) & data$outcome.measure.new == "Std. Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Std. Mean Difference", is.na(mean2)) %>%
		mutate(std.m.d = effect,
					 se.std.m.d = se,
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	data[is.na(data$sd2) & !is.na(data$outcome.measure.new) & data$outcome.measure.new == "Std. Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Std. Mean Difference", is.na(sd2)) %>%
		mutate(std.m.d = effect,
					 se.std.m.d = se,
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
					 fishersz.variance = 1/(total1 + total2 - 3),
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	data[is.na(data$sd1) & !is.na(data$outcome.measure.new) & data$outcome.measure.new == "Std. Mean Difference", ] <- data %>% 
		filter(outcome.measure.new == "Std. Mean Difference", is.na(sd1)) %>%
		mutate(std.m.d = effect,
					 se.std.m.d = se,
					 correlation = std.m.d/sqrt(std.m.d^2 + ((total1 + total2)^2)/(total2*total1)),
					 
					 pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))))
	
	#142850 have mean or std mean difference, 126620 can be used for calculations (92781 from 102315 md (have small ss or missing sds and effects))
	
	data[data$outcome.type == "surv" & !is.na(data$se) & !is.na(data$outcome.measure.new), ] <- data %>% filter(outcome.type == "surv" & !is.na(se)) %>% 
		mutate(pval.type = 2*(1-pnorm(abs((effect - 1) / se))))
	
	
	data[!is.na(data$fishersz) & !is.na(data$fishersz.variance),] <- data[!is.na(data$fishersz) & !is.na(data$fishersz.variance),] %>% 
		mutate(pval.fishersz = 2*(1-pnorm(abs(fishersz/fishersz.variance))))
	
	data[!is.na(data$pval.type),] <- data[!is.na(data$pval.type),] %>% 
		mutate(sig.type = ifelse(pval.type < 0.05, 1, 0))
	
	data[!is.na(data$pval.fishersz),] <- data[!is.na(data$pval.fishersz),] %>% 
		mutate(sig.fishersz = ifelse(pval.fishersz < 0.05, 1, 0))
	
	return(data)
}




mcor <- function(data){
	if(all(data$outcome.type == "bin")){
		mt <- metacor(cor = cor.phi, n = total1 + total2, studlab = study.name, data = data)
	}else{
		if(all(data$outcome.type == "cont")){
			mt <- metacor(cor = cor.pearson, n = total1 + total2, studlab = study.name, data = data)
		} else{
			if(all(data$outcome.type == "surv")){
				mt <- metagen(TE = effect, seTE =  se, studlab = study.name, data = data)
			} else{
				mt <- NA
			}
		}
	}
	
	return(mt)
}


meta.bin.complete2 = function(data, sig.level, min.study.number, sm1) {
	meta <- data %>% filter(outcome.type == "bin") %>% 
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
																		method = "linreg")$p.value, #P-value (of eggers test)
							egger.test = ifelse(pval.egger < sig.level, 1, 0),
							pval.schwarzer = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																				method = "count")$p.value, #Significance of p-value
							schwarzer.test = ifelse(pval.schwarzer < sig.level, 1, 0),
							pval.peter = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																		method = "peter")$p.value,
							peter.test = ifelse(pval.peter < sig.level, 1, 0),
							pval.harbord = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																			method = "score")$p.value,
							harbord.test = ifelse(pval.harbord < sig.level, 1, 0),
							pval.rucker = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = "ASD"), 
																		 method = "mm")$p.value,
							rucker.test = ifelse(pval.rucker < sig.level, 1, 0),
							stat.egger = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																		method.bias =  "linreg")$statistic, #Test statistic (of eggers test)
							stat.schwarzer = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																				method.bias ="count")$statistic,
							stat.peter = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																		method.bias ="peter")$statistic,
							stat.harbord = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = sm1), 
																			method.bias ="score")$statistic,
							stat.rucker = metabias(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = "ASD"), 
																		 method.bias ="mm")$statistic)
	return(meta)
}





meta.cont.complete2 = function(data, sig.level, min.study.number) {
	meta <- data %>% filter(outcome.type == "cont") %>% group_by(meta.id) %>% 
		mutate(n = n()) %>% filter(n >= 10) %>% filter(!all(mean1 == 0) & !all(mean2 == 0)) %>% 
		filter(!all(sd1 == 0) | !all(sd2 == 0)) %>% filter(meta.id < 42716 | meta.id > 42725) %>% #single patient data..
		summarize(doi = unique(doi),
							outcome.type = "cont",
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
																		method.bias ="linreg")$p.value,
							egger.test = ifelse(pval.egger < sig.level, 1, 0),
							pval.thompson = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																			 method.bias ="mm")$p.value,
							thompson.test = ifelse(pval.thompson < sig.level, 1, 0),
							pval.begg = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																	 method.bias ="rank")$p.value,
							begg.test = ifelse(pval.begg < sig.level, 1, 0),
							stat.egger = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																		method.bias ="linreg")$statistic,
							stat.thompson = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																			 method.bias ="mm")$statistic,
							stat.begg = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name), 
																	 method.bias ="rank")$statistic) 
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
																		method.bias ="linreg")$p.value,
							egger.test = ifelse(pval.egger < sig.level, 1, 0),
							pval.thompson = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																			 method.bias ="mm")$p.value,
							thompson.test = ifelse(pval.thompson < sig.level, 1, 0),
							pval.begg = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																	 method.bias ="rank")$p.value,
							begg.test = ifelse(pval.begg < sig.level, 1, 0),
							stat.egger = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																		method.bias ="linreg")$statistic,
							stat.thompson = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																			 method.bias ="mm")$statistic,
							stat.begg = metabias(metagen(TE = effect, seTE = se, studlab = study.name, sm = "HR"), 
																	 method.bias ="rank")$statistic) 
	return(meta)
}


