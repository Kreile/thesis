################################################################################################################
################################################################################################################
#Import data
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


data = pb.readData(path = PATH_DATA, file = FILE)
tmp = pb.clean(data)
data = tmp[[1]]
aliases = tmp[[2]]


################################################################################################################
################################################################################################################
# Transform Binary outcomes to odds ratios -> standardized mean difference -> correlation -> fisher scala correlation

mly.fishersz <- function(data){
  data <- data %>% mutate(odds.ratio = NA,
                          std.mean.d = NA,
                          correlation = NA,
                          fishersz = NA, 
                          fishersz.variance = NA,
                          log.OR = NA,
                          se.log.OR = NA,
                          std.m.d = NA,
                          se.std.m.d = NA, 
                          outcome.type = ifelse(outcome.measure.new == "Odds Ratio" | outcome.measure.new == "Risk Ratio", "bin", "cont"),
                          pval.type = NA,
                          pval.fishersz = NA,
                          sig.type = NA,
                          sig.fishersz = NA,
                          samplesize = ifelse(total1 > 12 & total2 > 12, 1, 0))
  
  binary.unprob <- which(data$events1 > 0 & data$events2 > 0 & data$total1 - data$events1 > 0 & data$total2 - data$events2 > 0 &
                           data$total1 + data$total2 > 3 & data$samplesize == 1)
  
  data[binary.unprob,] <- data[binary.unprob,] %>%
    mutate(log.OR = log( (events1 * (total2 - events2))/ (events2 * (total1 - events1)) ), 
           se.log.OR = sqrt( (1/events1) + 1/(total1 - events1) + (1/events2) + 1/(total2 - events2)),
           std.mean.d = log.OR * (sqrt(3)/pi),
           correlation = std.mean.d/sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
           fishersz.variance = 1/(total1 + total2 - 3),
           pval.type = 2*(1- pnorm(abs(log.OR/se.log.OR))),
           outcome.type = "bin")
  
  samplesize <- which(data$samplesize == 1)
  binary.zeros <- which(data$events1 == 0 | data$events2 == 0 & data$outcome.type == "bin")
  binary.zeros <- binary.zeros %in% samplesize
  binary.totals <- which(data$samplesize == 1 & data$total1 - data$events1 == 0 | data$total2 - data$events2 == 0 & data$outcome.type == "bin")
  
  data[binary.zeros,] <- data[binary.zeros,] %>%
    mutate(events1 = events1 + 0.5, events2 = events2 + 0.5,
      log.OR = log( (events1 * (total2 - events2))/ (events2 * (total1 - events1)) ), 
           se.log.OR = sqrt( (1/events1) + 1/(total1 - events1) + (1/events2) + 1/(total2 - events2)),
           std.mean.d = log.OR * (sqrt(3)/pi),
           correlation = std.mean.d/sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
           fishersz.variance = 1/(total1 + total2 - 3),
           pval.type = 2*(1- pnorm(abs(log.OR/se.log.OR))),
           outcome.type = "bin")
  
  data[binary.totals,] <- data[binary.totals,] %>%
    mutate(events1 = events1 - 0.5, events2 = events2 - 0.5,
           log.OR = log( (events1 * (total2 - events2))/ (events2 * (total1 - events1)) ), 
           se.log.OR = sqrt( (1/events1) + 1/(total1 - events1) + (1/events2) + 1/(total2 - events2)),
           std.mean.d = log.OR * (sqrt(3)/pi),
           correlation = std.mean.d/sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
           fishersz.variance = 1/(total1 + total2 - 3),
           pval.type = 2*(1- pnorm(abs(log.OR/se.log.OR))),
           outcome.type = "bin")
  
  
  data[!is.na(data$outcome.measure.new) & data$outcome.measure.new == "Mean Difference", ] <- data %>% 
    filter(outcome.measure.new == "Mean Difference") %>%
    mutate(se.std.mean.d = sqrt( (((total1 - 1)*sd1^2) + ((total2 - 1)*sd2^2))/(total1 + total2 - 2) ),
           std.m.d = (mean1 - mean2)/se.std.mean.d,
           correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
           fishersz.variance = 1/(total1 + total2 - 3),
           pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))),
           outcome.type = "cont"
    )
  
  data[!is.na(data$outcome.measure.new) & data$outcome.measure.new == "Std. Mean Difference", ] <- data %>% 
    filter(outcome.measure.new == "Std. Mean Difference") %>%
    mutate(se.std.m.d = sqrt( (((total1 - 1)*sd1^2) + ((total2 - 1)*sd2^2))/(total1 + total2 - 2) ),
           std.m.d = effect,
           correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
           fishersz.variance = 1/(total1 + total2 - 3),
           pval.type = 2*(1-pnorm(abs(std.m.d/se.std.m.d))),
           outcome.type = "cont")
  
  # #Filter out fishersz abs(correlation) > 1 and subgroups with fewer than 2 reproductions
  # data <- data %>% filter(fishersz < 1 & fishersz > -1) %>% filter(!is.na(fishersz))
  # data <- data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
  #   filter(fishersz < 1 & fishersz > -1) %>%
  #   mutate(counts = n()) %>% filter(counts > 1)
  data[!is.na(data$fishersz) & !is.na(data$fishersz.variance),] <- data[!is.na(data$fishersz) & !is.na(data$fishersz.variance),] %>% 
    mutate(pval.fishersz = 2*(1-pnorm(abs(fishersz/fishersz.variance))))
  data[!is.na(data$pval.type),] <- data[!is.na(data$pval.type),] %>% 
    mutate(sig.type = ifelse(pval.type < 0.05, 1, 0))
  data[!is.na(data$pval.fishersz),] <- data[!is.na(data$pval.fishersz),] %>% 
    mutate(sig.fishersz = ifelse(pval.fishersz < 0.05, 1, 0))
  return(data)
}

data.ext <- mly.fishersz(data)

mly.bin = function(data, sig.level, min.study.number) {
  meta.bin <- data %>% filter(outcome.measure.new == "Risk Ratio" | outcome.measure.new == "Odds Ratio" | outcome.measure.new == "Risk Difference") %>%  
    filter(file.nr != 3014) %>% #file.nr 3014 gibt eine seltsame Fehlermeldung, muss ausgeschlossen werden
    filter(events1 > 0 | events2 > 0) %>% #comparisons mit events  = 0 ausschliessen
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>% #comparisons mit events = total ausschliessen - vertraegt der alg. irgendwie nicht.
    filter(total1 > 12 & total2 > 12) %>% 
    group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% #Die nachfolgenden Rechnungen werden jeweils fuer diese Gruppen gemacht.
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
              lor.fixef.hkn.bin = #Pooled odds ratio estimate of fixed effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", hakn = T)$TE.fixed),
              se.lor.fixef.hkn.bin = #Pooled odds ratio estimate of fixed effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", hakn = T)$seTE.fixed),
              pval.fixef.hkn.bin = #P-value
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", hakn = T)$pval.fixed,
              sig.fixef.hkn.bin = ifelse(pval.fixef.bin > sig.level, 0, 1), #If significant
              lor.ranef.hkn.bin = #Pooled odds ratio estimate of random effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", hakn = T)$TE.random),
              se.lor.ranef.hkn.bin = #Pooled odds ratio estimate of random effects meta-analysis model
                (metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", hakn = T)$TE.random),
              pval.ranef.hkn.bin = 
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", hakn = T)$pval.random,
              sig.ranef.hkn.bin = ifelse(pval.ranef.bin > sig.level, 0, 1),
              
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
    group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% #Die nachfolgenden Rechnungen werden jeweils fuer diese Gruppen gemacht.
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
              est.fixef.hkn.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, hakn = T)$TE.fixed,
              se.fixef.hkn.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, hakn = T)$seTE.fixed,
              pval.fixef.hkn.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, hakn = T)$pval.fixed,
              sig.fixef.hkn.cont = ifelse(pval.fixef.cont > sig.level, 0, 1),
              est.ranef.hkn.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, hakn = T)$TE.random,
              se.ranef.hkn.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, hakn = T)$seTE.random,
              pval.ranef.hkn.cont = 
                metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name, hakn = T)$pval.random,
              sig.ranef.hkn.cont = ifelse(pval.ranef.cont > sig.level, 0, 1),
              
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
                                      est.fixef.hkn = ifelse(!is.na(lor.fixef.hkn.bin), lor.fixef.hkn.bin, est.fixef.hkn.cont),
                                      est.ranef.hkn = ifelse(!is.na(lor.ranef.hkn.bin), lor.ranef.hkn.bin, est.ranef.hkn.cont),
                                      pval.fixef.hkn = ifelse(!is.na(pval.fixef.hkn.bin), pval.fixef.hkn.bin, pval.fixef.hkn.cont),
                                      pval.ranef.hkn = ifelse(!is.na(pval.ranef.hkn.bin), pval.ranef.hkn.bin, pval.ranef.hkn.cont),
                                      sig.fixef.hkn = ifelse(!is.na(sig.fixef.hkn.bin), sig.fixef.hkn.bin, sig.fixef.hkn.cont),
                                      sig.ranef.hkn = ifelse(!is.na(sig.ranef.hkn.bin), sig.ranef.hkn.bin, sig.fixef.hkn.cont),
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
file = "mly.RData"
save(mly, file =  file.path(PATH_RESULTS, file))

print(data.ext %>% group_by(file.nr, comparison.nr, outcome.name, subgroup.nr) %>% distinct(outcome.type), n = 400)
print(data.ext %>% filter(file.nr == 2 & outcome.name == "Resuturing of wound - up to 3 months" & subgroup.nr == 1) %>% select(outcome.measure, events1, events2, total1, total2))
