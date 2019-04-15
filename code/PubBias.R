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

load(file.path(PATH_RESULTS, 'PubBias_results1.RData'))
load(file.path(PATH_RESULTS, 'PubBias_results2.RData'))


data = pb.readData(path = PATH_DATA, file = FILE)
tmp = pb.clean(data)
data = tmp[[1]]
aliases = tmp[[2]]

require(biostatUZH)
require(tidyverse)
require(meta)
require(xtable)



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
              trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$k0 / n(),
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
              trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$k0 / n(),
              Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name)$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() %>% select(-outcome.nr, -subgroup.nr)
  return(metadat)
}


meta.bin.complete <- function(data, min.study.number, sig.level){
  metadat <- data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% filter(file.nr != 3014) %>% 
    filter(events1 > 0 | events2 > 0) %>% 
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), n = n(),
              OR.fixef.bin = 
                exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.fixed),
              pval.fixef.bin = 
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.fixed,
              OR.ranef.bin = 
                exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.random),
              pval.fixef.bin = 
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.random,
              OR.trimfill.bin = 
                exp(trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE),
              pval.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                        method = "linreg")$p.val,
              egger.test = ifelse(pval.egger.bin < sig.level, 1, 0),
              pval.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                            method = "count")$p.val,
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
                                        method = "linreg")$statistic,
              stat.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                            method = "count")$statistic,
              stat.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                        method = "peter")$statistic,
              stat.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "score")$statistic,
              stat.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), 
                                         method = "mm")$statistic,
              trim.bin = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name))$k0 / n(),
              Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() %>% select(-outcome.nr, -subgroup.nr)
  return(metadat)
} 


#Proportion of pbias according to 5% significance level:
rejection.bin <- data.bin.ext %>% summarize(egger.rejection = mean(egger.test),
                                            schwarzer.rejection = mean(schwarzer.test),
                                            rucker.rejection = mean(rucker.test),
                                            harbord.rejection = mean(harbord.test),
                                            peter.rejection = mean(peter.test))
rejection.bin <- rejection.bin %>% gather(key = test.type)
rejection.bin %>% ggplot(aes(y = value, x = test.type)) + geom_col()
bin.spread <- data.bin.ext %>% select(egger.test, schwarzer.test, harbord.test, peter.test) %>% gather(key = test.type, value = "null.hypothesis")
bin.spread$value <- ifelse(bin.spread$value == 1, "reject", "accept")
bin.spread$value <- factor(bin.spread$value)
bin.spread %>% ggplot(aes(x = test.type, fill = value)) + geom_bar() + coord_flip()

rejection.cont <- data.cont.ext %>% summarize(egger.rejection = mean(egger.test),
                                              begg.rejection = mean(begg.test),
                                              thomson.rejection = mean(thomson.test))

cont.spread <- data.cont.ext %>% select(egger.test, begg.test, thomson.test) %>% gather(key = test.type, value = "null.hypothesis")
cont.spread$value <- ifelse(cont.spread$value == 1, "rejected", "accepted")
cont.spread$value <- factor(cont.spread$value)
cont.spread %>% ggplot(aes(x = test.type, fill = value)) + geom_bar() + coord_flip() + theme_bw()



#Agreement proportions of publication bias tests:

#Binary:
data.bin.ext <- data.bin.ext %>% mutate(egger.schwarzer = ifelse(egger.test == schwarzer.test, "agree", "disagree"),
                                        egger.peter = ifelse(egger.test == peter.test, "agree", "disagree"),
                                        egger.rucker = ifelse(egger.test == rucker.test, "agree", "disagree"),
                                        egger.harbord = ifelse(egger.test == harbord.test, "agree", "disagree"),
                                        schwarzer.peter = ifelse(schwarzer.test == peter.test, "agree", "disagree"),
                                        schwarzer.rucker = ifelse(schwarzer.test == rucker.test, "agree", "disagree"),
                                        schwarzer.harbord = ifelse(schwarzer.test == harbord.test, "agree", "disagree"),
                                        rucker.peter = ifelse(rucker.test == peter.test, "agree", "disagree"),
                                        rucker.harbord = ifelse(rucker.test == harbord.test, "agree", "disagree"),
                                        harbord.peter = ifelse(harbord.test == peter.test, "agree", "disagree"))

agreement.bin <- data.bin.ext %>% summarise(egger.schwarzer = sum(egger.schwarzer == "agree")/n(),
                                            egger.peter = sum(egger.peter == "agree")/n(),
                                            egger.rucker = sum(egger.rucker == "agree")/n(),
                                            egger.harbord = sum(egger.harbord == "agree")/n(),
                                            schwarzer.peter = sum(schwarzer.peter == "agree")/n(),
                                            schwarzer.rucker = sum(schwarzer.rucker == "agree")/n(),
                                            schwarzer.harbord = sum(schwarzer.harbord == "agree")/n(),
                                            rucker.peter = sum(rucker.peter == "agree")/n(),
                                            harbord.peter = sum(harbord.peter == "agree")/n())

correlation.bin <- data.bin.ext %>% summarise(egger.schwarzer = cor(pval.egger.bin, pval.schwarzer.bin),
                                              egger.peter = cor(pval.egger.bin, pval.peter.bin),
                                              egger.rucker = cor(pval.egger.bin, pval.rucker.bin),
                                              egger.harbord = cor(pval.egger.bin, pval.harbord.bin),
                                              schwarzer.peter = cor(pval.schwarzer.bin, pval.peter.bin),
                                              schwarzer.rucker = cor(pval.schwarzer.bin, pval.rucker.bin),
                                              schwarzer.harbord = cor(pval.schwarzer.bin, pval.harbord.bin),
                                              rucker.peter = cor(pval.rucker.bin, pval.peter.bin),
                                              harbord.peter = cor(pval.harbord.bin, pval.peter.bin))

binary.tests.agreement <- rbind(agreement.bin, correlation.bin)
rownames(binary.tests.agreement) <- c("Test Agreement","P-value Correlation")


#Continuous:

data.cont.ext <- data.cont.ext %>% mutate(thomson.egger = ifelse(thomson.test == egger.test, "agree", "disagree"),
                                          thomson.begg = ifelse(thomson.test == begg.test, "agree", "disagree"),
                                          egger.begg = ifelse(egger.test == begg.test, "agree", "disagree"))

agreement.cont <- data.cont.ext %>% summarise(thomson.egger = sum(thomson.egger == "agree")/n(),
                                              thomson.begg = sum(thomson.begg == "agree")/n(),
                                              egger.begg = sum(egger.begg == "agree")/n())

correlation.cont <- data.cont.ext %>% summarise(thomson.egger = cor(pval.thomson.cont, pval.egger.cont),
                                                thomson.begg = cor(pval.thomson.cont, pval.begg.cont),
                                                egger.begg = cor(pval.egger.cont, pval.begg.cont))

cont.tests.agreement <- rbind(agreement.cont, correlation.cont)
rownames(cont.tests.agreement) <- c("Test Agreement","P-value Correlation")




#Plots:
#Trimfill proportion and test results:
trimfill.cont.mean <- data.cont.ext %>% ungroup %>% summarise(mean = mean(missing.trim.cont)) %>% select(mean)
trimfill.cont.median <- data.cont.ext %>% ungroup %>% summarise(median = median(missing.trim.cont)) %>% select(median)

trimfill.bin.mean <- data.bin.ext %>% ungroup %>% summarise(mean = mean(missing.trim.bin)) %>% select(mean)
trimfill.bin.median <- data.bin.ext %>% ungroup %>% summarise(median = median(missing.trim.bin)) %>% select(median)

data.bin.ext %>% ggplot(aes(x = missing.trim.bin)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Fraction of Trimmed Comparisons (binary outcome)") + xlab("Fraction") + ylab("Frequency")

data.cont.ext %>% ggplot(aes(x = missing.trim.cont)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Fraction of Trimmed Comparisons (continuous outcome)") + xlab("Fraction") + ylab("Frequency")



cdplot(x = data.bin.ext$trim.bin, y = as.factor(ifelse(data.bin.ext$peter.test == 1, "yes", "no")), xlab = "Fraction of added studies", ylab = "Bias")
cdplot(x = data.bin.ext$trim.bin, y = as.factor(data.bin.ext$harbord.test))
cdplot(x = data.bin.ext$trim.bin, y = as.factor(data.bin.ext$schwarzer.test))

cdplot(x = data.cont.ext$trim.cont, y = as.factor(ifelse(data.cont.ext$thomson.test == 1, "yes", "no")), xlab = "Fraction of added studies", ylab = "Bias")


####################################################################################################
# Correction for small study effects
####################################################################################################




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
              # trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, studlab = study.name))$k0 / n(),
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



pb.corr.bin(data[c(1:1000),], 10)

ds <- data %>% filter(file.nr == 11) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  filter(outcome.measure == "Odds Ratio" | outcome.measure == "Risk Ratio") %>% filter(events1 > 0 & events2 > 0) %>% 
  mutate(counts = n()) %>% filter(counts == 18)
print(ds %>% select(effect, events1, events2, study.name, outcome.name, counts), n = 200)

meta.ds <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, data = ds, studlab = study.name, sm = "OR")
meta.ds



# copas.ds <- copas(meta.ds)
# 
# pb.corr.bin(ds, 10)
# 
# copas.ds$pval.rsb
# which.min(abs(0.05 - copas.ds$pval.rsb
#               [copas.ds$pval.rsb > 0.05]))


meta.bin.complete <- function(data, min.study.number, sig.level){
  metadat <- data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% filter(file.nr != 3014) %>% 
    filter(events1 > 0 | events2 > 0) %>% 
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
    mutate(study.year = ifelse(study.year > 1950, study.year, NA)) %>% 
    summarize(doi = unique(doi), n = n(),
              mean.samplesize = mean(total1 + total2, na.rm = T),
              total.samplesize = sum(total1 + total2),
              mean.publication.year = mean(study.year, na.rm = TRUE),
              first.publication.year = min(study.year, na.rm = T),
              OR.fixef.bin = 
                exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.fixed),
              pval.fixef.bin = 
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.fixed,
              sig.fixef.bin = ifelse(pval.fixef.bin > sig.level, 0, 1),
              OR.ranef.bin = 
                exp(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$TE.random),
              pval.ranef.bin = 
                metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR")$pval.random,
              sig.ranef.bin = ifelse(pval.ranef.bin > sig.level, 0, 1),
              OR.trimfill.fixef.bin = 
                trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE.fixed,
              OR.trimfill.random.bin = 
                trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE.random,
              OR.reg.random.bin = 
                exp(limitmeta(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR"))$TE.adjust),
              pval.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                        method = "linreg")$p.val,
              egger.test = ifelse(pval.egger.bin < sig.level, 1, 0),
              pval.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                            method = "count")$p.val,
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
                                        method = "linreg")$statistic,
              stat.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                            method = "count")$statistic,
              stat.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                        method = "peter")$statistic,
              stat.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR"), 
                                          method = "score")$statistic,
              stat.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "ASD"), 
                                         method = "mm")$statistic,
              trim.bin = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name))$k0 / n(),
              Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm = "OR")$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() %>% select(-outcome.nr, -subgroup.nr)
  return(metadat)
} 


meta.bin.complete <- meta.bin.complete(data, 10, 0.05)



meta.bin.complete %>% 
  filter(sig.fixef.bin == 1) %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis )) + geom_bar() + coord_flip() + theme_bw()

meta.bin.complete %>% 
  filter(sig.fixef.bin == 0) %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis )) + geom_bar() + coord_flip() + theme_bw()


meta.bin.complete %>% 
  filter(sig.ranef.bin == 1) %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis )) + geom_bar() + coord_flip() + theme_bw()

meta.bin.complete %>% 
  filter(sig.ranef.bin == 0) %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis )) + geom_bar() + coord_flip() + theme_bw()


meta.bin.complete %>% filter(mean.publication.year < 2019 & mean.publication.year > 1960) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef.bin))) + geom_histogram()

meta.bin.complete %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef.bin))) + geom_density(na.rm = T, position = "fill")

meta.bin.complete %>% filter(first.publication.year < 2013 & first.publication.year > 1990) %>% 
  ggplot(aes(x = first.publication.year, fill = factor(sig.fixef.bin))) + geom_density(na.rm = T, position = "fill")

meta.bin.complete %>% filter(first.publication.year < 2013 & first.publication.year > 1990) %>% 
  ggplot(aes(x = first.publication.year, fill = factor(peter.test))) + geom_density(na.rm = T, position = "fill")

meta.bin.complete %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(rucker.test))) + geom_density(na.rm = T, position = "fill")

data %>% filter(study.year < 2019 & study.year > 1980) %>% group_by(study.year) %>% summarize(samplesize = mean(total1 + total2, na.rm = T)) %>% 
  ggplot(aes(y = samplesize, x = study.year)) + geom_line()






egger.schwarzer <- data.bin.ext %>% group_by(egger.schwarzer) %>% count() %>% ggplot(aes(x = egger.schwarzer, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.peter <- data.bin.ext %>% group_by(egger.peter) %>% count() %>% ggplot(aes(x = egger.peter, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.rucker <- data.bin.ext %>% group_by(egger.rucker) %>% count() %>% ggplot(aes(x = egger.rucker, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.harbord <- data.bin.ext %>% group_by(egger.harbord) %>% count() %>% ggplot(aes(x = egger.harbord, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.rucker <- data.bin.ext %>% group_by(egger.rucker) %>% count() %>% ggplot(aes(x = egger.rucker, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.rucker <- data.bin.ext %>% group_by(egger.rucker) %>% count() %>% ggplot(aes(x = egger.rucker, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.rucker <- data.bin.ext %>% group_by(egger.rucker) %>% count() %>% ggplot(aes(x = egger.rucker, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")

library(grid)
library(gridExtra)
grid.arrange(egger.schwarzer, egger.peter, egger.rucker)

data.cont %>% ggplot(aes(y = pval.egger.cont, x = pval.thomson.cont)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
data.cont %>% ggplot(aes(x = pval.thomson.cont, y = pval.begg.cont)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
data.cont.ext %>% filter(n < 100) %>% ggplot(aes(y = pval.thomson.cont, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess", se = F) 
data.cont %>% filter(n < 100) %>% ggplot(aes(y = pval.egger.cont, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess", se = F) 
data.cont.ext %>% ggplot(aes(y = stat.egger.cont, x = stat.thomson.cont)) + geom_point(alpha = 0.3) + theme_bw() #+ geom_smooth(method = "lm")  
data.cont %>% ggplot(aes(x = pval.thomson.cont, y = trim.cont)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 

data.bin %>% ggplot(aes(x = pval.harbord.bin, y = pval.peter.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
data.bin.ext %>% filter(n < 150)%>% ggplot(aes(y = pval.peter.bin, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
data.bin %>% filter(n < 150)%>% ggplot(aes(y = pval.harbord.bin, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(y = pval.egger.bin, x = pval.peter.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(y = pval.egger.bin, x = pval.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(x = pval.schwarzer.bin, y = pval.peter.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(x = pval.schwarzer.bin, y = pval.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(y = pval.rucker.bin, y = pval.egger.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(y = pval.rucker.bin, y = pval.peter.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(y = pval.rucker.bin, y = pval.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 

data.bin.ext %>% ggplot(aes(x = stat.peter.bin, y = stat.schwarzer.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(x = stat.peter.bin, y = stat.egger.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(x = stat.peter.bin, y = stat.rucker.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(x = stat.peter.bin, y = stat.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess")

data.bin.ext %>% ggplot(aes(y = stat.egger.bin, x = stat.schwarzer.bin)) + geom_point(alpha = 0.5) + theme_bw()# + geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(y = stat.rucker.bin, x = stat.schwarzer.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(y = stat.harbord.bin, x = stat.schwarzer.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 

data.bin.ext %>% ggplot(aes(x = stat.rucker.bin, y = stat.harbord.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
data.bin.ext %>% ggplot(aes(x = stat.rucker.bin, y = stat.egger.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 

data.bin.ext %>% ggplot(aes(x = stat.harbord.bin, y = stat.egger.bin)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess")



# data.bin.plot <- data.bin %>% ungroup() %>% select(pval.peter.bin, pval.harbord.bin, trim.bin, Q, I2)
# data.cont.plot <- data.cont %>% ungroup() %>% select(pval.egger.cont, pval.thomson.cont, pval.begg.cont, trim.cont, Q, I2)
# scatterplotMatrix(data.bin.plot)
# scatterplotMatrix(data.cont.plot)
# pairs(data.bin.plot)

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

