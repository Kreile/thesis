rm(list = ls())
PATH_HOME = path.expand("~") # user home
PATH = file.path(PATH_HOME, 'Data/PubBias')
PATH2 = file.path(PATH_HOME, 'PubBias')
FILE = 'cochrane_2018-06-09.csv'
PATH_DATA = file.path(PATH, 'data')
PATH_CODE = file.path(PATH2, 'code')
PATH_RESULTS = file.path(PATH, 'results')
PATH_FIGURES = file.path(PATH_RESULTS, 'figures')

file_results = "pb.RData"

source(file.path(PATH_CODE, 'PubBias_functions.R'))

data = pb.readData(path = PATH_DATA, file = FILE)
tmp = pb.clean(data)
data = tmp[[1]]
aliases = tmp[[2]]

require(biostatUZH)
require(tidyverse)
require(meta)
require(xtable)
library(car)

####################################################################################################
#Look for sample size effect on effect size
####################################################################################################

data %>% filter(total1 + total2 < 100 & total1 + total2 > 5) %>%
  mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>%
  group_by(sample.size) %>% 
  summarize(mean.scaled.effect = mean(scaled.effect, na.rm = T)) %>% 
  ggplot(aes(x = sample.size, y = mean.scaled.effect)) + geom_point() + geom_smooth(method = "lm")

data %>%  filter(total1 + total2 < 100 & total1 + total2 > 5) %>%
  mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>%
  ggplot(aes(x = sample.size, y = scaled.effect)) + geom_point(alpha = 0.15, size = 0.5) + geom_smooth(method = "lm")

data %>% filter(total1 + total2 > 5 & total1 + total2 < 300) %>% 
  mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>%
  group_by(sample.size) %>% 
  summarize(mean.effect = mean(scaled.effect, na.rm = T)) %>% 
  ggplot(aes(x = sample.size, y = mean.effect)) + geom_point() + geom_smooth(method = "lm") + geom_smooth(method = "loess", color = "red")

# Effect sizes scaled meta-analyses reviews.
# data %>% group_by(file.nr, outcome.nr, subgroup.nr) %>% 
#   filter(n() > 2) %>% 
#   mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>% 
#   filter(sample.size < 500 & sample.size > 5) %>% 
#   ggplot(aes(x = sample.size, y = scaled.effect)) + geom_point(alpha = 0.3, size = 0.1) + geom_smooth(method = "lm") #+ geom_smooth(method = "loess")
# 
# 
# data %>% group_by(file.nr, outcome.nr, subgroup.nr) %>% 
#   filter(n() > 2) %>% 
#   mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>%
#   group_by(sample.size) %>% 
#   summarize(mean.effect = mean(scaled.effect, na.rm = T)) %>% 
#   filter(sample.size < 500) %>% 
#   ggplot(aes(x = sample.size, y = mean.effect)) + geom_point() + geom_smooth(method = "lm")




####################################################################################################
#Meta-analysis dataset creation function (p-values, heterogeneity statistics, trim-fill number of trimmed studies).
####################################################################################################

pb.bias.cont <- function(data, min.study.number){
  metadat <- data %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>%# filter(file.nr < 503) %>% 
    filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
    filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), n = n(),
              pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                         method = "linreg")$p.val,
              pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                           method = "mm")$p.val,
              pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                        method = "rank")$p.val,
              trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2))$k0 / n(),
              Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2)$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() %>% select(-outcome.nr, -subgroup.nr)
  return(metadat)
}

pb.bias.bin <- function(data, min.study.number){
  metadat <- data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% filter(file.nr != 3014) %>% 
    filter(events1 > 0 | events2 > 0) %>% 
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
    group_by(file.nr, outcome.nr, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), n = n(),
              pval.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "peter")$p.val,
              pval.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "score")$p.val,
              trim.bin = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2))$k0 / n(),
              Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR")$Q,
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
              pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                         method = "linreg")$p.val,
              egger.test = ifelse(pval.egger.cont < sig.level, 1, 0),
              pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                           method = "mm")$p.val,
              thomson.test = ifelse(pval.thomson.cont < sig.level, 1, 0),
              pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                        method = "rank")$p.val,
              begg.test = ifelse(pval.begg.cont < sig.level, 1, 0),
              stat.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                         method = "linreg")$statistic,
              stat.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                           method = "mm")$statistic,
              stat.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                        method = "rank")$statistic,
              trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2))$k0 / n(),
              Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2)$Q,
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
              pval.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "linreg")$p.val,
              egger.test = ifelse(pval.egger.bin < sig.level, 1, 0),
              pval.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "count")$p.val,
              schwarzer.test = ifelse(pval.schwarzer.bin < sig.level, 1, 0),
              pval.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "peter")$p.val,
              peter.test = ifelse(pval.peter.bin < sig.level, 1, 0),
              pval.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "score")$p.val,
              harbord.test = ifelse(pval.harbord.bin < sig.level, 1, 0),
              pval.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "ASD"), method = "mm")$p.val,
              rucker.test = ifelse(pval.rucker.bin < sig.level, 1, 0),
              stat.egger.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "linreg")$statistic,
              stat.schwarzer.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "count")$statistic,
              stat.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "peter")$statistic,
              stat.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "score")$statistic,
              stat.rucker.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "ASD"), method = "mm")$statistic,
              trim.bin = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2))$k0 / n(),
              Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR")$Q,
              I2  = max(0, (Q - n + 1)/Q)) %>% 
    ungroup() %>% select(-outcome.nr, -subgroup.nr)
  return(metadat)
} 

# data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% filter(file.nr != 3014) %>% 
#   filter(events1 > 0 | events2 > 0) %>% 
#   filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
#   group_by(file.nr, outcome.nr, subgroup.nr) %>% 
#   mutate(n = n()) %>% filter(n > 9) %>% 
#   summarize(doi = unique(doi),
#     pval.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "peter")$p.val,
#     pval.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "score")$p.val,
#     trim.bin = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2))$k0 / n(),
#     Q.bin = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR")$Q,
#     I2  = max(0, (Q - n + 1)/Q))

# data %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>%# filter(file.nr < 503) %>% 
#   filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
#   filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
#   group_by(file.nr, outcome.nr, subgroup.nr) %>% 
#   mutate(n = n()) %>% filter(n > 9) %>% 
#   summarize(doi = unique(doi),
#     pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
#                             method = "linreg")$p.val,
#             pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
#                             method = "mm")$p.val,
#             pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
#                             method = "rank")$p.val,
#             trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2))$k0 / n(),
#             Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2)$Q) 


####################################################################################################
#Investigation of similarity of tests and test results
####################################################################################################




#Make dataset to store all those measures for easier comparison

data.bin <- pb.bias.bin(data, 10)
data.cont <- pb.bias.cont(data, 10)  

data.bin.ext <- pb.bias.bin.extended(data, min.study.number = 10, sig.level = 0.05)
max(data.bin.ext$stat.peter.bin)
data.bin.ext <- data.bin.ext %>% filter(stat.peter.bin < 10000)
data.cont.ext <- pb.bias.cont.extended(data, min.study.number = 10, sig.level = 0.05)

save(data.bin.ext, file = "PubBias_results1.RData")
save(data.cont.ext, file = "PubBias_results2.RData")


#Proportion of pbias according to 5% significance level:
rejection.bin <- data.bin.ext %>% summarize(egger.rejection = mean(egger.test),
                                            schwarzer.rejection = mean(schwarzer.test),
                                            rucker.rejection = mean(rucker.test),
                                            harbord.rejection = mean(harbord.test),
                                            peter.rejection = mean(peter.test))
rejection.bin <- rejection.bin %>% gather(key = test.type)
rejection.bin %>% ggplot(aes(y = value, x = test.type)) + geom_col()
bin.spread <- data.bin.ext %>% select(egger.test, schwarzer.test, harbord.test, peter.test) %>% gather(key = test.type)
bin.spread$value <- ifelse(bin.spread$value == 1, "reject", "accept")
bin.spread$value <- factor(bin.spread$value)
bin.spread %>% ggplot(aes(x = test.type, fill = value)) + geom_bar() + coord_flip()

rejection.cont <- data.cont.ext %>% summarize(egger.rejection = mean(egger.test),
                                              begg.rejection = mean(begg.test),
                                              thomson.rejection = mean(thomson.test))

cont.spread <- data.cont.ext %>% select(egger.test, begg.test, thomson.test) %>% gather(key = test.type)
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
trimfill.cont <- data %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% #filter(file.nr < 50) %>% 
  filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
  filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
  group_by(file.nr, outcome.nr, subgroup.nr) %>% 
  mutate(n = n()) %>% filter(n > 9) %>% 
  summarize(trim = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2))$k0 / n()) 

trimfill.cont.mean <- trimfill.cont %>% ungroup %>% summarise(mean = mean(trim)) %>% select(mean)
trimfill.cont.median <- trimfill.cont %>% ungroup %>% summarise(median = median(trim)) %>% select(median)

trimfill.bin <- data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% #filter(file.nr < 503) %>% 
  filter(events1 > 0 | events2 > 0) %>% 
  filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
  group_by(file.nr, outcome.nr, subgroup.nr) %>% 
  mutate(n = n()) %>% filter(n > 9) %>% 
  summarize(trim = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2))$k0 / n())

trimfill.bin.mean <- trimfill.bin %>% ungroup %>% summarise(mean = mean(trim)) %>% select(mean)
trimfill.bin.median <- trimfill.bin %>% ungroup %>% summarise(median = median(trim)) %>% select(median)

trimfill.bin %>% ggplot(aes(x = trim)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Fraction of Trimmed Comparisons (binary outcome)") + xlab("Fraction") + ylab("Frequency")

trimfill.cont %>% ggplot(aes(x = trim)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Fraction of Trimmed Comparisons (continuous outcome)") + xlab("Fraction") + ylab("Frequency")



cdplot(x = data.bin.ext$trim.bin, y = as.factor(ifelse(data.bin.ext$peter.test == 1, "yes", "no")), xlab = "Fraction of added studies", ylab = "Bias")
cdplot(x = data.bin.ext$trim.bin, y = as.factor(data.bin.ext$harbord.test))
cdplot(x = data.bin.ext$trim.bin, y = as.factor(data.bin.ext$schwarzer.test))

cdplot(x = data.cont.ext$trim.cont, y = as.factor(ifelse(data.cont.ext$thomson.test == 1, "yes", "no")), xlab = "Fraction of added studies", ylab = "Bias")


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
#   summarize(pval.thomson.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
#                                          method = "mm")$p.val,
#             Q = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2)$Q, n = n(),
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
#   summarize(pval.peter.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "peter")$p.val,
#             Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR")$Q, n = n(),
#             I2  = max(0, (Q - n + 1)/Q)) %>%
#   ggplot(aes(x = I2, y = pval.peter.bin)) + geom_point(alpha = 0.7) + theme_bw() + geom_smooth(method = "lm")

