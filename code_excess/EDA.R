rm(list = ls())
PATH_HOME = path.expand("~") # user home
PATH = file.path(PATH_HOME, 'Data/PubBias')
PATH2 = file.path(PATH_HOME, 'PubBias')
FILE = 'cochrane_2019-07-04.csv'
PATH_DATA = file.path(PATH, 'data')
PATH_CODE = file.path(PATH2, 'code')
PATH_RESULTS = file.path(PATH2, 'results_new')
PATH_FIGURES = file.path(PATH_RESULTS, 'figures')

require(lme4)

source(file.path(PATH_CODE, 'PubBias_functions.R'))

load(file.path(PATH_RESULTS, "data2.RData"))
# data.ext2 <- pb.process3(data)
load(file.path(PATH_RESULTS, "data_used_for_analysis.RData"))

load(file.path(PATH_RESULTS, "meta_id_vector.RData"))


eda.dat <- data.ext2 %>% filter(meta.id %in% meta.id.vector) #Filter meta-analysis data
#--------------------------------------------------------------------------------------------------------------------#

#Results meta-analysis properties:
studies.per.meta <- c(quantile(x = meta.f$k, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), mean(meta.f$k))
totalsample.per.meta <- c(quantile(x = meta.f$total.samplesize, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), mean(meta.f$total.samplesize))
I2.per.meta <- c(quantile(x = meta.f$I2, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), mean(meta.f$I2))
tau.per.meta <- c(quantile(x = meta.f$tau, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), mean(meta.f$tau))
variance.ratio.per.meta <- c(quantile(x = meta.f$var.ratio, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), mean(meta.f$var.ratio))

total.events.per.meta <- c(quantile(x = meta.f$tau, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = T), mean(meta.f$tau, na.rm = T))

meta.properties <- rbind(studies.per.meta, totalsample.per.meta, total.events.per.meta, 
                         I2.per.meta, tau.per.meta, variance.ratio.per.meta)
colanmes(meta.properties) 
rownames(meta.properties) <- c("Study number", "Total sample size", "Total event counts", 
                               "I-squared", "Tau variance", "Ratio max/min variance")
#--------------------------------------------------------------------------------------------------------------------#

#Other properties:
z.score.excluded.results <- data.ext2 %>% group_by(meta.id) %>% mutate(n = n()) %>% filter(n > 9) %>%
  filter(outcome.flag != "IV") %>% filter(total1 + total2 < 4) %>%
  filter(!is.na(effect)) %>% filter(!is.na(se)) %>% ungroup() %>% count()

#--------------------------------------------------------------------------------------------------------------------#


#Exploratory plots:
hist(meta.f$mean.publication.year, nclass = 50)
meta.f %>% filter(mean.publication.year <= 2010 & mean.publication.year >= 1980) %>% 
  mutate(pb.test = factor(ifelse(pval.se < 0.1, "significant", "non-significant"))) %>% 
  ggplot(aes(stat(count), x = mean.publication.year, fill = pb.test)) + 
  geom_density(position = "fill")

hist(meta.f$k, nclass = 30)
meta.f %>% filter(k < 30) %>% 
  mutate(pb.test = factor(ifelse(pval.se < 0.1, "significant", "non-significant"))) %>% 
  ggplot(aes(stat(count), x = k, fill = pb.test)) + 
  geom_density(position = "fill")

hist(meta.f$est.d.fixef, nclass = 30)
meta.f %>% filter(abs(est.d.fixef) < 3) %>% 
  mutate(pb.test = factor(ifelse(pval.se < 0.1, "significant", "non-significant"))) %>% 
  ggplot(aes(stat(count), x = est.d.fixef, fill = pb.test)) + 
  geom_density(position = "fill")

hist(meta.f$n.sig.single2, nclass = 30)
meta.f %>% filter(n.sig.single < 10) %>% 
  mutate(pb.test = factor(ifelse(pval.se < 0.1, "significant", "non-significant"))) %>% 
  ggplot(aes(stat(count), x = n.sig.single2, fill = pb.test)) + 
  geom_density(position = "fill")














eda.dat %>% filter(outcome.flag == "DICH") %>% group_by(meta.id) %>%
  filter(all(events1 == 0) & all(events2 == 0)) %>% ungroup() %>% count()

eda.dat %>% filter(outcome.flag == "CONT") %>%
  filter(all(mean1 == 0) & all(mean2 == 0)) %>% ungroup() %>% count()

eda.dat %>% filter(outcome.flag == "CONT") %>%
  filter(all(sd1 == 0) | all(sd2 == 0)) %>% ungroup() %>% count()

eda.dat %>% filter(outcome.flag == "IV") %>%
  filter(all(effect == 0) | all(se == 0)) %>% ungroup() %>% count()

















