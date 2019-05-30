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


# data = pb.readData(path = PATH_DATA, file = FILE)
# tmp = pb.clean(data)
# data = tmp[[1]]
# aliases = tmp[[2]]

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


file.bin <- "pb.bin.RData"
if (file.exists(file.path(PATH_RESULTS, file.bin))) {
  load(file.path(PATH_RESULTS, file.bin))
} else {
  meta.bin <- meta.bin.complete(data.ext, min.study.number = 10, sig.level = 0.05, sm = "OR")
  save(meta.bin, file =  file.path(PATH_RESULTS, file.bin))
}

file.cont <- "pb.cont.RData"
if (file.exists(file.path(PATH_RESULTS, file.cont))) {
  load(file.path(PATH_RESULTS, file.cont))
} else {
  meta.cont <- meta.cont.complete(data.ext, min.study.number = 10, sig.level = 0.05)
  save(meta.cont, file =  file.path(PATH_RESULTS, file.cont))
}

file.meta <- "meta.RData"
if (file.exists(file.path(PATH_RESULTS, file.meta))) {
  load(file.path(PATH_RESULTS, file.meta))
} else {
  meta <- pb.meta.merge(meta.bin, meta.cont)
  save(meta, file =  file.path(PATH_RESULTS, file.meta))
}


file.cont <- "mly.cont.RData"
if (file.exists(file.path(PATH_RESULTS, file.cont))) {
  load(file.path(PATH_RESULTS, file.cont))
} else {
  data.cont <- mly.cont(data.ext, 0.05, min.study.number = 2)
  save(data.cont, file =  file.path(PATH_RESULTS, file.cont))
}

file.bin <- "mly.bin.RData"
if (file.exists(file.path(PATH_RESULTS, file.bin))) {
  load(file.path(PATH_RESULTS, file.bin))
} else {
  data.bin <- mly.bin(data.ext, 0.05, min.study.number = 2)
  save(data.bin, file =  file.path(PATH_RESULTS, file.bin))
}


load(file.path(PATH_RESULTS, file = "mly.RData"))
load(file.path(PATH_RESULTS, file = "data.processed.RData"))

require(biostatUZH)
require(tidyverse)
require(meta)
require(metasens)
require(gridExtra)


data %>% #distinct(study.name, .keep_all = T) %>%  
  filter(study.year < 2018 & study.year > 1980) %>% group_by(study.year) %>% 
  summarize(samplesize = median(total1 + total2, na.rm = T)) %>% 
  ggplot(aes(y = samplesize, x = study.year)) + geom_line()

########################################################################################################################
########################################################################################################################
#Cumulative Meta-analysis:
########################################################################################################################
########################################################################################################################

#Prepare data: Exclude n < 24, unknown study years, and study sample size < 3
data.cum <- data.ext %>% group_by(meta.id) %>% filter(n() > 2) %>% 
  filter(!is.na(study.year) & study.year < 2020 & study.year > 1945) %>% 
  filter(total1 > 11 & total2 > 11)

tmp <- data.cum %>% filter(meta.id == 312 | meta.id == 310) %>% select(study.name, study.year, events1, total1, events2, total2, study.id) 


#Create cumulative meta-analysis function: "Only do a meta-analysis for studies with study year <= study year of given study"
#(Issue with ties: will give same results for all studies with equal study years.)

#Binary outcomes
cum.bin.zval.fixed <- function(study.year, study.name, events1, events2, total1, total2, sm){
  results <- c()
  for(u in seq_along(study.name)){ #Loop through the studies
    year <- study.year[u] #Pick study year of given study
    
    if(length(which(study.year <= year)) < 2){ #If there are no earlier published studies
      results <- c(results, NA)                         #give NA
    } else{ #Else analyse the study together with studies published earlier by meta-analysis
      index.lower.year <- which(study.year <= year)
      TE <- metabin(event.e = events1[index.lower.year], n.e = total1[index.lower.year], 
                    event.c = events2[index.lower.year], n.c = total2[index.lower.year], 
                    studlab = study.name[index.lower.year], sm= sm)$pval.fixed #log of treatment effect of fixed effects m.a.
      results <- c(results, TE)
    }
  } 
  return(results)
}

# tmp <- data.cum %>% filter(meta.id == 86 | meta.id == 87) %>% select(study.name, study.year, events1, total1, events2, total2, effect) 
# tmp %>% mutate(zval = exp(cum.zval.fixed(study.year, study.name, events1, events2, total1, total2, sm = "RR")))

#Continuous outcomes
cum.cont.zval.fixed <- function(study.year, study.name, mean1, sd1, total1, mean2, sd2, total2){
  results <- c()
  for(u in seq_along(study.name)){ #Loop through the studies
    year <- study.year[u] #Pick study year of given study
    
    if(length(which(study.year <= year)) < 2){ #If there are no earlier published studies
      results <- c(results, NA)                         #give NA
    } else{ #Else analyse the study together with studies published earlier by meta-analysis
      index.lower.year <- which(study.year <= year)
      TE <- metacont(n.e = total1[index.lower.year], mean.e = mean1[index.lower.year], 
                     sd.e = sd1[index.lower.year], n.c = total2[index.lower.year], 
                     mean.c = mean2[index.lower.year], sd.c = sd2[index.lower.year], 
                     studlab = study.name[index.lower.year])$pval.fixed #log of treatment effect of fixed effects m.a.
      results <- c(results, TE)
    }
  } 
  return(results)
}
# tmp <- data.cum %>% filter(meta.id == 163) %>% select(study.name, study.year, mean1, sd1, total1, mean2, sd2, total2, effect)
# tmp %>% mutate(cum.effect = cum.meta.cont.effect(study.year, study.name, mean1, sd1, total1, mean2, sd2, total2))

#Apply to dataset
cum.cont <- data.cum %>% filter(outcome.measure.new == "Mean Difference" | outcome.measure.new == "Std. Mean Difference" | 
  outcome.measure.new == "Hedges' G") %>% 
  mutate(zval.cont = cum.cont.zval.fixed(study.year, study.name, mean1, sd1, total1, mean2, sd2, total2))

cum.bin <- data.cum %>% filter(outcome.measure.new == "Risk Ratio" | outcome.measure.new == "Odds Ratio" | 
  outcome.measure.new == "Risk Difference") %>%
  mutate(zval.bin = cum.bin.zval.fixed(study.year, study.name, events1, events2, total1, total2, sm = "RR"))

cum.cont %>% ggplot(aes(x = study.year, y = abs(zval.cont))) +
  geom_smooth(method = "lm") + geom_point()

cum.cont %>% mutate(scaled.study.year = scale(study.year)) %>% filter(abs(zval.cont) < 20 & zval.cont) %>% ggplot(aes(x = scaled.study.year, y = abs(zval.cont))) +
  geom_smooth(method = "lm") + geom_point(size = .5)

cum.cont %>% mutate(scaled.study.year = scale(study.year)) %>% filter(abs(zval.cont) < 50 & zval.cont) %>% ggplot(aes(x = abs(zval.cont))) +
  geom_histogram()

cum.bin %>% mutate(scaled.study.year = scale(study.year)) %>% ggplot(aes(x = scaled.study.year, y = abs(zval.bin))) +
  geom_smooth(method = "lm") + geom_point()

cum.bin %>% ggplot(aes(x = (zval.bin))) +
  geom_histogram(bins = 100)

# print(cum.bin %>% group_by(zval.bin) %>% count() %>% summarize(max(n)),n = 200)
# 
# print(cum.cont %>% group_by(zval.cont) %>% count() %>% summarize(max(n)),n = 200)


m.cont <- lm(abs(zval.cont) ~ scale(study.year), cum.cont)
anova(m.cont)
m2.cont <- lm(abs(zval.cont) ~ scale(study.year) * scale(study.year, scale = F), cum.cont)
anova(m.cont, m2.cont)

m.bin <- lm(abs(zval.bin) ~ scale(study.year), cum.bin)
anova(m.bin)
m2.bin <- lm(abs(zval.bin) ~ scale(study.year) * scale(study.year, scale = F), cum.bin)
anova(m.bin, m2.bin)

