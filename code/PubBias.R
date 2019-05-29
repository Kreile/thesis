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
  meta.bin <- meta.bin.complete(data, min.study.number = 10, sig.level = 0.05)
  save(meta.bin, file =  file.path(PATH_RESULTS, file.bin))
}

file.cont <- "pb.cont.RData"
if (file.exists(file.path(PATH_RESULTS, file.cont))) {
  load(file.path(PATH_RESULTS, file.cont))
} else {
  meta.cont <- meta.cont.complete(data, min.study.number = 10, sig.level = 0.05)
  save(meta.cont, file =  file.path(PATH_RESULTS, file.cont))
}

file.meta <- "meta.RData"
if (file.exists(file.path(PATH_RESULTS, file.meta))) {
  load(file.path(PATH_RESULTS, file.meta))
} else {
  meta <- pb.meta.merge(meta.bin, meta.cont)
  save(meta, file =  file.path(PATH_RESULTS, file.meta))
}



load(file.path(PATH_RESULTS, file = "mly.RData"))
load(file.path(PATH_RESULTS, file = "data.processed.RData"))

require(biostatUZH)
require(tidyverse)
require(meta)
require(metasens)
require(gridExtra)
require(xtable)

# meta.cont <- meta.cont %>% mutate(sig.change.ranef.reg = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.reg.ranef.cont),
#                                             sig.change.ranef.copas = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.copas.cont),
#                                             sig.change.ranef.trimfill = sig.change(sig.before = sig.ranef.cont, sig.after =  sig.trimfill.cont),
#                                             sig.change.fixef.reg = sig.change(sig.before = sig.fixef.cont, sig.after =  sig.reg.ranef.cont),
#                                             sig.change.fixef.copas = sig.change(sig.before = sig.fixef.cont, sig.after =  sig.copas.cont),
#                                             sig.change.fixef.trimfill = sig.change(sig.before = sig.fixef.cont, sig.after =  sig.trimfill.cont))
# 
# meta.bin <- meta.bin %>% mutate(sig.copas = ifelse(pval.copas > 0.05, 0, 1),
#   sig.change.ranef.reg = sig.change(sig.before = sig.ranef, sig.after =  sig.reg.ranef), 
#                                                 sig.change.ranef.copas = sig.change(sig.before = sig.ranef, sig.after =  sig.copas),
#                                                 sig.change.ranef.trimfill = sig.change(sig.before = sig.ranef, sig.after =  sig.trimfill),
#                                                 sig.change.fixef.reg = sig.change(sig.before = sig.fixef, sig.after =  sig.reg.ranef),
#                                                 sig.change.fixef.copas = sig.change(sig.before = sig.fixef, sig.after =  sig.copas),
#                                                 sig.change.fixef.trimfill = sig.change(sig.before = sig.fixef, sig.after =  sig.trimfill))


# save(meta.cont, file =  file.path(PATH_RESULTS, file.cont))
# 
# file.bin = "pb.bin.RData"
# save(meta.bin, file =  file.path(PATH_RESULTS, file = "pb.bin.RData"))

     
     
meta <- pb.meta.merge(meta.bin, meta.cont)

require(biostatUZH)
require(tidyverse)
require(meta)
require(metasens)
require(gridExtra)
require(xtable)



########################################################################################################################################
########################################################################################################################################
#Significance of pbbias based on variance of total ss:
meta %>% filter(var.samplesize < 60000) %>% ggplot(aes(stat(count), x = var.samplesize, fill = factor(thomson.test))) + 
  geom_density(position = "fill")



## Publication test results: Binary.
test.bin <- meta.bin %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                          schwarzer.test = mean(schwarzer.test),
                                                          rucker.test = mean(rucker.test),
                                                          harbord.test = mean(harbord.test),
                                                          peter.test = mean(peter.test))
test.bin <- test.bin %>% gather(key = "test.type", value = "mean")

p.bin <- meta.bin %>% filter(sig.fixef.bin == 1) %>% ungroup() %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>%  
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Significant Pooled Effects")+
  annotate("text", x = test.bin$test.type, y = 200, 
           label = paste(round(test.bin$mean, 2)*100, "% rejected"), 
           color = "white", size = 5)



test.sig.bin <- meta.bin %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                              schwarzer.test = mean(schwarzer.test),
                                                              rucker.test = mean(rucker.test),
                                                              harbord.test = mean(harbord.test),
                                                              peter.test = mean(peter.test))
test.sig.bin <- test.sig.bin %>% gather(key = "test.type", value = "mean")

p1 <- meta.bin %>% ungroup() %>% 
select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + xlab(label = NULL) +
  annotate("text", x = test.sig.bin$test.type, y = 200, 
           label = paste(round(test.sig.bin$mean, 2)*100, "% rejected"), 
           color = "white", size = 5)

test.nonsig.bin <- meta.bin %>% filter(sig.fixef.bin == 0) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                                                                schwarzer.test = mean(schwarzer.test),
                                                                                                rucker.test = mean(rucker.test),
                                                                                                harbord.test = mean(harbord.test),
                                                                                                peter.test = mean(peter.test))
test.nonsig.bin <- test.nonsig.bin %>% gather(key = "test.type", value = "mean")

p2 <- meta.bin %>% 
  filter(sig.fixef.bin == 0) %>% ungroup() %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Non-significant Pooled Effects")+
  annotate("text", x = test.nonsig.bin$test.type, y = 200, 
           label = paste(round(test.nonsig.bin$mean, 2)*100, "% rejected"), 
           color = "white", size = 5)


grid.arrange(p1, p2, ncol = 1)

#Test Results: Continuous
test.cont <- meta.cont %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                            begg.test = mean(begg.test),
                                                            thomson.test = mean(thomson.test))

test.cont <- test.cont %>% gather(key = "test.type", value = "mean")

p.cont <- meta.cont %>% ungroup() %>% 
  select(egger.test, thomson.test, begg.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Significant Pooled Effects")+
  annotate("text", x = test.cont$test.type, y = 200, 
           label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
           color = "white", size = 5)


test.sig.cont <- meta.cont %>% filter(sig.fixef.cont == 1) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                                                                begg.test = mean(begg.test),
                                                                                                thomson.test = mean(thomson.test))

test.sig.cont <- test.sig.cont %>% gather(key = "test.type", value = "mean")

p3 <- meta.cont %>% 
  filter(sig.fixef.cont == 1) %>% ungroup() %>% 
  select(egger.test, thomson.test, begg.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Significant Pooled Effects")+
  annotate("text", x = test.sig.cont$test.type, y = 200, 
           label = paste(round(test.sig.cont$mean, 2)*100, "% rejected"), 
           color = "white", size = 5)

test.nonsig.cont <- meta.cont %>% filter(sig.fixef.cont == 0) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                                                                   begg.test = mean(begg.test),
                                                                                                   thomson.test = mean(thomson.test))

test.nonsig.cont <- test.nonsig.cont %>% gather(key = "test.type", value = "mean")

p4 <- meta.cont %>% 
  filter(sig.fixef.cont == 0) %>% ungroup() %>% 
  select(egger.test, thomson.test, begg.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Non-significant Pooled Effects")+
  annotate("text", x = test.nonsig.cont$test.type, y = 100, 
           label = paste(round(test.nonsig.cont$mean, 2)*100, "% rejected"), 
           color = "white", size = 5)

grid.arrange(p3, p4, ncol = 1)

grid.arrange(p.bin, p.cont, ncol = 1)

########################################################################################################################################
########################################################################################################################################
# Exploratory Data analysis
########################################################################################################################################
########################################################################################################################################

#Mean and Total study sample size
meta %>% filter(mean.samplesize < 1500) %>%  ggplot(aes(x = mean.samplesize)) + geom_histogram()
meta %>% filter(total.samplesize < 50000) %>% ggplot(aes(x = total.samplesize)) + geom_histogram()
meta %>% filter(total.samplesize < 1000) %>% ggplot(aes(x = total.samplesize, fill = factor(thomson.test))) + 
  geom_histogram()
meta %>% filter(mean.samplesize < 100) %>% ggplot(aes(x = mean.samplesize, fill = factor(thomson.test))) + 
  geom_histogram()

#Meta analysis sample size
meta %>% filter(n < 100) %>% ggplot(aes(x = n, fill = factor(thomson.test))) + 
  geom_histogram()
meta %>% filter(n < 40) %>% ggplot(aes(stat(count), x = n, fill = factor(thomson.test))) + 
  geom_density(position = "fill")

meta %>% filter(n < 100) %>% ggplot(aes(x = n, fill = factor(sig.fixef))) + 
  geom_histogram()
meta %>% filter(n < 40) %>% ggplot(aes(stat(count), x = n, fill = factor(sig.fixef))) + 
  geom_density(position = "fill")

#Dependence of significance of publication bias tests on mean and total samplesize
meta %>% filter(mean.samplesize < 750) %>% ggplot(aes(stat(count), x = mean.samplesize, fill = factor(thomson.test))) + 
  geom_density(position = "fill")

meta %>% filter(total.samplesize < 10000) %>% ggplot(aes(stat(count), x = total.samplesize, fill = factor(thomson.test))) + 
  geom_density(position = "fill")

#Change in significance
meta %>% filter(!is.na(sig.change)) %>% ggplot(aes(x = sig.change, stat = "count")) + geom_histogram(stat = "count")

#Trimfill proportion and test results:
trimfill.cont.mean <- meta.cont %>% ungroup() %>% summarise(mean = mean(missing.trim.cont)) %>% select(mean)
trimfill.cont.median <- meta.cont %>% ungroup() %>% summarise(median = median(missing.trim.cont)) %>% select(median)

trimfill.bin.mean <- meta.bin %>% ungroup() %>% summarise(mean = mean(missing.trim)) %>% select(mean)
trimfill.bin.median <- meta.bin %>% ungroup() %>% summarise(median = median(missing.trim)) %>% select(median)

meta.bin %>% ggplot(aes(x = missing.trim)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Fraction of Trimmed Comparisons (binary outcome)") + xlab("Fraction") + ylab("Frequency")

meta.cont %>% ggplot(aes(x = missing.trim.cont)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Fraction of Trimmed Comparisons (continuous outcome)") + xlab("Fraction") + ylab("Frequency")

meta.bin %>% ggplot(aes(x = missing.trim, fill = factor(peter.test), stat(count))) + 
  geom_density(na.rm = T, position = "fill") + 
  labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))

meta.bin %>% ggplot(aes(x = missing.trim, fill = factor(harbord.test), stat(count))) + geom_density(na.rm = T, position = "fill")

meta.cont %>% filter(missing.trim.cont < 0.5) %>% 
  ggplot(aes(x = missing.trim.cont, fill = factor(thomson.test), stat(count))) + 
  geom_density(na.rm = T, position = "fill") + 
  labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))

#Time trends in significance of tests
meta.bin %>% filter(mean.publication.year < 2019 & mean.publication.year > 1960) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef))) + geom_histogram()

meta.bin %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef), stat(count))) + geom_density(na.rm = T, position = "fill")

meta.bin %>% filter(first.publication.year < 2013 & first.publication.year > 1990) %>% 
  ggplot(aes(x = first.publication.year, fill = factor(sig.fixef))) + geom_density(na.rm = T, position = "fill")

meta.bin %>% filter(first.publication.year < 2013 & first.publication.year > 1990) %>% 
  ggplot(aes(x = first.publication.year, fill = factor(peter.test))) + geom_density(na.rm = T, position = "fill")

meta.bin %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(rucker.test))) + geom_density(na.rm = T, position = "fill")


meta.cont %>% filter(mean.publication.year < 2019 & mean.publication.year > 1960) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef.cont))) + geom_histogram()

meta.cont %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef.cont))) + geom_density(na.rm = T, position = "fill")

meta.cont %>% filter(first.publication.year < 2013 & first.publication.year > 1990) %>% 
  ggplot(aes(x = first.publication.year, fill = factor(sig.fixef.cont))) + geom_density(na.rm = T, position = "fill")

meta.cont %>% filter(first.publication.year < 2013 & first.publication.year > 1990) %>% 
  ggplot(aes(x = first.publication.year, fill = factor(thomson.test))) + geom_density(na.rm = T, position = "fill")

meta.cont %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(egger.test))) + geom_density(na.rm = T, position = "fill")


meta %>% filter(mean.publication.year < 2019 & mean.publication.year > 1960) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef))) + geom_histogram() + 
  labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))

meta %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef))) + geom_density(na.rm = T, position = "identity") + 
  labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))

meta %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(thomson.test), stat(count))) + geom_density(na.rm = T, position = "fill") + 
  labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))

cdplot(x = meta$mean.publication.year[!is.na(meta$mean.publication.year)], y = factor(meta$thomson.test[!is.na(meta$mean.publication.year)]))

meta %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
  ggplot(aes(x = mean.publication.year, fill = factor(egger.test))) + geom_density(na.rm = T, position = "fill")

#Sample size median over the years
data %>% distinct(study.name, .keep_all = T) %>%  filter(study.year < 2018 & study.year > 1980) %>% group_by(study.year) %>% 
  summarize(samplesize = median(total1 + total2, na.rm = T)) %>% 
  ggplot(aes(y = samplesize, x = study.year)) + geom_line()


########################################################################################################################################
########################################################################################################################################
# Pubbias test statistics 
########################################################################################################################################
########################################################################################################################################
#Agreement proportions of publication bias tests:

#Binary:
meta.bin <- meta.bin %>%ungroup() %>%  mutate(egger.schwarzer = ifelse(egger.test == schwarzer.test, "agree", "disagree"),
                                              egger.peter = ifelse(egger.test == peter.test, "agree", "disagree"),
                                              egger.rucker = ifelse(egger.test == rucker.test, "agree", "disagree"),
                                              egger.harbord = ifelse(egger.test == harbord.test, "agree", "disagree"),
                                              schwarzer.peter = ifelse(schwarzer.test == peter.test, "agree", "disagree"),
                                              schwarzer.rucker = ifelse(schwarzer.test == rucker.test, "agree", "disagree"),
                                              schwarzer.harbord = ifelse(schwarzer.test == harbord.test, "agree", "disagree"),
                                              rucker.peter = ifelse(rucker.test == peter.test, "agree", "disagree"),
                                              rucker.harbord = ifelse(rucker.test == harbord.test, "agree", "disagree"),
                                              harbord.peter = ifelse(harbord.test == peter.test, "agree", "disagree"))

agreement.bin <- meta.bin %>% ungroup() %>% summarise(egger.schwarzer = sum(egger.schwarzer == "agree")/n(),
                                                      egger.peter = sum(egger.peter == "agree")/n(),
                                                      egger.rucker = sum(egger.rucker == "agree")/n(),
                                                      egger.harbord = sum(egger.harbord == "agree")/n(),
                                                      schwarzer.peter = sum(schwarzer.peter == "agree")/n(),
                                                      schwarzer.rucker = sum(schwarzer.rucker == "agree")/n(),
                                                      schwarzer.harbord = sum(schwarzer.harbord == "agree")/n(),
                                                      rucker.peter = sum(rucker.peter == "agree")/n(),
                                                      harbord.peter = sum(harbord.peter == "agree")/n())

correlation.bin <- meta.bin %>%ungroup() %>%  summarise(egger.schwarzer = cor(pval.egger, pval.schwarzer),
                                                        egger.peter = cor(pval.egger, pval.peter),
                                                        egger.rucker = cor(pval.egger, pval.rucker),
                                                        egger.harbord = cor(pval.egger, pval.harbord),
                                                        schwarzer.peter = cor(pval.schwarzer, pval.peter),
                                                        schwarzer.rucker = cor(pval.schwarzer, pval.rucker),
                                                        schwarzer.harbord = cor(pval.schwarzer, pval.harbord),
                                                        rucker.peter = cor(pval.rucker, pval.peter),
                                                        harbord.peter = cor(pval.harbord, pval.peter))

binary.tests.agreement <- rbind(agreement.bin, correlation.bin)
rownames(binary.tests.agreement) <- c("Test Agreement","P-value Correlation")


#Continuous:

meta.cont <- meta.cont %>% ungroup() %>% mutate(thomson.egger = ifelse(thomson.test == egger.test, "agree", "disagree"),
                                                thomson.begg = ifelse(thomson.test == begg.test, "agree", "disagree"),
                                                egger.begg = ifelse(egger.test == begg.test, "agree", "disagree"))

agreement.cont <- meta.cont %>% ungroup() %>%  summarise(thomson.egger = sum(thomson.egger == "agree")/n(),
                                                         thomson.begg = sum(thomson.begg == "agree")/n(),
                                                         egger.begg = sum(egger.begg == "agree")/n())

correlation.cont <- meta.cont %>% ungroup() %>% summarise(thomson.egger = cor(pval.thomson.cont, pval.egger.cont),
                                                          thomson.begg = cor(pval.thomson.cont, pval.begg.cont),
                                                          egger.begg = cor(pval.egger.cont, pval.begg.cont))

cont.tests.agreement <- rbind(agreement.cont, correlation.cont)
rownames(cont.tests.agreement) <- c("Test Agreement","P-value Correlation")

egger.schwarzer <- meta.bin %>% group_by(egger.schwarzer) %>% count() %>% ggplot(aes(x = egger.schwarzer, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.peter <- meta.bin %>% group_by(egger.peter) %>% count() %>% ggplot(aes(x = egger.peter, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.rucker <- meta.bin %>% group_by(egger.rucker) %>% count() %>% ggplot(aes(x = egger.rucker, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.harbord <- meta.bin %>% group_by(egger.harbord) %>% count() %>% ggplot(aes(x = egger.harbord, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.rucker <- meta.bin %>% group_by(egger.rucker) %>% count() %>% ggplot(aes(x = egger.rucker, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.rucker <- meta.bin %>% group_by(egger.rucker) %>% count() %>% ggplot(aes(x = egger.rucker, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")
egger.rucker <- meta.bin %>% group_by(egger.rucker) %>% count() %>% ggplot(aes(x = egger.rucker, y = nn)) +
  geom_col(col = "gray15", fill = "dodgerblue")

library(grid)
library(gridExtra)
grid.arrange(egger.schwarzer, egger.peter, egger.rucker)





########################################################################################################################################
########################################################################################################################################
# Training dataset
########################################################################################################################################
########################################################################################################################################



meta1 <- data %>% filter(file.nr == 8 & comparison.nr == 1 & outcome.nr == 1 & subgroup.nr == 1)
meta2 <- data %>% filter(file.nr == 76 & comparison.nr == 2 & outcome.nr == 2 & subgroup.nr == 4)
meta3 <- data %>% filter(file.nr == 953 & comparison.nr == 2 & outcome.nr == 4 & subgroup.nr == 1)
meta4 <- data %>% filter(file.nr == 944 & comparison.nr == 12 & outcome.nr == 6 & subgroup.nr == 4)
meta5 <- data %>% filter(file.nr == 544 & comparison.nr == 1 & outcome.nr == 1 & subgroup.nr == 3)

metay.1 <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", data = meta1)
metay.2 <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", data = meta2)
metay.3 <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", data = meta3)
metay.4 <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", data = meta4)
metay.5 <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, sm= "OR", data = meta5)

meta.g <- rbind(meta1, meta2, meta3, meta4, meta5)
meta.g.sum <- meta.bin.complete(meta.g, min.study.number = 10, sig.level = 0.05)

listm = list(m1 = metay.1, m2 = metay.2, m3 = metay.3, m4 = metay.4, m5 = metay.5)

for(u in 1:5){
  funnel(listm[[u]], main = paste("m", u))
}


########################################################################################################################################
########################################################################################################################################
#########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
#########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
#########################################################################################################################################
#Meta-analysis Functions

#Copas selection model automatic estimate and std. error and N.unpubl extraction:
#The estimate with smallest N.unpubl and a p-value larger than sig.level + 0.05 is chosen.
auto.copas <- function(meta.obj, sig.level){
  sig.level <- sig.level + 0.05
  gamma0 <- -1.7 #analog to P(select|small trial w. sd = 0.4) = 0.1 and P(select|large trial w. sd  = 0.05) = 0.9
  gamma1 <- 0.16 #from limitmeta paper (Rücker 2011): "small range" procedure - if no nonsignificance - "broad range"
  copas <- copas(meta.obj, gamma0.range = c(gamma0, 2), gamma1.range = c(0, gamma1))
  pval.rsb <- copas$pval.rsb
  N.unpubl <- copas$N.unpubl
  if(all(pval.rsb < sig.level)){
    copas <- copas(meta.obj, , gamma0.range = c(2*gamma0 - 2, 2), gamma1.range = c(0, 2*gamma1))
    pval.rsb <- copas$pval.rsb
    N.unpubl <- copas$N.unpubl
    if(all(pval.rsb < sig.level)){
      corr.est <- NA
      se.corr.est <- NA
      N.unpubl <- NA
    } else{ 
      ind.nonsig <- which(pval.rsb > sig.level)
      ind.estimate <- which.min(N.unpubl[ind.nonsig])
      corr.est <- copas$TE.slope[ind.nonsig[ind.estimate]]
      se.corr.est <- copas$seTE.slope[ind.nonsig[ind.estimate]]
      N.unpubl <- copas$N.unpubl[ind.nonsig[ind.estimate]]
    }
  } else{
    if(all(pval.rsb > sig.level)){
      corr.est <- NA
      se.corr.est <- NA
      N.unpubl <- NA
    } else{
      ind.nonsig <- which(pval.rsb > sig.level)
      ind.estimate <- which.min(N.unpubl[ind.nonsig])
      corr.est <- copas$TE.slope[ind.nonsig[ind.estimate]]
      se.corr.est <- copas$seTE.slope[ind.nonsig[ind.estimate]]
      N.unpubl <- copas$N.unpubl[ind.nonsig[ind.estimate]]
    }
    return(c(corr.est, se.corr.est, N.unpubl))
  }
}


#Function to classify how significance changed after correction
sig.change <- function(sig.before, sig.after){
  if(!is.na(sig.before & !is.na(sig.after))){
    if(sig.before == 0){
      results <- ifelse(sig.after == 0, 0, 1)
    } else{
      results <- ifelse(sig.after == 1, 2, 3)}
  } else{
    results <- NA
  }
  return(results)
}



meta.bin.complete <- function(data, min.study.number, sig.level){
  metadat <- data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>%  
    filter(file.nr != 3014) %>% #file.nr 3014 gibt eine seltsame Fehlermeldung, muss ausgeschlossen werden
    filter(events1 > 0 | events2 > 0) %>% #comparisons mit events  = 0 ausschliessen
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>% #comparisons mit events = total ausschliessen - vertraegt der alg. irgendwie nicht.
    group_by(file.nr, comparison.nr, outcome.measure.new, subgroup.nr) %>% #Die nachfolgenden Rechnungen werden jeweils fuer diese Gruppen gemacht.
    mutate(n = n()) %>% filter(n > 9) %>% #Nur meta-analysen mit mehr als 9 comparisons behalten
    mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
    mutate(study.year = ifelse(study.year > 1950, study.year, NA)) %>% 
    summarize(doi = unique(doi), #doi
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
  return(metadat)
} 


meta.cont.complete <- function(data, min.study.number, sig.level){
  metadat <- data %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>%
    filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
    filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
    mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
    mutate(study.year = ifelse(study.year > 1950, study.year, NA)) %>% 
    group_by(file.nr, comparison.nr, outcome.mesaure.new, subgroup.nr) %>% 
    mutate(n = n()) %>% filter(n >= min.study.number) %>% 
    summarize(doi = unique(doi), #doi
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
  return(metadat)
}

########################################################################################################################################
########################################################################################################################################
#########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
#########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
#########################################################################################################################################


meta.bin <- meta.bin.complete(data, min.study.number = 10, sig.level = 0.05)
save(meta.bin, file = "PubBias_results3.RData")

meta.cont <- meta.cont.complete(data, min.study.number = 10, sig.level = 0.05)
save(meta.cont, file = "PubBias_results4.RData", file.p)
save(meta.cont, file.path(PATH_RESULTS, 'PubBias_results1.RData'))

meta <- bind_rows(meta.bin, meta.cont)
meta <- meta %>% mutate(pval.thomson = ifelse(!is.na(pval.rucker.bin), pval.rucker.bin, pval.thomson.cont),
                        thomson.test = ifelse(!is.na(thomson.test), thomson.test, rucker.test),
                        sig.ranef = ifelse(!is.na(sig.ranef.bin), sig.ranef.bin, sig.ranef.cont),
                        sig.fixef = ifelse(!is.na(sig.fixef.bin), sig.fixef.bin, sig.fixef.cont),
                        pval.fixef = ifelse(!is.na(pval.fixef.bin), pval.fixef.bin, pval.fixef.cont),
                        pval.ranef = ifelse(!is.na(pval.ranef.bin), pval.ranef.bin, pval.ranef.cont),
                        pval.copas = ifelse(!is.na(pval.copas.bin), pval.copas.bin, pval.copas.cont),
                        pval.reg.ranef = ifelse(!is.na(pval.reg.ranef.bin), pval.reg.ranef.bin, pval.reg.ranef.cont),
                        pval.trimfill = ifelse(!is.na(pval.timfill.bin), pval.trimfill.bin, pval.trimfill.cont),
                        missing.trim = ifelse(!is.na(missing.trim.bin), missing.trim.bin, missing.trim.cont),
                        missing.copas = ifelse(!is.na(missing.copas.bin), missing.copas.bin, missing.copas.cont),
                        sig.copas = ifelse(!is.na(sig.copas.bin), sig.copas.bin, sig.copas.cont),
                        sig.trimfill = ifelse(!is.na(sig.trimfill.bin), sig.trimfill.bin, sig.trimfill.cont),
                        sig.reg.ranef = ifelse(!is.na(sig.reg.ranef.bin), sig.reg.ranef.bin, sig.reg.ranef.cont),
                        Q = ifelse(!is.na(Q.bin), Q.bin, Q.cont),
                        pval.Q = ifelse(!is.na(pval.Q.bin), pval.Q.bin, pval.Q.cont),
                        sig.Q = ifelse(!is.na(sig.Q.bin), sig.Q.bin, sig.Q.cont))
# = ifelse(!is.na(), .bin, .cont),
# = ifelse(!is.na(), .bin, .cont),
# = ifelse(!is.na(), .bin, .cont),
# = ifelse(!is.na(), .bin, .cont),
# = ifelse(!is.na(), .bin, .cont),
# val.fixef = ifelse(!is.na(pval.fixef.bin), pval.fixef.bin, pval.fixef.cont),)


