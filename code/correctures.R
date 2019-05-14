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

load(file.path(PATH_RESULTS, file = "mly.RData"))
load(file.path(PATH_RESULTS, file = "data.processed.RData"))

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


require(biostatUZH)
require(tidyverse)
require(meta)
require(metasens)
require(gridExtra)


#Meta filtering: 
#Meta filtering: 
metac.bin <- meta.bin %>% filter(n.sig.type.bin > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac.cont <- meta.cont %>% filter(n.sig.type.cont > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac <- meta %>% filter(n.sig.type > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)

#Publication Bias Test Agreement:
meta.bin %>% mutate(n.sig = peter.test + rucker.test + egger.test + harbord.test + schwarzer.test) %>% 
  group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
  ggplot(aes(y = nn, x = n.sig)) + geom_col() 
meta.cont %>% mutate(n.sig = egger.test + thomson.test + begg.test) %>% 
  group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
  ggplot(aes(y = nn, x = n.sig)) + geom_col() 


########################################################################################################################
########################################################################################################################
#Check if proportions of pbbias change if only meta-analyses from different reviews are analyzed:
########################################################################################################################
########################################################################################################################

#Overall test:
test.bin <- meta.bin %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                 schwarzer.test = mean(schwarzer.test),
                                                 rucker.test = mean(rucker.test),
                                                 harbord.test = mean(harbord.test),
                                                 peter.test = mean(peter.test))
test.bin <- test.bin %>% gather(key = "test.type", value = "mean")

p.bin <- meta.bin %>% ungroup() %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>%  
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + xlab(label = NULL) +
  annotate("text", x = test.bin$test.type, y = 1000, 
           label = paste(round(test.bin$mean, 2)*100, "% rejected"), 
           color = "white")



test.sig.bin <- meta.bin %>% filter(sig.fixef.bin == 1) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                                                    schwarzer.test = mean(schwarzer.test),
                                                                                    rucker.test = mean(rucker.test),
                                                                                    harbord.test = mean(harbord.test),
                                                                                    peter.test = mean(peter.test))
test.sig.bin <- test.sig.bin %>% gather(key = "test.type", value = "mean")

p1 <- meta.bin %>% ungroup() %>% filter(sig.fixef.bin == 1) %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
  annotate("text", x = test.sig.bin$test.type, y = 1750, 
           label = paste(round(test.sig.bin$mean, 2)*100, "% rejected"), 
           color = "white")

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
  theme_bw() + ggtitle("Non-significant Pooled Effects")+ xlab(label = NULL) +theme(legend.position="none") +
  annotate("text", x = test.nonsig.bin$test.type, y = 650, 
           label = paste(round(test.nonsig.bin$mean, 2)*100, "% rejected"), 
           color = "white")

range.pb.difference.bin <- range(test.sig.bin$mean - test.nonsig.bin$mean)

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
  theme_bw() + xlab(label = NULL) +
  annotate("text", x = test.cont$test.type, y = 750, 
           label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
           color = "white")


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
  theme_bw() + ggtitle("Significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
  annotate("text", x = test.sig.cont$test.type, y = 600, 
           label = paste(round(test.sig.cont$mean, 2)*100, "% rejected"), 
           color = "white")

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
  theme_bw() + ggtitle("Non-significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
  annotate("text", x = test.nonsig.cont$test.type, y = 100, 
           label = paste(round(test.nonsig.cont$mean, 2)*100, "% rejected"), 
           color = "white")

########################################################################################################################
########################################################################################################################
#Selective test
########################################################################################################################
########################################################################################################################


meta.bin <- meta.bin %>% distinct(file.nr, .keep_all = T)
meta.cont <- meta.cont %>% distinct(file.nr, .keep_all = T)

test.bin <- meta.bin %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                 schwarzer.test = mean(schwarzer.test),
                                                 rucker.test = mean(rucker.test),
                                                 harbord.test = mean(harbord.test),
                                                 peter.test = mean(peter.test))
test.bin <- test.bin %>% gather(key = "test.type", value = "mean")

p.bins <- meta.bin %>% ungroup() %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>%  
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + xlab(label = NULL) +
  annotate("text", x = test.bin$test.type, y = 1000, 
           label = paste(round(test.bin$mean, 2)*100, "% rejected"), 
           color = "white")



test.sig.bin <- meta.bin %>% filter(sig.fixef.bin == 1) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                                                    schwarzer.test = mean(schwarzer.test),
                                                                                    rucker.test = mean(rucker.test),
                                                                                    harbord.test = mean(harbord.test),
                                                                                    peter.test = mean(peter.test))
test.sig.bin <- test.sig.bin %>% gather(key = "test.type", value = "mean")

p1s <- meta.bin %>% ungroup() %>% filter(sig.fixef.bin == 1) %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
  annotate("text", x = test.sig.bin$test.type, y = 1750, 
           label = paste(round(test.sig.bin$mean, 2)*100, "% rejected"), 
           color = "white")

test.nonsig.bin <- meta.bin %>% filter(sig.fixef.bin == 0) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                                                       schwarzer.test = mean(schwarzer.test),
                                                                                       rucker.test = mean(rucker.test),
                                                                                       harbord.test = mean(harbord.test),
                                                                                       peter.test = mean(peter.test))
test.nonsig.bin <- test.nonsig.bin %>% gather(key = "test.type", value = "mean")

p2s <- meta.bin %>% 
  filter(sig.fixef.bin == 0) %>% ungroup() %>% 
  select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Non-significant Pooled Effects")+ xlab(label = NULL) +theme(legend.position="none") +
  annotate("text", x = test.nonsig.bin$test.type, y = 650, 
           label = paste(round(test.nonsig.bin$mean, 2)*100, "% rejected"), 
           color = "white")

range.pb.difference.bin <- range(test.sig.bin$mean - test.nonsig.bin$mean)

#Test Results: Continuous
test.cont <- meta.cont %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                   begg.test = mean(begg.test),
                                                   thomson.test = mean(thomson.test))

test.cont <- test.cont %>% gather(key = "test.type", value = "mean")

p.conts <- meta.cont %>% ungroup() %>% 
  select(egger.test, thomson.test, begg.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + xlab(label = NULL) +
  annotate("text", x = test.cont$test.type, y = 150, 
           label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
           color = "white")


test.sig.cont <- meta.cont %>% filter(sig.fixef.cont == 1) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                                                       begg.test = mean(begg.test),
                                                                                       thomson.test = mean(thomson.test))

test.sig.cont <- test.sig.cont %>% gather(key = "test.type", value = "mean")

p3s <- meta.cont %>% 
  filter(sig.fixef.cont == 1) %>% ungroup() %>% 
  select(egger.test, thomson.test, begg.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
  annotate("text", x = test.sig.cont$test.type, y = 100, 
           label = paste(round(test.sig.cont$mean, 2)*100, "% rejected"), 
           color = "white")

test.nonsig.cont <- meta.cont %>% filter(sig.fixef.cont == 0) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
                                                                                          begg.test = mean(begg.test),
                                                                                          thomson.test = mean(thomson.test))

test.nonsig.cont <- test.nonsig.cont %>% gather(key = "test.type", value = "mean")

p4s <- meta.cont %>% 
  filter(sig.fixef.cont == 0) %>% ungroup() %>% 
  select(egger.test, thomson.test, begg.test) %>% 
  gather(key = "test.type", value = "null.hypothesis") %>% 
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + ggtitle("Non-significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
  annotate("text", x = test.nonsig.cont$test.type, y = 100, 
           label = paste(round(test.nonsig.cont$mean, 2)*100, "% rejected"), 
           color = "white")

grid.arrange(p.bin, p.bins, ncol = 2)
grid.arrange(p.cont, p.conts, ncol = 2)
grid.arrange(p1, p1s, ncol = 2)
grid.arrange(p2, p2s, ncol = 2)
grid.arrange(p3, p3s, ncol = 2)
grid.arrange(p4, p4s, ncol = 2)
