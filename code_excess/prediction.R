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

yodata <- yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  mutate(counts = n()) %>% filter(counts > 2)

#Fishersz transformation applied has probably rounding errors.. :

#Random effects prediction intervals, empirical vs. nominal coverage
loo.metacor.upr <- function(cor, n, lev){
  result <- c()
  for ( i in seq_along(cor) ) {
    result[i] <- metacor(cor = cor[-i], n = n[-i], sm = "COR", prediction = T, comb.random = T, level.predict = lev)$upper.predict
  }
  return(result)
}

loo.metacor.lwr <- function(cor, n, lev){
  result <- c()
  for ( i in seq_along(cor) ) {
    result[i] <- metacor(cor = cor[-i], n = n[-i], sm = "COR", prediction = T, comb.random = T, level.predict = lev)$lower.predict
  }
  return(result)
}

#Empirical vs nominal coverage comparison for fixed effect prediction intervals
#Loop:
levels <- seq(0.5, 0.99, by = 0.1)
e.cs <- c()
temp.loop <- yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  
  filter(file.nr < 500)
for(level in levels){
  empirical.coverage <-  temp.loop %>% 
    mutate(upr.pred = loo.metacor.upr(cor = fishersz, n = total1 + total2, lev = level),
           lwr.pred = loo.metacor.lwr(cor = fishersz, n = total1 + total2, lev = level)) %>% 
    filter(!is.na(lwr.pred) & !is.na(upr.pred)) %>% arrange(study.year) %>% 
    mutate(covered = ifelse(fishersz < upr.pred & fishersz > lwr.pred, 1, 0), id = row_number()) %>%
    filter(row_number() > 1) %>% 
    ungroup() %>% summarise(sum(covered)/length(covered))
  e.cs <- c(e.cs, empirical.coverage)
}

plot(levels, as.numeric(e.cs), ylab = "empirical coverage", xlab = "nominal coverage", ylim = c(.5, 1), xlim = c(.5,1))
lines(c(0,1),c(0,1), col = 11)

#Random effects prediction intervals, empirical vs. nominal coverage
loo.metacor.upr <- function(cor, n, lev){
  result <- c()
  for ( i in seq_along(cor) ) {
    result[i] <- metacor(cor = cor[-i], n = n[-i], sm = "COR", prediction = T, comb.random = T, level.predict = lev, hakn = T)$upper.predict
  }
  return(result)
}

loo.metacor.lwr <- function(cor, n, lev){
  result <- c()
  for ( i in seq_along(cor) ) {
    result[i] <- metacor(cor = cor[-i], n = n[-i], sm = "COR", prediction = T, comb.random = T, level.predict = lev, hakn = T)$lower.predict
  }
  return(result)
}

#Empirical vs nominal coverage comparison for fixed effect prediction intervals (Hartung-Knapp)
#Loop:
levels <- seq(0.5, 0.99, by = 0.1)
e.cs.h <- c()
for(level in levels){
  empirical.coverage <- yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  
    filter(file.nr < 100) %>% 
    mutate(upr.pred = loo.metacor.upr(cor = fishersz, n = total1 + total2, lev = level),
           lwr.pred = loo.metacor.lwr(cor = fishersz, n = total1 + total2, lev = level)) %>% 
    arrange(study.year) %>% filter(!is.na(lwr.pred) & !is.na(upr.pred)) %>% 
    mutate(covered = ifelse(fishersz < upr.pred & fishersz > lwr.pred, 1, 0), id = row_number()) %>%
    filter(row_number() > 1) %>% 
    ungroup() %>% summarise(sum(covered)/length(covered))
  e.cs.h <- c(e.cs.h, empirical.coverage)
}

plot(levels, as.numeric(e.cs.h), ylab = "empirical coverage", xlab = "nominal coverage", ylim = c(.5, 1), xlim = c(.5,1))
lines(c(0,1),c(0,1), col = 11)

#In some cases, metacor$lower.predict gives NA ->
# #Check1
# yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
#   filter(file.nr == 21) %>% 
#   mutate(
#     l.pred = loo.metacor.lwr(cor = fishersz, n = total1 + total2),
#     u.pred = loo.metacor.upr(cor = fishersz, n = total1 + total2)) %>% select(fishersz, total1, l.pred, u.pred)
# 
# yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  
#   filter(file.nr == 22) %>% 
#   mutate(
#     l.pred = loo.metacor.upr(cor = fishersz, n = total1 + total2),
#     u.pred = loo.metacor.lwr(cor = fishersz, n = total1 + total2)) %>% select(fishersz, total1, l.pred, u.pred)
# print(metacor(cor = c(0, 0.106), n = c(41, 37), sm = "COR", prediction = T))
# meta <- metacor(cor = c(0, 0.106), n = c(41, 37), sm = "COR", prediction = T) 
# tp <- yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  
#   filter(file.nr == 22, outcome.nr == 2, comparison.nr == 1) %>% select(fishersz, total1, total2)
# tp2 <- yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  
#   filter(file.nr == 21, outcome.nr == 2, comparison.nr == 1) %>% select(fishersz, total1, total2)
# 
# 
# loo.metacor.lwr(tp$fishersz, tp$total1 + tp$total2)
# loo.metacor.lwr(tp2$fishersz, tp2$total1 + tp2$total2)
# tp %>% mutate(n = total1 + total2) %>%  mutate(lwr = loo.metacor.lwr(cor = fishersz, n = n)) %>% select(lwr)
# tp %>% ungroup() %>% mutate(n = total1 + total2) %>%  mutate(lwr = loo.metacor.lwr(cor = fishersz, n = n)) %>% select(lwr)

