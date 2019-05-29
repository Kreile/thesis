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





#New Data process function

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

if (file.exists(file.path(PATH_RESULTS, "data.RData"))) {
  load(file.path(PATH_RESULTS, "data.RData"))
} else {
  data = pb.readData(path = PATH_DATA, file = FILE)
  tmp = pb.clean(data)
  data = tmp[[1]]
  aliases = tmp[[2]]
  save(data, file =  file.path(PATH_RESULTS, "data.RData"))
}

if (file.exists(file.path(PATH_RESULTS, "data.ext.RData"))) {
  load(file.path(PATH_RESULTS, "data.ext.RData"))
} else {
  data.ext <- pb.process(data)
  save(data.ext, file =  file.path(PATH_RESULTS, "data.ext.RData"))
}



#Meta-analysis function for any type "bin", "cont" and "surv":
#data has to be processed by pb.readData() and cleaned by pb.clean(), 
#processed by pb.process() (all functions in global environment from "PubBias_functions.R" script)

mtd <- data.ext %>% filter(meta.id == 101436)

meta.cor <- metacor(cor = fishersz, n = total1 + total2, studlab = study.name, sm = "ZCOR", 
                    comb.fixed = TRUE, method.tau = "REML", null.effect = 0, data = mtd)
meta.bin <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, method = "Inverse", 
                    method.tau = "REML", data = mtd, sm = "RR")

forest(meta.cor)
forest(meta.bin)

mtd.c <- data.ext %>% filter(meta.id == 122943)

meta.cor.c <- metacor(cor = fishersz, n = total1 + total2, studlab = study.name, sm = "ZCOR", 
                      comb.fixed = TRUE, method.tau = "REML", null.effect = 0, data = mtd.c)

meta.cont.c <- metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", 
                        studlab = study.name, data = mtd.c)

forest(meta.cor.c)
forest(meta.cont.c)

mtd %>%
  mutate(log.OR = log( (events1 * (total2 - events2))/ (events2 * (total1 - events1)) ),
         se.log.OR = sqrt( (1/events1) + 1/(total1 - events1) + (1/events2) + 1/(total2 - events2)),
         std.mean.d = log.OR * (sqrt(3)/pi),
         correlation = std.mean.d/(sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2 *total1))),
         fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ),
         fishersz.variance = 1/(total1 + total2 - 3),
         pval.type = 2*(1- pnorm(abs(log.OR/se.log.OR)))) %>% select(correlation)
