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

load(file.path(PATH_RESULTS, "meta_analyses_summary_complete.RData"))
load(file.path(PATH_RESULTS, "meta.bin.RData"))
load(file.path(PATH_RESULTS, "meta.cont.RData"))
load(file.path(PATH_RESULTS, "meta.surv.RData"))
load(file.path(PATH_RESULTS, "data_used_for_analysis.RData"))




#Binary data examples:
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



#Examplary m.a. for test purpose: 134944 (bin), 158355 (cont)

bin <- data.ext2 %>% filter(meta.id == 134944)
cont <- data.ext2 %>% filter(meta.id == 158355)
bin2 <- data.ext2 %>% filter(meta.id == 16060)


forest(meta.id(134944))
forest(metacor(cor = cor.phi, n = total1 + total2, data = bin))
forest(metacor(cor = cor.fisher, n = total1 + total2, data = bin))


forest(rma.uni(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2,  
							 measure = "RR", method = "FE", data = bin)) #For comparison
forest(rma.uni(yi = smd.pbit, vi = var.smd.pbit, measure = "PBIT", method = "FE", 
							 data = bin)) 
forest(rma.uni(yi = cor.phi, vi = var.cor.phi, measure = "PHI", method = "FE", 
							 data = bin)) 
forest(rma.uni(yi = cor.fisher, vi = var.cor.fisher, measure = "ZCOR", method = "FE", 
							 data = bin)) 

forest(rma.uni(yi = effect, sei = se, measure = "MD", method = "FE", data = cont))
forest(rma.uni(m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, n1i = total1, n2i = total2, 
							 measure = "MD", method = "FE", data = cont))
forest(rma.uni(yi = smd, vi = var.smd, measure = "SMD", method = "FE", data = cont))
forest(rma.uni(yi = cor.pearson, vi = var.cor.pearson, measure = "COR", method = "FE", data = cont))
forest(rma.uni(yi = cor.fisher, vi = var.cor.fisher, measure = "ZCOR", method = "FE", data = cont))

