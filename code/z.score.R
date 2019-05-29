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





####################################################################################
#Comparison of all adjusted z-scores:
sig.zcor <- meta.wcor %>% mutate(z.fixef = est.z.fixef/se.est.z.fixef,
																 z.ranef = est.z.ranef/se.est.z.ranef,
																 z.reg = est.z.reg/se.est.z.reg,
																 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	group_by(method) %>% 
	summarise(significant = length(which(abs(fisher.z) > 1.96)),
						p.significant = significant/length(fisher.z))
sig.zcor <- data.frame(method = sig.zcor$method,
											 label = paste("n = ", sig.zcor$significant, ", ", round(sig.zcor$p.significant,3)*100, "% significant", sep = ""))

meta.wcor %>% mutate(z.fixef = est.z.fixef/se.est.z.fixef,
										 z.ranef = est.z.ranef/se.est.z.ranef,
										 z.reg = est.z.reg/se.est.z.reg,
										 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	# filter(abs(fisher.z) < 10) %>% 
	ggplot(aes(x = abs(fisher.z))) + geom_histogram() + theme_bw() + facet_wrap(~method, ncol = 2) + 
	geom_vline(xintercept = 1.96, color = "red") + 
	geom_text(data = sig.zcor, aes(x = 30, y = 500, label = label), position = "dodge")

#Comparison of m.a. with copas correction method applied:
sig.zcor <- meta.wcor %>% filter(!is.na(est.z.copas)) %>% 
	mutate(z.fixef = est.z.fixef/se.est.z.fixef,
				 z.ranef = est.z.ranef/se.est.z.ranef,
				 z.reg = est.z.reg/se.est.z.reg,
				 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	group_by(method) %>% 
	summarise(significant = length(which(abs(fisher.z) > 1.96)),
						p.significant = significant/length(fisher.z))
sig.zcor <- data.frame(method = sig.zcor$method,
											 label = paste("n = ", sig.zcor$significant, ", ", round(sig.zcor$p.significant,3)*100, "% significant", sep = ""))

meta.wcor %>% filter(!is.na(est.z.copas)) %>% 
	mutate(z.fixef = est.z.fixef/se.est.z.fixef,
				 z.ranef = est.z.ranef/se.est.z.ranef,
				 z.reg = est.z.reg/se.est.z.reg,
				 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	# filter(abs(fisher.z) < 10) %>% 
	ggplot(aes(x = abs(fisher.z))) + geom_histogram() + theme_bw() + facet_wrap(~method, ncol = 2) + 
	geom_vline(xintercept = 1.96, color = "red") + 
	geom_text(data = sig.zcor, aes(x = 10, y = 125, label = label), position = "dodge")

