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

load(file.path(PATH_RESULTS, file = "data.processed.RData"))

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

data.ext2 <- pb.process2(data)

metac.bin <- meta.bin %>% filter(n.sig.type.bin > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac.cont <- meta.cont %>% filter(n.sig.type.cont > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac <- meta %>% filter(n.sig.type > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)

load(file.path(PATH_RESULTS, "meta.complete.RData"))

#Significant proportion
test.bin <- metac.bin %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																									schwarzer.test = mean(schwarzer.test),
																									rucker.test = mean(rucker.test),
																									harbord.test = mean(harbord.test),
																									peter.test = mean(peter.test))
test.bin <- test.bin %>% gather(key = "test.type", value = "mean")

p.bin <- metac.bin %>% ungroup() %>% 
	select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>%  
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Binary Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.bin$test.type, y = 500, 
					 label = paste(round(test.bin$mean, 2)*100, "% rejected"), 
					 color = "white")

test.cont <- metac.cont %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																										begg.test = mean(begg.test),
																										thomson.test = mean(thomson.test))

test.cont <- test.cont %>% gather(key = "test.type", value = "mean")

p.cont <- metac.cont %>% ungroup() %>% 
	select(egger.test, thomson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Continuous Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.cont$test.type, y = 150, 
					 label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
					 color = "white")

dat_text <- data.frame(
	label = paste(c(sum(metac.cont$begg.test), sum(metac.cont$egger.test), sum(metac.cont$thomson.test)), "< 0.1,",
								round(c(mean(metac.cont$begg.test), mean(metac.cont$egger.test), mean(metac.cont$thomson.test)),2)*100, "%"),
	test.type   = c("pval.begg.cont", "pval.egger.cont", "pval.thomson.cont"))

labels <- c(pval.begg.cont = "Begg Mazumdar", pval.egger.cont = "Egger", pval.thomson.cont = "Thompson Sharp")
p.dist.cont <- metac.cont %>% ungroup() %>% 
	select(pval.egger.cont, pval.thomson.cont, pval.begg.cont) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(bins = 20) + theme_bw() + 
	facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 25, label = label),  color = "black") + 
	theme(strip.text.x = element_text(size=15))


dat_text <- data.frame(
	label = paste(c(sum(metac.bin$harbord.test), sum(metac.bin$peter.test), 
									sum(metac.bin$rucker.test), sum(metac.bin$schwarzer.test)), "< 0.1,",
								round(c(mean(metac.bin$harbord.test), mean(metac.bin$peter.test), 
												mean(metac.bin$rucker.test), mean(metac.bin$schwarzer.test)),2)*100, "%"),
	test.type   = c("pval.harbord", "pval.peter", "pval.rucker", "pval.schwarzer"))

labels <- c(pval.harbord = "Harbord", pval.peter = "Peter", pval.rucker = "Rucker", pval.schwarzer = "Schwarzer")

p.dist.bin <- metac.bin %>% ungroup() %>% 
	select(pval.harbord, pval.peter, pval.rucker, pval.schwarzer) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(bins = 20) + theme_bw() + facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 100, label = label), color = "black") + 
	theme(strip.text.x = element_text(size=15))


#Test agreement
agree.bin <- metac.bin %>% mutate(n.sig = peter.test + rucker.test + harbord.test + schwarzer.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + xlab("Number of significant tests") + ylab("count") + ggtitle("Binary Outcomes")

agree.cont <- metac.cont %>% mutate(n.sig = egger.test + thomson.test + begg.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + xlab("Number of significant tests") + ylab("count") + ggtitle("Continuous Outcomes")


####################################################################################
####################################################################################
####################################################################################
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



####################################################################################
####################################################################################
####################################################################################
#####Only one systematic review restriction:
####################################################################################
####################################################################################
####################################################################################

metac.bin <- metac.bin %>% distinct(file.nr, .keep_all = T)
metac.cont <- metac.cont %>% distinct(file.nr, .keep_all = T)
metac <- metac %>% distinct(file.nr, .keep_all = T)

#Significant proportion
test.bin <- metac.bin %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																									schwarzer.test = mean(schwarzer.test),
																									rucker.test = mean(rucker.test),
																									harbord.test = mean(harbord.test),
																									peter.test = mean(peter.test))
test.bin <- test.bin %>% gather(key = "test.type", value = "mean")

p.bin <- metac.bin %>% ungroup() %>% 
	select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>%  
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Binary Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.bin$test.type, y = 500, 
					 label = paste(round(test.bin$mean, 2)*100, "% rejected"), 
					 color = "white")

test.cont <- metac.cont %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																										begg.test = mean(begg.test),
																										thomson.test = mean(thomson.test))

test.cont <- test.cont %>% gather(key = "test.type", value = "mean")

p.cont <- metac.cont %>% ungroup() %>% 
	select(egger.test, thomson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Continuous Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.cont$test.type, y = 150, 
					 label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
					 color = "white")

dat_text <- data.frame(
	label = paste(c(sum(metac.cont$begg.test), sum(metac.cont$egger.test), sum(metac.cont$thomson.test)), "< 0.1,",
								round(c(mean(metac.cont$begg.test), mean(metac.cont$egger.test), mean(metac.cont$thomson.test)),2)*100, "%"),
	test.type   = c("pval.begg", "pval.egger", "pval.thomson"))

labels <- c(pval.begg = "Begg Mazumdar", pval.egger = "Egger", pval.thomson = "Thompson Sharp")
p.dist.cont <- metac.cont %>% ungroup() %>% 
	select(pval.egger, pval.thomson, pval.begg) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(bins = 20) + theme_bw() + 
	facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 25, label = label),  color = "black") + 
	theme(strip.text.x = element_text(size=15))


dat_text <- data.frame(
	label = paste(c(sum(metac.bin$harbord.test), sum(metac.bin$peter.test), 
									sum(metac.bin$rucker.test), sum(metac.bin$schwarzer.test)), "< 0.1,",
								round(c(mean(metac.bin$harbord.test), mean(metac.bin$peter.test), 
												mean(metac.bin$rucker.test), mean(metac.bin$schwarzer.test)),2)*100, "%"),
	test.type   = c("pval.harbord", "pval.peter", "pval.rucker", "pval.schwarzer"))

labels <- c(pval.harbord = "Harbord", pval.peter = "Peter", pval.rucker = "Rucker", pval.schwarzer = "Schwarzer")

p.dist.bin <- metac.bin %>% ungroup() %>% 
	select(pval.harbord, pval.peter, pval.rucker, pval.schwarzer) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(bins = 20) + theme_bw() + facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 100, label = label), color = "black") + 
	theme(strip.text.x = element_text(size=15))


#Test agreement
agree.bin <- metac.bin %>% mutate(n.sig = peter.test + rucker.test + egger.test + harbord.test + schwarzer.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + xlab("Number of significant tests") + ylab("count") + ggtitle("Binary Outcomes")

agree.cont <- metac.cont %>% mutate(n.sig = egger.test + thomson.test + begg.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + xlab("Number of significant tests") + ylab("count") + ggtitle("Continuous Outcomes")



meta.wcor <- meta.wcor %>% distinct(file.nr, .keep_all = T)

####################################################################################
####################################################################################
####################################################################################
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
	geom_text(data = sig.zcor, aes(x = 30, y = 150, label = label), position = "dodge")

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
	geom_text(data = sig.zcor, aes(x = 10, y = 35, label = label), position = "dodge")

