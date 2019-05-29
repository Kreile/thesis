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





#Meta filtering: 
metac.bin <- meta.bin %>% filter(n.sig.type. > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac.cont <- meta.cont %>% filter(n.sig.type.cont > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
metac <- meta %>% filter(n.sig.type > 1) %>% filter((se.max^2)/(se.min^2) > 4) %>% filter(I2 < 0.5)
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



#Adjustment example plots regression:
reg.1 <- qplot(x = biased.rev$se, y = biased.rev$effect, xlim = c(0, 1.2)) + geom_point() + 
  stat_smooth(method="lm", fullrange=TRUE, se = F) + coord_flip() + theme_bw() +
  ylab("Std. Mean Difference") + xlab("Standard error")

reg.2 <- qplot(x = unbiased.rev$se, y = unbiased.rev$effect, xlim = c(0, 6.5)) + geom_point() + 
  stat_smooth(method="lm", fullrange=TRUE, se = F) + coord_flip() + theme_bw() +
  ylab("Std. Mean Difference") + xlab("Standard error")

################################################################################################
################################################################################################
################################################################################################
######## Histogramms of effect size changes of adjustments:
#Adjustment results: Effect difference histograms
effect.diff.tmp <- metac %>% mutate(
  est.fixef = ifelse(outcome.type == "bin", est.fixef - 1, est.fixef),
  est.ranef = ifelse(outcome.type == "bin", est.ranef - 1, est.ranef),
  est.trimfill.fixef = ifelse(outcome.type == "bin", est.trimfill.fixef - 1, est.trimfill.fixef),
  est.copas = ifelse(outcome.type == "bin", est.copas - 1, est.copas),
  est.reg.ranef = ifelse(outcome.type == "bin", est.reg.ranef - 1, est.reg.ranef),
  fixef.ranef = ifelse(abs(est.fixef) > abs(est.ranef), "Reduction", "Amplification"),
  fixef.trimfill = ifelse(abs(est.fixef) > abs(est.trimfill.fixef), "Reduction", "Amplification"),
  fixef.copas = ifelse(abs(est.fixef) > abs(est.copas), "Reduction", "Amplification"),
  fixef.reg = ifelse(abs(est.fixef) > abs(est.reg.ranef), "Reduction", "Amplification"), 
  outcome.type = factor(ifelse(outcome.type == "bin", "Risk Ratios", "Mean Differences")),
  outcome.measure = factor(ifelse(outcome.type == "bin", "Risk Ratio", 
  												 ifelse(outcome.measure.new == "Std. Mean Difference", "Std. Mean Difference",
  												 			 "Mean Difference"))))
  # fixef.trimfill.s = ifelse(sign(est.fixef) == sign(est.trimfill.fixef), "Unchanged", "Reversed"),
  # fixef.copas.s = ifelse(sign(est.fixef) == sign(est.copas), "Unchanged", "Reversed"),
  # fixef.reg.s = ifelse(sign(est.fixef) == sign(est.reg.ranef), "Unchanged", "Reversed")) 


effect.diff <- effect.diff.tmp %>% 
	mutate(change.fixef.trimfill = ifelse(fixef.trimfill == "Reduction", -abs(est.fixef - est.trimfill.fixef),
																				abs(est.fixef - est.trimfill.fixef)),
				 change.fixef.copas = ifelse(fixef.copas == "Reduction", -abs(est.fixef - est.copas),
				 														abs(est.fixef - est.copas)),
				 change.fixef.reg = ifelse(fixef.ranef == "Reduction", -abs(est.fixef - est.reg.ranef),
				 													abs(est.fixef - est.reg.ranef)))

#Trim-and-Fill
trimfill.miss <- effect.diff %>% group_by(outcome.type) %>% mutate(dif = abs(est.fixef - est.trimfill.fixef)) %>% 
	filter(dif > 2) %>% mutate(label = round(dif, 1)) %>% arrange(label) %>% select(label, meta.id, outcome.measure.new)
trimfill.label <- data.frame(label = paste("missing:", c(paste(trimfill.miss$label[trimfill.miss$outcome.type == "Risk Ratios"], collapse = ", "),
																												 paste(trimfill.miss$label[trimfill.miss$outcome.type == "Mean Differences"], collapse = ", "))),
														 outcome.type = c("Risk Ratios", "Mean Differences"), 
														 outcome.measure = unique(trimfill.miss$outcome.measure.new))
trimfill.label$label <- as.character(trimfill.label$label)
trimfill.label$label[2] <- sub("2.3,", "2.3,\n", trimfill.label$label[2])

effect.diff %>% ggplot(aes(x = change.fixef.trimfill, fill = outcome.measure)) + 
	facet_wrap(outcome.type~., scales = "free", labeller = labeller(outcome.type = NULL)) +
  geom_histogram() + xlim(-2, 2) + xlab("Effect change") + 
  theme_bw() + theme(legend.position = "bottom", strip.text = element_blank(), legend.title = element_blank())  +
	geom_text(data = trimfill.label, 
						mapping = aes(x = 0, y = c(460, 100), label = label, vjust = "center"), position = "dodge") #+ ylim(0, 250)

#Copas selection model:
copas.miss <- effect.diff %>% group_by(outcome.type) %>% mutate(dif = abs(est.fixef - est.copas)) %>% 
	filter(dif > 2) %>% mutate(label = round(dif, 1)) %>% arrange(label) %>% select(label, meta.id, outcome.measure.new)
copas.label <- data.frame(label = paste("missing:", c(paste(copas.miss$label[copas.miss$outcome.type == "Risk Ratios"], collapse = ", "),
																												 paste(copas.miss$label[copas.miss$outcome.type == "Mean Differences"], collapse = ", "))),
														 outcome.type = c("Risk Ratios", "Mean Differences"), 
														 outcome.measure = unique(copas.miss$outcome.measure.new))
copas.label$label <- as.character(copas.label$label)
copas.label$label[2] <- sub("2.8,", "2.8,\n", copas.label$label[2])

effect.diff %>% ggplot(aes(x = change.fixef.copas, fill = outcome.measure)) + 
	facet_wrap(outcome.type~., scales = "free", labeller = labeller(outcome.type = NULL)) +
	geom_histogram() + xlim(-2, 2) + xlab("Effect change") + 
	theme_bw() + theme(legend.position = "bottom", strip.text = element_blank(), legend.title = element_blank())  +
	geom_text(data = copas.label, 
						mapping = aes(x = 0, y = c(220, 47.4), label = label, vjust = "center"), position = "dodge") #+ ylim(0, 250)

#Regression model:
reg.miss <- effect.diff %>% group_by(outcome.type) %>% mutate(dif = abs(est.fixef - est.reg.ranef)) %>% 
	filter(dif > 2) %>% mutate(label = round(dif, 1)) %>% arrange(label) %>% select(label, meta.id, outcome.measure.new)
reg.label <- data.frame(label = paste("missing:", c(paste(reg.miss$label[reg.miss$outcome.type == "Risk Ratios"], collapse = ", "),
																											paste(reg.miss$label[reg.miss$outcome.type == "Mean Differences"], collapse = ", "))),
													outcome.type = c("Risk Ratios", "Mean Differences"), 
													outcome.measure = unique(reg.miss$outcome.measure.new)[1:2])
reg.label$label <- as.character(reg.label$label)
reg.label$label[2] <- sub("2.2,","2.2,\n ", reg.label$label[2])
reg.label$label[2] <- sub("2.6,","2.6,\n ", reg.label$label[2])
reg.label$label[2] <- sub(", 3.1,",", 3.1,\n ", reg.label$label[2])
reg.label$label[2] <- sub("4.4, 4.4,", "4.4, 4.4,\n", reg.label$label[2])
reg.label$label[2] <- sub(", 9.5,", ", 9.5,\n", reg.label$label[2])

reg.label$label[1] <- sub(", 3.3,", ", 3.3,\n", reg.label$label[1])
reg.label$label[1] <- sub(", 7.6,", ", 7.6,\n", reg.label$label[1])

effect.diff %>% ggplot(aes(x = change.fixef.reg, fill = outcome.measure)) + 
	facet_wrap(outcome.type~., scales = "free", labeller = labeller(outcome.type = NULL)) +
	geom_histogram() + xlim(-2, 2) + xlab("Effect change") + 
	theme_bw() + theme(legend.position = "bottom", strip.text = element_blank(), legend.title = element_blank())  +
	geom_text(data = reg.label, 
						mapping = aes(x = 0, y = c(330, 50), label = label, vjust = "center"), position = "dodge") #+ ylim(0, 250)

################################################################################################
################################################################################################
################################################################################################

#Scatterplots of metac-analysis and corrected log estimates
ranef.fixef.sc.effect <- metac %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.ranef)))) + geom_point(size = 0.8) + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
  xlab("Fixed effects log estimate") + ylab("Random effects log estimate") + ggtitle("Fixed effects to Random effects")


#Scatterplots of metac-analysis and corrected log(abs(abs estimates
trimfill.fixef.sc.effect <- metac %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.trimfill.fixef)))) + geom_point(size = 0.8) + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
  xlab("Fixed effects log estimate") + ylab("Trimfill log estimate") + ggtitle("Fixed effects to Trimfill")

trimfill.ranef.sc.effect <- metac %>% ggplot(aes(x = log(abs((est.ranef))), y = log(abs(est.trimfill.fixef)))) + geom_point(size = 0.8) + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
  xlab("Random effects log estimate") + ylab("Trimfill log estimate") + ggtitle("Random effects to Trimfill")


#Scatterplots of metac-analysis and corrected log estimates
copas.fixef.sc.effect <- metac %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.copas)))) + geom_point(size = 0.8) + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
  xlab("Fixed effects log estimate") + ylab("Copas log estimate") + ggtitle("Fixed effects to Copas") 

copas.ranef.sc.effect <- metac %>% ggplot(aes(x = log(abs((est.ranef))), y = log(abs(est.copas)))) + geom_point(size = 0.8) + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE, show.legend = T) + 
  xlab("Random effects log estimate") + ylab("Copas log estimate") + ggtitle("Random effects to Copas") 


#Scatterplots of metac-analysis and corrected log estimates
reg.fixef.sc.effect <- metac %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.reg.ranef)))) + geom_point(size = 0.8) + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
  xlab("Fixed effects log estimate") + ylab("Regression log estimate") + ggtitle("Fixed effects to Regression")

reg.ranef.sc.effect <- metac %>% ggplot(aes(x = log(abs((est.ranef))), y = log(abs(est.reg.ranef)))) + geom_point(size = 0.8) + 
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
  xlab("Random effects log estimate") + ylab("Regression log estimate") + ggtitle("Random effects to Regression")


#Scatterplots of test- statistics:
ranef.fixef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.ranef))) + geom_point(size = 0.8) +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
  ggtitle("Fixed effects and Random effects z statistics") + ylab("Random effects statistic") + xlab("Fixed effects statistic")

metac %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.trimfill.fixef))) + geom_point(size = 0.8) +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
  ggtitle("Fixed effects and Random effects z statistics") + ylab("Random effects statistic") + xlab("Fixed effects statistic")


#Scatterplots of test- statistics:
trimfill.fixef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.trimfill.fixef))) + geom_point(size = 0.8) +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
  ggtitle("Fixed effects and trimfill z statistics") + ylab("Trimfill  statistic") + xlab("Fixed effects  statistic")

trimfill.ranef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.ranef), y = abs(zval.trimfill.fixef))) + geom_point(size = 0.8) +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
  ggtitle("Random effects and trimfill z statistics") + ylab("Trimfill statistic") + xlab("Random effects  statistic")

#Scatterplots of test- statistics:
copas.fixef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.copas))) + geom_point(size = 0.8) +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + xlim(0, 25) +
  ggtitle("Fixed effects and Copas z statistics") + ylab("Copas  statistic") + xlab("Fixed effects  statistic")

copas.ranef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.ranef), y = abs(zval.copas))) + geom_point(size = 0.8) +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + xlim(0, 20) +
  ggtitle("Random effects and Copas z statistics") + ylab("Copas  statistic") + xlab("Random effects  statistic")

#Scatterplots of test- statistics:
reg.fixef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.reg.ranef))) + geom_point(size = 0.8) +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
  ggtitle("Fixed effects and Regression z statistics") + ylab("Regression  statistic") + xlab("Fixed effects  statistic")

reg.ranef.sc.zval <- metac %>% ggplot(aes(x = abs(zval.ranef), y = abs(zval.reg.ranef))) + geom_point(size = 0.8) +
  geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
  ggtitle("Random effects and Regression z statistics") + ylab("Regression  statistic") + xlab("Random effects  statistic")

#Missing study proportions:
p.missing.copas <- metac %>% ggplot(aes(x = missing.copas/n)) + geom_histogram(bins = 20) + 
  xlim(0, 0.5) + ylim(0, 200) + xlab("Proportion of missing studies") + ggtitle("Copas selection model") + theme_bw()
p.missing.trim <- metac %>% ggplot(aes(x = missing.trim)) + geom_histogram(bins = 20) + 
  xlim(0, 0.5) + ylim(0, 200) + xlab("Proportion of missing studies") + ggtitle("Trim-and-fill") + theme_bw()