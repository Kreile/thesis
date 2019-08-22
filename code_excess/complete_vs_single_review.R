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



#Explanation of coding: 0 = (unsig, unsig), 1 = (unsig, sig), 2 = (sig, sig), 3 = (sig, unsig)
sig.level <- 0.05
meta.f <- meta.f %>% rowwise() %>% 
	mutate(
		pval.copas = 2*(1-pnorm(abs(est.copas/se.est.copas))),
		
		sig.fixef = ifelse(pval.fixef > sig.level, 0, 1),
		sig.ranef = ifelse(pval.ranef > sig.level, 0, 1),
		sig.reg = ifelse(pval.reg > sig.level, 0, 1),
		sig.copas = ifelse(pval.copas > sig.level, 0, 1),
		
		sig.change.ranef.reg = sig.change(sig.before = sig.ranef, sig.after =  sig.reg), #Function to register if significance of estimate has changed after correction
		sig.change.ranef.copas = sig.change(sig.before = sig.ranef, sig.after =  sig.copas),
		sig.change.fixef.reg = sig.change(sig.before = sig.fixef, sig.after =  sig.reg),
		sig.change.fixef.copas = sig.change(sig.before = sig.fixef, sig.after =  sig.copas))


sig.level <- 0.1
meta.f <- meta.f %>% rowwise() %>% 
	mutate(egger.test = ifelse(pval.egger < sig.level, 1, 0),
				 thompson.test = ifelse(pval.thompson < sig.level, 1, 0),
				 begg.test = ifelse(pval.begg < sig.level, 1, 0),
				 
				 schwarzer.test = ifelse(pval.schwarzer < sig.level, 1, 0),
				 rucker.test = ifelse(pval.rucker < sig.level, 1, 0),
				 harbord.test = ifelse(pval.harbord < sig.level, 1, 0),
				 peter.test = ifelse(pval.peter < sig.level, 1, 0))

meta.bin <- meta.f %>% filter(outcome.type == "bin")
meta.f <- meta.f %>% filter(outcome.type == "cont")
meta.surv <- meta.f %>% filter(outcome.type == "surv")

#Comparison of "original effect sizes":
effect.diff <- meta.f %>% mutate(
	est.fixef = ifelse(outcome.type == "bin", exp(est.fixef), est.fixef),
	est.ranef = ifelse(outcome.type == "bin", exp(est.ranef), est.ranef),
	est.copas = ifelse(outcome.type == "bin", exp(est.copas), est.copas),
	est.reg = ifelse(outcome.type == "bin", exp(est.reg), est.reg)
)

effect.diff <- effect.diff %>% mutate(
	est.fixef = ifelse(outcome.type == "bin", est.fixef - 1, est.fixef),
	est.ranef = ifelse(outcome.type == "bin", est.ranef - 1, est.ranef),
	est.copas = ifelse(outcome.type == "bin", est.copas - 1, est.copas),
	est.reg = ifelse(outcome.type == "bin", est.reg - 1, est.reg),
	fixef.ranef = ifelse(abs(est.ranef) > abs(est.fixef), "Reduction", "Amplification"),
	fixef.copas = ifelse(abs(est.fixef) > abs(est.copas), "Reduction", "Amplification"),
	fixef.reg = ifelse(abs(est.fixef) > abs(est.reg), "Reduction", "Amplification"),
	fixef.copas.s = ifelse(sign(est.fixef) == sign(est.copas), "Unchanged", "Reversed"),
	fixef.reg.s = ifelse(sign(est.fixef) == sign(est.reg), "Unchanged", "Reversed")) 








#Significant proportion
test.bin <- meta.f %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																									schwarzer.test = mean(schwarzer.test),
																									rucker.test = mean(rucker.test),
																									harbord.test = mean(harbord.test),
																									peter.test = mean(peter.test))
test.bin <- test.bin %>% gather(key = "test.type", value = "mean")

p.bin <- meta.f %>% ungroup() %>% 
	select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>%  
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Binary Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.bin$test.type, y = 500, 
					 label = paste(round(test.bin$mean, 2)*100, "% rejected"), 
					 color = "white")

test.cont <- meta.f %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																										begg.test = mean(begg.test),
																										thompson.test = mean(thompson.test))

test.cont <- test.cont %>% gather(key = "test.type", value = "mean")

p.cont <- meta.f %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Continuous Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.cont$test.type, y = 150, 
					 label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
					 color = "white")

dat_text <- data.frame(
	label = paste(c(sum(meta.f$begg.test), sum(meta.f$egger.test), sum(meta.f$thompson.test)), "< 0.1,",
								round(c(mean(meta.f$begg.test), mean(meta.f$egger.test), mean(meta.f$thompson.test)),2)*100, "%"),
	test.type   = c("pval.begg.cont", "pval.egger", "pval.thompson.cont"))

labels <- c(pval.begg.cont = "Begg Mazumdar", pval.egger = "Egger", pval.thompson.cont = "Thompson Sharp")
p.dist.cont <- meta.f %>% ungroup() %>% 
	select(pval.egger, pval.thompson.cont, pval.begg.cont) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(bins = 20) + theme_bw() + 
	facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 25, label = label),  color = "black") + 
	theme(strip.text.x = element_text(size=15))


dat_text <- data.frame(
	label = paste(c(sum(meta.f$harbord.test), sum(meta.f$peter.test), 
									sum(meta.f$rucker.test), sum(meta.f$schwarzer.test)), "< 0.1,",
								round(c(mean(meta.f$harbord.test), mean(meta.f$peter.test), 
												mean(meta.f$rucker.test), mean(meta.f$schwarzer.test)),2)*100, "%"),
	test.type   = c("pval.harbord", "pval.peter", "pval.rucker", "pval.schwarzer"))

labels <- c(pval.harbord = "Harbord", pval.peter = "Peter", pval.rucker = "Rucker", pval.schwarzer = "Schwarzer")

p.dist.bin <- meta.f %>% ungroup() %>% 
	select(pval.harbord, pval.peter, pval.rucker, pval.schwarzer) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(bins = 20) + theme_bw() + facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 100, label = label), color = "black") + 
	theme(strip.text.x = element_text(size=15))


#Test agreement
agree.bin <- meta.f %>% mutate(n.sig = peter.test + rucker.test + harbord.test + schwarzer.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + xlab("Number of significant tests") + ylab("count") + ggtitle("Binary Outcomes")

agree.cont <- meta.f %>% mutate(n.sig = egger.test + thompson.test + begg.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + xlab("Number of significant tests") + ylab("count") + ggtitle("Continuous Outcomes")


####################################################################################
####################################################################################
####################################################################################
####################################################################################
#Comparison of all adjusted z-scores:
sig.zcor <- meta.f %>% mutate(z.fixef = est.z.fixef/se.est.z.fixef,
																 z.ranef = est.z.ranef/se.est.z.ranef,
																 z.reg = est.z.reg/se.est.z.reg,
																 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	group_by(method) %>% 
	summarise(significant = length(which(abs(fisher.z) > 1.96)),
						p.significant = significant/length(fisher.z))
sig.zcor <- data.frame(method = sig.zcor$method,
											 label = paste("n = ", sig.zcor$significant, ", ", round(sig.zcor$p.significant,3)*100, "% significant", sep = ""))

meta.f %>% mutate(z.fixef = est.z.fixef/se.est.z.fixef,
										 z.ranef = est.z.ranef/se.est.z.ranef,
										 z.reg = est.z.reg/se.est.z.reg,
										 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	# filter(abs(fisher.z) < 10) %>% 
	ggplot(aes(x = abs(fisher.z))) + geom_histogram() + theme_bw() + facet_wrap(~method, ncol = 2) + 
	geom_vline(xintercept = 1.96, color = "red") + 
	geom_text(data = sig.zcor, aes(x = 30, y = 500, label = label), position = "dodge")

#Comparison of m.a. with copas correction method applied:
sig.zcor <- meta.f %>% filter(!is.na(est.z.copas)) %>% 
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

meta.f %>% filter(!is.na(est.z.copas)) %>% 
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

meta.f <- meta.f %>% distinct(file.nr, .keep_all = T)
meta.f <- meta.f %>% distinct(file.nr, .keep_all = T)
meta.f <- meta.f %>% distinct(file.nr, .keep_all = T)

#Significant proportion
test.bin <- meta.f %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																									schwarzer.test = mean(schwarzer.test),
																									rucker.test = mean(rucker.test),
																									harbord.test = mean(harbord.test),
																									peter.test = mean(peter.test))
test.bin <- test.bin %>% gather(key = "test.type", value = "mean")

p.bin <- meta.f %>% ungroup() %>% 
	select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>%  
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Binary Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.bin$test.type, y = 500, 
					 label = paste(round(test.bin$mean, 2)*100, "% rejected"), 
					 color = "white")

test.cont <- meta.f %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																										begg.test = mean(begg.test),
																										thompson.test = mean(thompson.test))

test.cont <- test.cont %>% gather(key = "test.type", value = "mean")

p.cont <- meta.f %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Continuous Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.cont$test.type, y = 150, 
					 label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
					 color = "white")

dat_text <- data.frame(
	label = paste(c(sum(meta.f$begg.test), sum(meta.f$egger.test), sum(meta.f$thompson.test)), "< 0.1,",
								round(c(mean(meta.f$begg.test), mean(meta.f$egger.test), mean(meta.f$thompson.test)),2)*100, "%"),
	test.type   = c("pval.begg", "pval.egger", "pval.thompson"))

labels <- c(pval.begg = "Begg Mazumdar", pval.egger = "Egger", pval.thompson = "Thompson Sharp")
p.dist.cont <- meta.f %>% ungroup() %>% 
	select(pval.egger, pval.thompson, pval.begg) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(bins = 20) + theme_bw() + 
	facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 25, label = label),  color = "black") + 
	theme(strip.text.x = element_text(size=15))


dat_text <- data.frame(
	label = paste(c(sum(meta.f$harbord.test), sum(meta.f$peter.test), 
									sum(meta.f$rucker.test), sum(meta.f$schwarzer.test)), "< 0.1,",
								round(c(mean(meta.f$harbord.test), mean(meta.f$peter.test), 
												mean(meta.f$rucker.test), mean(meta.f$schwarzer.test)),2)*100, "%"),
	test.type   = c("pval.harbord", "pval.peter", "pval.rucker", "pval.schwarzer"))

labels <- c(pval.harbord = "Harbord", pval.peter = "Peter", pval.rucker = "Rucker", pval.schwarzer = "Schwarzer")

p.dist.bin <- meta.f %>% ungroup() %>% 
	select(pval.harbord, pval.peter, pval.rucker, pval.schwarzer) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(bins = 20) + theme_bw() + facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 100, label = label), color = "black") + 
	theme(strip.text.x = element_text(size=15))


#Test agreement
agree.bin <- meta.f %>% mutate(n.sig = peter.test + rucker.test + egger.test + harbord.test + schwarzer.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + xlab("Number of significant tests") + ylab("count") + ggtitle("Binary Outcomes")

agree.cont <- meta.f %>% mutate(n.sig = egger.test + thompson.test + begg.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + xlab("Number of significant tests") + ylab("count") + ggtitle("Continuous Outcomes")



meta.f <- meta.f %>% distinct(file.nr, .keep_all = T)

####################################################################################
####################################################################################
####################################################################################
####################################################################################
#Comparison of all adjusted z-scores:
sig.zcor <- meta.f %>% mutate(z.fixef = est.z.fixef/se.est.z.fixef,
																 z.ranef = est.z.ranef/se.est.z.ranef,
																 z.reg = est.z.reg/se.est.z.reg,
																 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	group_by(method) %>% 
	summarise(significant = length(which(abs(fisher.z) > 1.96)),
						p.significant = significant/length(fisher.z))
sig.zcor <- data.frame(method = sig.zcor$method,
											 label = paste("n = ", sig.zcor$significant, ", ", round(sig.zcor$p.significant,3)*100, "% significant", sep = ""))

meta.f %>% mutate(z.fixef = est.z.fixef/se.est.z.fixef,
										 z.ranef = est.z.ranef/se.est.z.ranef,
										 z.reg = est.z.reg/se.est.z.reg,
										 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	# filter(abs(fisher.z) < 10) %>% 
	ggplot(aes(x = abs(fisher.z))) + geom_histogram() + theme_bw() + facet_wrap(~method, ncol = 2) + 
	geom_vline(xintercept = 1.96, color = "red") + 
	geom_text(data = sig.zcor, aes(x = 30, y = 150, label = label), position = "dodge")

#Comparison of m.a. with copas correction method applied:
sig.zcor <- meta.f %>% filter(!is.na(est.z.copas)) %>% 
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

meta.f %>% filter(!is.na(est.z.copas)) %>% 
	mutate(z.fixef = est.z.fixef/se.est.z.fixef,
				 z.ranef = est.z.ranef/se.est.z.ranef,
				 z.reg = est.z.reg/se.est.z.reg,
				 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	# filter(abs(fisher.z) < 10) %>% 
	ggplot(aes(x = abs(fisher.z))) + geom_histogram() + theme_bw() + facet_wrap(~method, ncol = 2) + 
	geom_vline(xintercept = 1.96, color = "red") + 
	geom_text(data = sig.zcor, aes(x = 10, y = 35, label = label), position = "dodge")



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



test.sig.bin <- meta.bin %>% filter(sig.fixef == 1) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																																								schwarzer.test = mean(schwarzer.test),
																																								rucker.test = mean(rucker.test),
																																								harbord.test = mean(harbord.test),
																																								peter.test = mean(peter.test))
test.sig.bin <- test.sig.bin %>% gather(key = "test.type", value = "mean")

p1 <- meta.bin %>% ungroup() %>% filter(sig.fixef == 1) %>% 
	select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + ggtitle("Significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
	annotate("text", x = test.sig.bin$test.type, y = 1750, 
					 label = paste(round(test.sig.bin$mean, 2)*100, "% rejected"), 
					 color = "white")

test.nonsig.bin <- meta.bin %>% filter(sig.fixef == 0) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																																									 schwarzer.test = mean(schwarzer.test),
																																									 rucker.test = mean(rucker.test),
																																									 harbord.test = mean(harbord.test),
																																									 peter.test = mean(peter.test))
test.nonsig.bin <- test.nonsig.bin %>% gather(key = "test.type", value = "mean")

p2 <- meta.bin %>% 
	filter(sig.fixef == 0) %>% ungroup() %>% 
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
test.cont <- meta.f %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																									 begg.test = mean(begg.test),
																									 thompson.test = mean(thompson.test))

test.cont <- test.cont %>% gather(key = "test.type", value = "mean")

p.cont <- meta.f %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) +
	annotate("text", x = test.cont$test.type, y = 750, 
					 label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
					 color = "white")


test.sig.cont <- meta.f %>% filter(sig.fixef == 1) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																																									begg.test = mean(begg.test),
																																									thompson.test = mean(thompson.test))

test.sig.cont <- test.sig.cont %>% gather(key = "test.type", value = "mean")

p3 <- meta.f %>% 
	filter(sig.fixef == 1) %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + ggtitle("Significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
	annotate("text", x = test.sig.cont$test.type, y = 600, 
					 label = paste(round(test.sig.cont$mean, 2)*100, "% rejected"), 
					 color = "white")

test.nonsig.cont <- meta.f %>% filter(sig.fixef == 0) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																																										 begg.test = mean(begg.test),
																																										 thompson.test = mean(thompson.test))

test.nonsig.cont <- test.nonsig.cont %>% gather(key = "test.type", value = "mean")

p4 <- meta.f %>% 
	filter(sig.fixef == 0) %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
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
meta.f <- meta.f %>% distinct(file.nr, .keep_all = T)

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



test.sig.bin <- meta.bin %>% filter(sig.fixef == 1) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																																								schwarzer.test = mean(schwarzer.test),
																																								rucker.test = mean(rucker.test),
																																								harbord.test = mean(harbord.test),
																																								peter.test = mean(peter.test))
test.sig.bin <- test.sig.bin %>% gather(key = "test.type", value = "mean")

p1s <- meta.bin %>% ungroup() %>% filter(sig.fixef == 1) %>% 
	select(egger.test, schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + ggtitle("Significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
	annotate("text", x = test.sig.bin$test.type, y = 1750, 
					 label = paste(round(test.sig.bin$mean, 2)*100, "% rejected"), 
					 color = "white")

test.nonsig.bin <- meta.bin %>% filter(sig.fixef == 0) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																																									 schwarzer.test = mean(schwarzer.test),
																																									 rucker.test = mean(rucker.test),
																																									 harbord.test = mean(harbord.test),
																																									 peter.test = mean(peter.test))
test.nonsig.bin <- test.nonsig.bin %>% gather(key = "test.type", value = "mean")

p2s <- meta.bin %>% 
	filter(sig.fixef == 0) %>% ungroup() %>% 
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
test.cont <- meta.f %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																									 begg.test = mean(begg.test),
																									 thompson.test = mean(thompson.test))

test.cont <- test.cont %>% gather(key = "test.type", value = "mean")

p.conts <- meta.f %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) +
	annotate("text", x = test.cont$test.type, y = 150, 
					 label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
					 color = "white")


test.sig.cont <- meta.f %>% filter(sig.fixef == 1) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																																									begg.test = mean(begg.test),
																																									thompson.test = mean(thompson.test))

test.sig.cont <- test.sig.cont %>% gather(key = "test.type", value = "mean")

p3s <- meta.f %>% 
	filter(sig.fixef == 1) %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "rejected", "not rejected"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + ggtitle("Significant Pooled Effects")+ xlab(label = NULL) + theme(legend.position="none") +
	annotate("text", x = test.sig.cont$test.type, y = 100, 
					 label = paste(round(test.sig.cont$mean, 2)*100, "% rejected"), 
					 color = "white")

test.nonsig.cont <- meta.f %>% filter(sig.fixef == 0) %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																																										 begg.test = mean(begg.test),
																																										 thompson.test = mean(thompson.test))

test.nonsig.cont <- test.nonsig.cont %>% gather(key = "test.type", value = "mean")

p4s <- meta.f %>% 
	filter(sig.fixef == 0) %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
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

