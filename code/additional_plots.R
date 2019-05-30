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
meta.f <- meta.f %>% 
	mutate(
		pval.copas = 2*(1-pnorm(abs(est.copas/se.est.copas))),
		
		sig.fixef = ifelse(pval.fixef > sig.level, 0, 1),
		sig.ranef = ifelse(pval.ranef > sig.level, 0, 1),
		sig.reg.ranef = ifelse(pval.reg > sig.level, 0, 1),
		sig.copas = ifelse(pval.copas > sig.level, 0, 1),
		
		sig.change.ranef.reg = sig.change(sig.before = sig.ranef, sig.after =  sig.reg), #Function to register if significance of estimate has changed after correction
		sig.change.ranef.copas = sig.change(sig.before = sig.ranef, sig.after =  sig.copas),
		sig.change.fixef.reg = sig.change(sig.before = sig.fixef, sig.after =  sig.reg),
		sig.change.fixef.copas = sig.change(sig.before = sig.fixef, sig.after =  sig.copas))


sig.level <- 0.1
meta.f <- meta.f %>% 
	mutate(egger.test = ifelse(pval.peter < sig.level, 1, 0),
				 thompson.test = ifelse(pval.thompson < sig.level, 1, 0),
				 begg.test = ifelse(pval.begg < sig.level, 1, 0),
				 
				 schwarzer.test = ifelse(pval.schwarzer < sig.level, 1, 0),
				 rucker.test = ifelse(pval.rucker < sig.level, 1, 0),
				 harbord.test = ifelse(pval.harbord < sig.level, 1, 0),
				 peter.test = ifelse(pval.peter < sig.level, 1, 0))

meta.bin <- meta.f %>% filter(outcome.type == "bin")
meta.cont <- meta.f %>% filter(outcome.type == "cont")
meta.surv <- meta.f %>% filter(outcome.type == "surv")

#Comparison of "original effect sizes":
effect.diff <- meta.f %>% mutate(
	est.fixef = ifelse(outcome.type == "bin", exp(est.fixef), est.fixef),
	est.ranef = ifelse(outcome.type == "bin", exp(est.ranef), est.ranef),
	est.copas = ifelse(outcome.type == "bin", exp(est.copas), est.copas),
	est.reg = ifelse(outcome.type == "bin", exp(est.reg.ranef), est.reg.ranef)
)

effect.diff <- effect.diff %>% mutate(
	est.fixef = ifelse(outcome.type == "bin", est.fixef - 1, est.fixef),
	est.ranef = ifelse(outcome.type == "bin", est.ranef - 1, est.ranef),
	est.copas = ifelse(outcome.type == "bin", est.copas - 1, est.copas),
	est.reg = ifelse(outcome.type == "bin", est.reg.ranef - 1, est.reg.ranef),
	fixef.ranef = ifelse(abs(est.ranef) > abs(est.fixef), "Reduction", "Amplification"),
	fixef.copas = ifelse(abs(est.fixef) > abs(est.copas), "Reduction", "Amplification"),
	fixef.reg = ifelse(abs(est.fixef) > abs(est.reg), "Reduction", "Amplification"),
	fixef.copas.s = ifelse(sign(est.fixef) == sign(est.copas), "Unchanged", "Reversed"),
	fixef.reg.s = ifelse(sign(est.fixef) == sign(est.reg), "Unchanged", "Reversed")) 





####################################################################################################
#Plots that have been created before and are potentially interesting also with the new data:
####################################################################################################

#Significance of pbbias based on variance of total ss:
meta %>% filter(var.samplesize < 60000) %>% ggplot(aes(stat(count), x = var.samplesize, fill = factor(thomson.test))) + 
	geom_density(position = "fill")
#--------------------------------------------------------------------------------------------------------------------#

#Mean and Total study sample size
meta %>% filter(mean.samplesize < 1500) %>%  ggplot(aes(x = mean.samplesize)) + geom_histogram()
meta %>% filter(total.samplesize < 50000) %>% ggplot(aes(x = total.samplesize)) + geom_histogram()
meta %>% filter(total.samplesize < 1000) %>% ggplot(aes(x = total.samplesize, fill = factor(thomson.test))) + 
	geom_histogram()
meta %>% filter(mean.samplesize < 100) %>% ggplot(aes(x = mean.samplesize, fill = factor(thomson.test))) + 
	geom_histogram()
#--------------------------------------------------------------------------------------------------------------------#

#Meta analysis sample size
meta %>% filter(n < 100) %>% ggplot(aes(x = n, fill = factor(thompson.test))) + 
	geom_histogram()
meta %>% filter(n < 40) %>% ggplot(aes(stat(count), x = n, fill = factor(thompson.test))) + 
	geom_density(position = "fill")

meta %>% filter(n < 100) %>% ggplot(aes(x = n, fill = factor(sig.fixef))) + 
	geom_histogram()
meta %>% filter(n < 40) %>% ggplot(aes(stat(count), x = n, fill = factor(sig.fixef))) + 
	geom_density(position = "fill")
#--------------------------------------------------------------------------------------------------------------------#

#Dependence of significance of publication bias tests on mean and total samplesize
meta %>% filter(mean.samplesize < 750) %>% ggplot(aes(stat(count), x = mean.samplesize, fill = factor(thompson.test))) + 
	geom_density(position = "fill")

meta %>% filter(total.samplesize < 10000) %>% ggplot(aes(stat(count), x = total.samplesize, fill = factor(thompson.test))) + 
	geom_density(position = "fill")
#--------------------------------------------------------------------------------------------------------------------#

# #Change in significance
# meta %>% filter(!is.na(sig.change)) %>% ggplot(aes(x = sig.change, stat = "count")) + geom_histogram(stat = "count")
#--------------------------------------------------------------------------------------------------------------------#

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
	ggplot(aes(x = first.publication.year, fill = factor(thompson.test))) + geom_density(na.rm = T, position = "fill")

meta.cont %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
	ggplot(aes(x = mean.publication.year, fill = factor(egger.test))) + geom_density(na.rm = T, position = "fill")


meta %>% filter(mean.publication.year < 2019 & mean.publication.year > 1960) %>% 
	ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef))) + geom_histogram() + 
	labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))

meta %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
	ggplot(aes(x = mean.publication.year, fill = factor(sig.fixef))) + geom_density(na.rm = T, position = "identity") + 
	labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))

meta %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
	ggplot(aes(x = mean.publication.year, fill = factor(thompson.test), stat(count))) + geom_density(na.rm = T, position = "fill") + 
	labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))

cdplot(x = meta$mean.publication.year[!is.na(meta$mean.publication.year)], y = factor(meta$thompson.test[!is.na(meta$mean.publication.year)]))

meta %>% filter(mean.publication.year < 2013 & mean.publication.year > 1990) %>% 
	ggplot(aes(x = mean.publication.year, fill = factor(egger.test))) + geom_density(na.rm = T, position = "fill")
#--------------------------------------------------------------------------------------------------------------------#

#Sample size median over the years
data %>% distinct(study.name, .keep_all = T) %>%  filter(study.year < 2018 & study.year > 1980) %>% group_by(study.year) %>% 
	summarize(samplesize = median(total1 + total2, na.rm = T)) %>% 
	ggplot(aes(y = samplesize, x = study.year)) + geom_line()
#--------------------------------------------------------------------------------------------------------------------#

require(GGally)
topairsplot.pb.stats.bin <- metac %>% select(stat.egger, stat.harbord, stat.peter, 
																						 stat.schwarzer, stat.rucker, harbord.test)
#topairsplot.pb.stats.bin <- abs(topairsplot.pb.stats.bin)
ggpairs(topairsplot.pb.stats.bin, columns = 1:5, aes(colour = factor(harbord.test), alpha = 0.5))
ggpairs(topairsplot.pb.stats.bin, columns = 1:5, aes(alpha = 0.5))

topairsplot.pb.stats.cont <- metac %>% select(stat.egger, stat.thomson, stat.begg,
																							egger.test)
#topairsplot.pb.stats.bin <- abs(topairsplot.pb.stats.bin)
ggpairs(topairsplot.pb.stats.cont, columns = 1:3, aes(colour = factor(egger.test)))
ggpairs(topairsplot.pb.stats.cont, columns = 1:3, aes(alpha = 0.7))
#--------------------------------------------------------------------------------------------------------------------#

#$Complete dataset pairsplot:
topairsplot <- metac %>% select(comparison.nr, pval.Q, pval.thomson, n.sig.type,	missing.copas,
																G.squared, n, var.samplesize,mean.publication.year
)
topairsplot <- topairsplot %>% filter(var.samplesize < 100000 & comparison.nr < 20)

ggpairs(topairsplot, aes(alpha = 0.1))
#--------------------------------------------------------------------------------------------------------------------#

#Comparing publication bias tests:
meta %>% ggplot(aes(y = pval.egger, x = pval.thomson)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta %>% ggplot(aes(x = pval.thomson, y = pval.begg)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta %>% filter(n < 100) %>% ggplot(aes(y = pval.thomson, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess", se = F) 
meta %>% filter(n < 100) %>% ggplot(aes(y = pval.egger, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess", se = F) 
meta %>% ggplot(aes(y = stat.egger, x = stat.thomson)) + geom_point(alpha = 0.3) + theme_bw() #+ geom_smooth(method = "lm")  
meta %>% ggplot(aes(x = pval.thomson, y = missing.trim)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 

meta %>% ggplot(aes(x = pval.harbord, y = pval.peter)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta %>% filter(n < 150)%>% ggplot(aes(y = pval.peter, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta %>% filter(n < 150)%>% ggplot(aes(y = pval.harbord, x = n)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta %>% ggplot(aes(y = pval.egger, x = pval.peter)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta %>% ggplot(aes(y = pval.egger, x = pval.harbord)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta %>% ggplot(aes(x = pval.schwarzer, y = pval.peter)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta %>% ggplot(aes(x = pval.schwarzer, y = pval.harbord)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta %>% ggplot(aes(y = pval.rucker, y = pval.egger)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta %>% ggplot(aes(y = pval.rucker, y = pval.peter)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 
meta %>% ggplot(aes(y = pval.rucker, y = pval.harbord)) + geom_point(alpha = 0.5) + theme_bw() + geom_smooth(method = "loess") 

meta %>% ggplot(aes(x = stat.peter, y = stat.schwarzer)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta %>% ggplot(aes(x = stat.peter, y = stat.egger)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta %>% ggplot(aes(x = stat.peter, y = stat.rucker)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta %>% ggplot(aes(x = stat.peter, y = stat.harbord)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess")

meta %>% ggplot(aes(y = stat.egger, x = stat.schwarzer)) + geom_point(alpha = 0.5) + theme_bw()# + geom_smooth(method = "loess") 
meta %>% ggplot(aes(y = stat.rucker, x = stat.schwarzer)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta %>% ggplot(aes(y = stat.harbord, x = stat.schwarzer)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 

meta %>% ggplot(aes(x = stat.rucker, y = stat.harbord)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 
meta %>% ggplot(aes(x = stat.rucker, y = stat.egger)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess") 

meta %>% ggplot(aes(x = stat.harbord, y = stat.egger)) + geom_point(alpha = 0.5) + theme_bw() #+ geom_smooth(method = "loess")
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#


effect.diff <- meta.f %>% mutate(
	est.fixef = ifelse(outcome.type == "bin", exp(est.fixef), est.fixef),
	est.ranef = ifelse(outcome.type == "bin", exp(est.ranef), est.ranef),
	est.copas = ifelse(outcome.type == "bin", exp(est.copas), est.copas),
	est.reg = ifelse(outcome.type == "bin", exp(est.reg.ranef), est.reg.ranef)
)

effect.diff <- meta.f %>% mutate(
	est.fixef = ifelse(outcome.type == "bin", est.fixef - 1, est.fixef),
	est.ranef = ifelse(outcome.type == "bin", est.ranef - 1, est.ranef),
	est.copas = ifelse(outcome.type == "bin", est.copas - 1, est.copas),
	est.reg = ifelse(outcome.type == "bin", est.reg.ranef - 1, est.reg.ranef),
	fixef.ranef = ifelse(abs(est.ranef) > abs(est.fixef), "Reduction", "Amplification"),
	fixef.copas = ifelse(abs(est.fixef) > abs(est.copas), "Reduction", "Amplification"),
	fixef.reg = ifelse(abs(est.fixef) > abs(est.reg), "Reduction", "Amplification"),
	fixef.copas.s = ifelse(sign(est.fixef) == sign(est.copas), "Unchanged", "Reversed"),
	fixef.reg.s = ifelse(sign(est.fixef) == sign(est.reg), "Unchanged", "Reversed")) 


#Scatterplots of effect.diff-analysis and corrected log(abs(abs estimates
ranef.fixef.sc.effect <- effect.diff %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.ranef)))) + geom_point(size = 0.8) + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	xlab("Fixed effects log estimate") + ylab("Random effects log estimate") + ggtitle("Fixed effects to Random effects")

# effect.diff %>% ggplot(aes(y = log(abs((est.ranef))), x = log(abs(est.trimfill.fixef)))) + geom_point() + 
# 	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
# 	xlab("Random effects log estimate") + ylab("Trimfill log estimate") + ggtitle("Random effects to Trimfill")


# trimfill.fixef.sc.effect <- effect.diff %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.trimfill.fixef)))) + geom_point(size = 0.8) + 
# 	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
# 	xlab("Fixed effects log estimate") + ylab("Trimfill log estimate") + ggtitle("Fixed effects to Trimfill")
# 
# trimfill.ranef.sc.effect <- effect.diff %>% ggplot(aes(x = log(abs((est.ranef))), y = log(abs(est.trimfill.fixef)))) + geom_point(size = 0.8) + 
# 	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
# 	xlab("Random effects log estimate") + ylab("Trimfill log estimate") + ggtitle("Random effects to Trimfill")

copas.fixef.sc.effect <- effect.diff %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.copas)))) + geom_point(size = 0.8) + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	xlab("Fixed effects log estimate") + ylab("Copas log estimate") + ggtitle("Fixed effects to Copas") 

copas.ranef.sc.effect <- effect.diff %>% ggplot(aes(x = log(abs((est.ranef))), y = log(abs(est.copas)))) + geom_point(size = 0.8) + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE, show.legend = T) + 
	xlab("Random effects log estimate") + ylab("Copas log estimate") + ggtitle("Random effects to Copas") 

reg.fixef.sc.effect <- effect.diff %>% ggplot(aes(x = log(abs((est.fixef))), y = log(abs(est.reg.ranef)))) + geom_point(size = 0.8) + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	xlab("Fixed effects log estimate") + ylab("Regression log estimate") + ggtitle("Fixed effects to Regression")

reg.ranef.sc.effect <- effect.diff %>% ggplot(aes(x = log(abs((est.ranef))), y = log(abs(est.reg.ranef)))) + geom_point(size = 0.8) + 
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + 
	xlab("Random effects log estimate") + ylab("Regression log estimate") + ggtitle("Random effects to Regression")
#--------------------------------------------------------------------------------------------------------------------#


#Scatterplots of test- statistics:
ranef.fixef.sc.zval <- effect.diff %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.ranef))) + geom_point(size = 0.8) +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Fixed effects and Random effects z statistics") + ylab("Random effects statistic") + xlab("Fixed effects statistic")

effect.diff %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.trimfill.fixef))) + geom_point(size = 0.8) +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Fixed effects and Random effects z statistics") + ylab("Random effects statistic") + xlab("Fixed effects statistic")

# trimfill.fixef.sc.zval <- effect.diff %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.trimfill.fixef))) + geom_point(size = 0.8) +
# 	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
# 	ggtitle("Fixed effects and trimfill z statistics") + ylab("Trimfill  statistic") + xlab("Fixed effects  statistic")
# 
# trimfill.ranef.sc.zval <- effect.diff %>% ggplot(aes(x = abs(zval.ranef), y = abs(zval.trimfill.fixef))) + geom_point(size = 0.8) +
# 	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
# 	ggtitle("Random effects and trimfill z statistics") + ylab("Trimfill statistic") + xlab("Random effects  statistic")

copas.fixef.sc.zval <- effect.diff %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.copas))) + geom_point(size = 0.8) +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + xlim(0, 25) +
	ggtitle("Fixed effects and Copas z statistics") + ylab("Copas  statistic") + xlab("Fixed effects  statistic")

copas.ranef.sc.zval <- effect.diff %>% ggplot(aes(x = abs(zval.ranef), y = abs(zval.copas))) + geom_point(size = 0.8) +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + xlim(0, 20) +
	ggtitle("Random effects and Copas z statistics") + ylab("Copas  statistic") + xlab("Random effects  statistic")

reg.fixef.sc.zval <- effect.diff %>% ggplot(aes(x = abs(zval.fixef), y = abs(zval.reg.ranef))) + geom_point(size = 0.8) +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Fixed effects and Regression z statistics") + ylab("Regression  statistic") + xlab("Fixed effects  statistic")

reg.ranef.sc.zval <- effect.diff %>% ggplot(aes(x = abs(zval.ranef), y = abs(zval.reg.ranef))) + geom_point(size = 0.8) +
	geom_abline(slope = 1, color = "red") + geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	ggtitle("Random effects and Regression z statistics") + ylab("Regression  statistic") + xlab("Random effects  statistic")
#--------------------------------------------------------------------------------------------------------------------#

#Missing study proportion plots:
p.missing.copas <- meta.f %>% ggplot(aes(x = miss.copas/n)) + geom_histogram(bins = 20) + 
	xlim(0, 0.5) + ylim(0, 200) + xlab("Proportion of missing studies") + ggtitle("Copas selection model") + theme_bw()
# p.missing.trim <- meta.f %>% ggplot(aes(x = missing.trim)) + geom_histogram(bins = 20) + 
# 	xlim(0, 0.5) + ylim(0, 200) + xlab("Proportion of missing studies") + ggtitle("Trim-and-fill") + theme_bw()
#--------------------------------------------------------------------------------------------------------------------#

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
	ggplot(aes(x = missing.trim.cont, fill = factor(thompson.test), stat(count))) + 
	geom_density(na.rm = T, position = "fill") + 
	labs(fill = "Significance") + scale_fill_discrete(labels= c("Yes", "No"))
#--------------------------------------------------------------------------------------------------------------------#

