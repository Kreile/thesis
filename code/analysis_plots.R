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


# data.ext2 <- pb.process2(data)

################################################################################################
################################################################################################
#PROPORTION OF SIGNIFICANT PUBLICATION BIAS TESTS
################################################################################################
################################################################################################

sig.level <- 0.1
meta.f <- meta.f %>% 
	mutate(egger.test = ifelse(pval.egger < sig.level, 1, 0),
				 thompson.test = ifelse(pval.thompson < sig.level, 1, 0),
				 begg.test = ifelse(pval.begg < sig.level, 1, 0),
				 
				 schwarzer.test = ifelse(pval.schwarzer < sig.level, 1, 0),
				 rucker.test = ifelse(pval.rucker < sig.level, 1, 0),
				 harbord.test = ifelse(pval.harbord < sig.level, 1, 0),
				 peter.test = ifelse(pval.peter < sig.level, 1, 0))

meta.bin <- meta.f %>% filter(outcome.type == "bin")
meta.cont <- meta.f %>% filter(outcome.type == "cont")
meta.surv <- meta.f %>% filter(outcome.type == "surv")

#Overall proportions, to be finsihed yet:
# test <- meta.f %>% group_by(outcome.type) %>% summarize(
#   schwarzer.test = mean(schwarzer.test),
#   rucker.test = mean(rucker.test),
#   harbord.test = mean(harbord.test),
#   peter.test = mean(peter.test),
#   egger.test = mean(egger.test),
#   begg.test = mean(begg.test),
#   thompson.test = mean(thompson.test))
# test <- test %>% gather(key = "test.type", value = "mean", schwarzer.test:thompson.test) 
# # test$mean[which(!is.na(test$mean))] <- round(test$mean[which(!is.na(test$mean))], 2)
# 
# #Compare proportions of total concordance vs. proportion of at least single positive test results among all
# meta.f %>% ungroup() %>% 
#   select(outcome.type, schwarzer.test, rucker.test, harbord.test, peter.test, egger.test, begg.test, thompson.test) %>% 
#   gather(key = "test.type", value = "null.hypothesis", schwarzer.test:thompson.test) %>%  
#   mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
#   ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
#   theme_bw() + xlab(label = NULL) + facet_wrap(~outcome.type, scales = "free") + theme(legend.position = "top") +
#   guides(fill=guide_legend(title=NULL)) 
# 
# # + annotate("text", x = test$test.type, y = 10,
# #            label = paste(round(test$mean, 2,)*100, "% rejected"),
# #            color = "white")




test.bin <- meta.bin %>% ungroup() %>% summarize(
																								 schwarzer.test = mean(schwarzer.test),
																								 rucker.test = mean(rucker.test),
																								 harbord.test = mean(harbord.test),
																								 peter.test = mean(peter.test))
test.bin <- test.bin %>% gather(key = "test.type", value = "mean")

p.bin <- meta.bin %>% ungroup() %>% 
	select(schwarzer.test, rucker.test, harbord.test, peter.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>%  
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Binary Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.bin$test.type, y = 500, 
					 label = paste(round(test.bin$mean, 2)*100, "% rejected"), 
					 color = "white")
#--------------------------------------------------------------------------------------------------------------------#

dat_text <- data.frame(
	label = paste(c(sum(meta.bin$harbord.test), sum(meta.bin$peter.test), 
									sum(meta.bin$rucker.test), sum(meta.bin$schwarzer.test)), "< 0.1,",
								round(c(mean(meta.bin$harbord.test), mean(meta.bin$peter.test), 
												mean(meta.bin$rucker.test), mean(meta.bin$schwarzer.test)),2)*100, "%"),
	test.type   = c("pval.harbord", "pval.peter", "pval.rucker", "pval.schwarzer"))

labels <- c(pval.harbord = "Harbord", pval.peter = "Peter", pval.rucker = "Rucker", pval.schwarzer = "Schwarzer")

p.dist.bin <- meta.bin %>% ungroup() %>% 
	select(pval.harbord, pval.peter, pval.rucker, pval.schwarzer) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(boundary = 0, binwidth = 0.1) + theme_bw() + facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 170, label = label), color = "black") + 
	theme(strip.text.x = element_text(size=15))
#--------------------------------------------------------------------------------------------------------------------#

#Survival Outcomes:
test.cont <- meta.cont %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																									 begg.test = mean(begg.test),
																									 thompson.test = mean(thompson.test))

test.cont <- test.cont %>% gather(key = "test.type", value = "mean")

p.cont <- meta.cont %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Continuous Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.cont$test.type, y = 150, 
					 label = paste(round(test.cont$mean, 2)*100, "% rejected"), 
					 color = "white")
#--------------------------------------------------------------------------------------------------------------------#

dat_text <- data.frame(
	label = paste(c(sum(meta.cont$begg.test), sum(meta.cont$egger.test), sum(meta.cont$thompson.test)), "< 0.1,",
								round(c(mean(meta.cont$begg.test), mean(meta.cont$egger.test), mean(meta.cont$thompson.test)),2)*100, "%"),
	test.type   = c("pval.begg", "pval.egger", "pval.thompson"))

labels <- c(pval.begg = "Begg Mazumdar", pval.egger = "Egger", pval.thompson = "Thompson Sharp")
p.dist.cont <- meta.cont %>% ungroup() %>% 
	select(pval.egger, pval.thompson, pval.begg) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(boundary = 0, binwidth = 0.1) + theme_bw() + 
	facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 50, label = label),  color = "black") + 
	theme(strip.text.x = element_text(size=15))
#--------------------------------------------------------------------------------------------------------------------#

#Survival Outcomes:
test.surv <- meta.surv %>% ungroup() %>% summarize(egger.test = mean(egger.test),
																									 begg.test = mean(begg.test),
																									 thompson.test = mean(thompson.test))

test.surv <- test.surv %>% gather(key = "test.type", value = "mean")

p.surv <- meta.surv %>% ungroup() %>% 
	select(egger.test, thompson.test, begg.test) %>% 
	gather(key = "test.type", value = "null.hypothesis") %>% 
	mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
	ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
	theme_bw() + xlab(label = NULL) + ggtitle("Continuous Outcomes") + theme(legend.position = "top") +
	guides(fill=guide_legend(title=NULL))+
	annotate("text", x = test.surv$test.type, y = 10, 
					 label = paste(round(test.surv$mean, 2)*100, "% rejected"), 
					 color = "white")
#--------------------------------------------------------------------------------------------------------------------#

dat_text <- data.frame(
	label = paste(c(sum(meta.surv$begg.test), sum(meta.surv$egger.test), sum(meta.surv$thompson.test)), 
	              "< 0.1,",
								round(c(mean(meta.surv$begg.test), mean(meta.surv$egger.test), mean(meta.surv$thompson.test)),2)*100, 
								"%"),
	test.type   = c("pval.begg", "pval.egger", "pval.thompson"))

labels <- c(pval.begg = "Begg Mazumdar", pval.egger = "Egger", pval.thompson = "Thompson Sharp")
p.dist.surv <- meta.surv %>% ungroup() %>% 
	select(pval.egger, pval.thompson, pval.begg) %>% 
	gather(key = "test.type", value = "p.value") %>% 
	ggplot(aes(x = p.value)) + geom_histogram(boundary = 0, binwidth = 0.1) + theme_bw() + 
	facet_wrap(~test.type, labeller = labeller(test.type = labels)) + 
	geom_text(data = dat_text, mapping = aes(x = 0.5, y = 13, label = label),  color = "black") + 
	theme(strip.text.x = element_text(size=15))
#--------------------------------------------------------------------------------------------------------------------#

#Test agreement
agree.bin <- meta.bin %>% mutate(n.sig = peter.test + rucker.test + harbord.test + schwarzer.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + 
  xlab("Number of significant tests") + ylab("count") + ggtitle("Binary Outcomes")
#--------------------------------------------------------------------------------------------------------------------#

agree.cont <- meta.cont %>% mutate(n.sig = egger.test + thompson.test + begg.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + 
  xlab("Number of significant tests") + ylab("count") + ggtitle("Continuous Outcomes")
#--------------------------------------------------------------------------------------------------------------------#

agree.surv <- meta.surv %>% mutate(n.sig = egger.test + thompson.test + begg.test) %>% 
	group_by(n.sig) %>% count %>% filter(n.sig > 0) %>% 
	ggplot(aes(y = nn, x = n.sig)) + theme_bw() + geom_col() + 
  xlab("Number of significant tests") + ylab("count") + ggtitle("Continuous Outcomes")
#--------------------------------------------------------------------------------------------------------------------#

################################################################################################
################################################################################################
#CHANGE OF P VALUES AFTER ADJUSTMENT
################################################################################################
################################################################################################

#Comparison of z-scores:
sig.zcor <- meta.f %>% mutate(z.fixef = est.z.fixef/se.est.z.fixef,
                              z.ranef = est.z.ranef/se.est.z.ranef,
                              z.reg = est.z.reg/se.est.z.reg,
                              z.copas = est.z.copas/se.est.z.copas) %>%  
  select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
  mutate(p.fisher.z = 2*(1-pnorm(abs(fisher.z)))) %>% 
  group_by(method) %>% 
  summarise(significant = length(which(abs(fisher.z) > 1.96)),
            p.significant = significant/length(fisher.z))
sig.zcor <- data.frame(method = sig.zcor$method,
                       label = paste(round(sig.zcor$p.significant,3)*100, "% significant", ", (n = ", sig.zcor$significant, ")", sep = ""))

method_names <- c(
  z.fixef = "Fixed Effects",
  z.ranef = "Random Effects",
  z.reg = "Regression",
  z.copas = "Copas"
)

adjustment.p.z <- meta.f %>% 
  mutate(z.fixef = est.z.fixef/se.est.z.fixef,
         z.ranef = est.z.ranef/se.est.z.ranef,
         z.reg = est.z.reg/se.est.z.reg,
         z.copas = est.z.copas/se.est.z.copas) %>%  
  select(z.fixef, z.ranef, z.reg, z.copas) %>% 
  gather(key = "method", value = "fisher.z") %>% 
  mutate(method = factor(method, levels = c("z.fixef", "z.ranef", "z.reg", "z.copas"))) %>% 
  mutate(p.fisher.z = 2*(1-pnorm(abs(fisher.z)))) %>% 
  ggplot(aes(x = p.fisher.z)) + 
  geom_histogram(boundary = 0, bins = 20) + 
  theme_bw() + 
  facet_wrap(~method, ncol = 2, labeller = as_labeller(method_names)) +
  geom_text(data = sig.zcor, aes(x = 0.5, y = 750, label = label), position = "dodge") + 
  xlab(expression(paste(italic(p),"-value of test statistic of the ", italic(z),"-score"))) +
  ggtitle(expression(paste(italic(z), " Score ", italic(p),"-value of Meta-Analysis and Adjustment Method")))
#--------------------------------------------------------------------------------------------------------------------#

#Comparison of SMD's:
sig.d <- meta.f %>% mutate(d.fixef = est.d.fixef/se.est.d.fixef,
                           d.ranef = est.d.ranef/se.est.d.ranef,
                           d.reg = est.d.reg/se.est.d.reg,
                           d.copas = est.d.copas/se.est.d.copas) %>%  
  select(d.fixef, d.ranef, d.reg, d.copas) %>% gather(key = "method", value = "smd") %>% 
  mutate(p.smd = 2*(1-pnorm(abs(smd)))) %>% 
  group_by(method) %>% 
  summarise(significant = length(which(abs(smd) > 1.96)),
            p.significant = significant/length(smd))
sig.d <- data.frame(method = sig.d$method,
                    label = paste(round(sig.d$p.significant,3)*100, "% significant", ", (n = ", sig.d$significant, ")", sep = ""))

method_names <- c(
  d.fixef = "Fixed Effects",
  d.ranef = "Random Effects",
  d.reg = "Regression",
  d.copas = "Copas"
)

adjustment.p.d <- meta.f %>% 
  mutate(d.fixef = est.d.fixef/se.est.d.fixef,
         d.ranef = est.d.ranef/se.est.d.ranef,
         d.reg = est.d.reg/se.est.d.reg,
         d.copas = est.d.copas/se.est.d.copas) %>%  
  select(d.fixef, d.ranef, d.reg, d.copas) %>% 
  gather(key = "method", value = "smd") %>% 
  mutate(method = factor(method, levels = c("d.fixef", "d.ranef", "d.reg", "d.copas"))) %>% 
  mutate(p.smd = 2*(1-pnorm(abs(smd)))) %>% 
  ggplot(aes(x = p.smd)) + 
  geom_histogram(boundary = 0, bins = 20) + 
  theme_bw() + 
  facet_wrap(~method, ncol = 2, labeller = as_labeller(method_names)) +
  geom_text(data = sig.d, aes(x = 0.5, y = 750, label = label), position = "dodge") + 
  xlab(expression(paste(italic(p),"-value of test statistic of Hedges ", italic(g)))) +
  ggtitle(expression(paste("Hedges ", italic(g), " ", italic(p),"-value of Meta-Analysis and Adjustment Method")))
#--------------------------------------------------------------------------------------------------------------------#

#Comparison of log hazard ratios:
sig.d <- meta.f %>% filter(outcome.type == "surv") %>% 
  mutate(d.fixef = est.fixef/se.est.fixef,
         d.ranef = est.ranef/se.est.ranef,
         d.reg = est.reg/se.est.reg,
         d.copas = est.copas/se.est.copas) %>%  
  select(d.fixef, d.ranef, d.reg, d.copas) %>% gather(key = "method", value = "smd") %>% 
  group_by(method) %>% 
  summarise(significant = length(which(abs(smd) > 1.96)),
            p.significant = (significant + sum(is.na(smd)))/(length(smd)))
sig.d <- data.frame(method = sig.d$method,
                    label = paste(round(sig.d$p.significant,3)*100, "% significant", ", (n = ", sig.d$significant, ")", sep = ""))

method_names <- c(
  d.fixef = "Fixed Effects",
  d.ranef = "Random Effects",
  d.reg = "Regression",
  d.copas = "Copas"
)

adjustment.p.log.hazard.ratio <- meta.f %>% filter(outcome.type == "surv") %>% 
  mutate(d.fixef = est.fixef/se.est.fixef,
         d.ranef = est.ranef/se.est.ranef,
         d.reg = est.reg/se.est.reg,
         d.copas = est.copas/se.est.copas) %>%  
  select(d.fixef, d.ranef, d.reg, d.copas) %>% gather(key = "method", value = "smd") %>% 
  mutate(method = factor(method, levels = c("d.fixef", "d.ranef", "d.reg", "d.copas"))) %>% 
  mutate(p.smd = 2*(1-pnorm(abs(smd)))) %>% 
  ggplot(aes(x = p.smd)) + 
  geom_histogram(boundary = 0, bins = 20) + 
  theme_bw() + 
  facet_wrap(~method, ncol = 2, labeller = as_labeller(method_names)) +
  geom_text(data = sig.d, aes(x = 0.0002, y = 30, label = label), position = "dodge") + 
  xlab(expression(paste(italic(p),"-value of test statistic of Hedges ", italic(g)))) +
  ggtitle(expression(paste("Log hazard ratio", italic(p),"-value of Meta-Analysis and Adjustment Method")))
#--------------------------------------------------------------------------------------------------------------------#


################################################################################################
################################################################################################
#CHANGE OF TEST STATISTICS AFTER ADJUSTMENT
################################################################################################
################################################################################################


#HISTOGRAMS: 
#Comparison of z-scores:
sig.zcor <- meta.f %>% mutate(z.fixef = est.z.fixef/se.est.z.fixef,
																 z.ranef = est.z.ranef/se.est.z.ranef,
																 z.reg = est.z.reg/se.est.z.reg,
																 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	group_by(method) %>% 
	summarise(significant = length(which(abs(fisher.z) > 1.96)),
						p.significant = significant/length(fisher.z))
sig.zcor <- data.frame(method = sig.zcor$method,
											 label = paste(round(sig.zcor$p.significant,3)*100, "% significant", ", (n = ", sig.zcor$significant, ")", sep = ""))

method_names <- c(
  z.fixef = "Fixed Effects",
  z.ranef = "Random Effects",
  z.reg = "Regression",
  z.copas = "Copas"
)

meta.f %>% mutate(z.fixef = est.z.fixef/se.est.z.fixef,
										 z.ranef = est.z.ranef/se.est.z.ranef,
										 z.reg = est.z.reg/se.est.z.reg,
										 z.copas = est.z.copas/se.est.z.copas) %>%  
	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
	# filter(abs(fisher.z) < 10) %>% 
	ggplot(aes(x = abs(fisher.z))) + 
  geom_histogram(boundary = 0, binwidth = 1) + 
  theme_bw() + facet_wrap(~method, ncol = 2, labeller = as_labeller(method_names)) + 
	geom_vline(xintercept = 1.96, color = "red") + xlim(c(0,100)) +
	geom_text(data = sig.zcor, aes(x = 50, y = 1000, label = label), position = "dodge")
#--------------------------------------------------------------------------------------------------------------------#

#Comparison of SMD's:
sig.d <- meta.f %>% mutate(d.fixef = est.d.fixef/se.est.d.fixef,
                              d.ranef = est.d.ranef/se.est.d.ranef,
                              d.reg = est.d.reg/se.est.d.reg,
                              d.copas = est.d.copas/se.est.d.copas) %>%  
  select(d.fixef, d.ranef, d.reg, d.copas) %>% gather(key = "method", value = "smd") %>% 
  group_by(method) %>% 
  summarise(significant = length(which(abs(smd) > 1.96)),
            p.significant = significant/length(smd))
sig.d <- data.frame(method = sig.d$method,
                       label = paste(round(sig.d$p.significant,3)*100, "% significant", ", (n = ", sig.d$significant, ")", sep = ""))

method_names <- c(
  d.fixef = "Fixed Effects",
  d.ranef = "Random Effects",
  d.reg = "Regression",
  d.copas = "Copas"
)

meta.f %>% mutate(d.fixef = est.d.fixef/se.est.d.fixef,
                  d.ranef = est.d.ranef/se.est.d.ranef,
                  d.reg = est.d.reg/se.est.d.reg,
                  d.copas = est.d.copas/se.est.d.copas) %>%  
  select(d.fixef, d.ranef, d.reg, d.copas) %>% gather(key = "method", value = "smd") %>% 
  filter(abs(smd) < 100) %>% 
  ggplot(aes(x = abs(smd))) + geom_histogram(boundary = 0, binwidth = 1) + 
  theme_bw() + facet_wrap(~method, ncol = 2, labeller = as_labeller(method_names)) + 
  geom_vline(xintercept = 1.96, color = "red") + #xlim(c(0,100)) +
  geom_text(data = sig.d, aes(x = 50, y = 600, label = label), position = "dodge") +
  scale_colour_manual(name="Line Color",
                      values=c(myline1="red"))
#--------------------------------------------------------------------------------------------------------------------#

#Comparison of log hazard ratios:
sig.d <- meta.f %>% filter(outcome.type == "surv") %>% 
  mutate(d.fixef = est.fixef/se.est.fixef,
                           d.ranef = est.ranef/se.est.ranef,
                           d.reg = est.reg/se.est.reg,
                           d.copas = est.copas/se.est.copas) %>%  
  select(d.fixef, d.ranef, d.reg, d.copas) %>% gather(key = "method", value = "smd") %>% 
  group_by(method) %>% 
  summarise(significant = length(which(abs(smd) > 1.96)),
            p.significant = significant/length(smd))
sig.d <- data.frame(method = sig.d$method,
                    label = paste(round(sig.d$p.significant,3)*100, "% significant", ", (n = ", sig.d$significant, ")", sep = ""))

method_names <- c(
  d.fixef = "Fixed Effects",
  d.ranef = "Random Effects",
  d.reg = "Regression",
  d.copas = "Copas"
)

meta.f %>% filter(outcome.type == "surv") %>% 
  mutate(d.fixef = est.fixef/se.est.fixef,
         d.ranef = est.ranef/se.est.ranef,
         d.reg = est.reg/se.est.reg,
         d.copas = est.copas/se.est.copas) %>%  
  select(d.fixef, d.ranef, d.reg, d.copas) %>% gather(key = "method", value = "smd") %>% 
  filter(abs(smd) < 100) %>% 
  ggplot(aes(x = abs(smd))) + geom_histogram(boundary = 0, binwidth = 1) + 
  theme_bw() + facet_wrap(~method, ncol = 2, labeller = as_labeller(method_names)) + 
  geom_vline(xintercept = 1.96, color = "red") + #xlim(c(0,100)) +
  geom_text(data = sig.d, aes(x = 32, y = 10, label = label), position = "dodge") +
  scale_colour_manual(name="Line Color",
                      values=c(myline1="red"))
#--------------------------------------------------------------------------------------------------------------------#


#Meta-Analysis and adjusted treatment effect estimate difference:

#Z-score:
diff.z.fixef <- meta.f %>% 
  mutate(fixef = est.z.fixef, ranef = est.z.ranef,
         copas = est.z.copas, regression = est.z.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(fixef) ~  abs(copas),
                           sign(copas) != sign(fixef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(fixef) ~  abs(regression),
                                sign(regression) != sign(fixef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  select(fixef, ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.z.score", copas:regression) %>% 
  mutate(difference = fixef - adjusted.z.score) %>% 
  filter(difference > -.75 & difference < .75) %>% 
  ggplot(aes(x = difference)) + geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + 
  xlab(expression(paste("Fixed effects - adjusted ", z, "-score")))
#--------------------------------------------------------------------------------------------------------------------#

diff.z.ranef <- meta.f %>% 
  mutate(fixef = est.z.fixef, ranef = est.z.ranef,
         copas = est.z.copas, regression = est.z.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(ranef) ~  abs(copas),
                           sign(copas) != sign(ranef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(ranef) ~  abs(regression),
                                sign(regression) != sign(ranef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  select(fixef, ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.z.score", copas:regression) %>% 
  mutate(difference = ranef - adjusted.z.score) %>% 
  filter(difference > -.75 & difference < .75) %>% 
  ggplot(aes(x = difference)) + geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + 
  xlab(expression(paste("Random effects - adjusted ", z, "-score")))
#--------------------------------------------------------------------------------------------------------------------#

#SMD:
diff.d.fixef <- meta.f %>% 
  mutate(fixef = est.d.fixef, ranef = est.d.ranef,
         copas = est.d.copas, regression = est.d.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(fixef) ~  abs(copas),
                           sign(copas) != sign(fixef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(fixef) ~  abs(regression),
                                sign(regression) != sign(fixef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  select(fixef, ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.smd", copas:regression) %>% 
  mutate(difference = fixef - adjusted.smd) %>% 
  filter(difference > -1 & difference < 1) %>% 
  ggplot(aes(x = difference)) + geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + 
  xlab(expression(paste("Fixed effects - adjusted Hedges ", g)))
#--------------------------------------------------------------------------------------------------------------------#

diff.d.ranef <- meta.f  %>% 
  mutate(fixef = est.d.fixef, ranef = est.d.ranef,
         copas = est.d.copas, regression = est.d.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(ranef) ~  abs(copas),
                           sign(copas) != sign(ranef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(ranef) ~  abs(regression),
                                sign(regression) != sign(ranef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  gather(key = "method", value = "adjusted.smd", copas:regression) %>% 
  mutate(difference = ranef - adjusted.smd) %>% 
  filter(difference > -1 & difference < 1) %>% 
  ggplot(aes(x = difference)) + geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + 
  xlab(expression(paste("Random effects - adjusted Hedges ", g)))
#--------------------------------------------------------------------------------------------------------------------#

#Log hazard ratios:
diff.log.hazard.ratio.fixef <- meta.surv %>% 
  mutate(fixef = est.fixef, ranef = est.ranef,
         copas = est.copas, regression = est.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(fixef) ~  abs(copas),
                           sign(copas) != sign(fixef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(fixef) ~  abs(regression),
                                sign(regression) != sign(fixef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  gather(key = "method", value = "adjusted.log.hazard.ratio", copas:regression)  %>% 
  ggplot(aes(x = (fixef - adjusted.log.hazard.ratio))) + 
  geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + xlab("Fixed effects - adjusted log hazard ratio")
#--------------------------------------------------------------------------------------------------------------------#

diff.log.hazard.ratio.ranef <- meta.surv %>% 
  mutate(fixef = est.fixef, ranef = est.ranef,
         copas = est.copas, regression = est.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(ranef) ~  abs(copas),
                           sign(copas) != sign(ranef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(ranef) ~  abs(regression),
                                sign(regression) != sign(ranef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  gather(key = "method", value = "adjusted.log.hazard.ratio", copas:regression)  %>% 
  ggplot(aes(x = (ranef - adjusted.log.hazard.ratio))) + 
  geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + xlab("Random effects - adjusted log hazard ratio")
#--------------------------------------------------------------------------------------------------------------------#

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

effect.diff <- effect.diff %>% 
  mutate(change.fixef.copas = ifelse(fixef.copas == "Reduction", -abs(est.fixef - est.copas),
                                     abs(est.fixef - est.copas)),
         change.fixef.reg = ifelse(fixef.ranef == "Reduction", -abs(est.fixef - est.reg),
                                   abs(est.fixef - est.reg)))

# hist.ranef <- effect.diff %>% ggplot(aes(x = abs(est.fixef - est.ranef))) +
# 	geom_histogram() + facet_grid(~factor(fixef.ranef)) + xlim(0, 2) + xlab("Effect difference") + 
# 	theme_bw() + theme(legend.position = "bottom", strip.text.x = element_text(size=15)) + ylim(0, 250)

#Copas selection model:
copas.miss <- effect.diff %>% group_by(outcome.type) %>% mutate(dif = abs(est.fixef - est.copas)) %>% 
	filter(dif > 2) %>% mutate(label = round(dif, 1)) %>% arrange(label) %>% select(label, meta.id, outcome.measure.new)
copas.label <- data.frame(label = paste("missing:", c(paste(copas.miss$label[copas.miss$outcome.type == "bin"], collapse = ", "),
																											paste(copas.miss$label[copas.miss$outcome.type == "cont"], collapse = ", "),
																											paste(copas.miss$label[copas.miss$outcome.type == "surv"], collapse = ", "))),
													outcome.type = c("bin", "cont", "surv"), 
													outcome.measure.new = c("Risk Ratio", "Mean Difference", "Hazard Ratio"))
copas.label$label <- as.character(copas.label$label)
copas.label$label[2] <- sub("2.8,", "2.8,\n", copas.label$label[2])

effect.diff %>% ggplot(aes(x = change.fixef.copas, fill = outcome.measure.new)) + 
	facet_wrap(outcome.type~., scales = "free", labeller = labeller(outcome.type = NULL)) +
	geom_histogram() + xlim(-2, 2) + xlab("Effect change") + 
	theme_bw() + theme(legend.position = "bottom", strip.text = element_blank(), legend.title = element_blank())  +
	geom_text(data = copas.label, 
						mapping = aes(x = 0, y = c(220, 47.4, 10), label = label, vjust = "center"), position = "dodge") #+ ylim(0, 250)
#--------------------------------------------------------------------------------------------------------------------#

#Regression model:
reg.miss <- effect.diff %>% group_by(outcome.type) %>% mutate(dif = abs(est.fixef - est.reg)) %>% 
	filter(dif > 2) %>% mutate(label = round(dif, 1)) %>% arrange(label) %>% select(label, meta.id, outcome.measure.new)
reg.label <- data.frame(label = paste("missing:", c(paste(reg.miss$label[reg.miss$outcome.type == "bin"], collapse = ", "),
																										paste(reg.miss$label[reg.miss$outcome.type == "cont"], collapse = ", "),
																										paste(reg.miss$label[reg.miss$outcome.type == "surv"], collapse = ", "))),
												outcome.type = c("bin", "cont", "surv"), 
												outcome.measure.new = c("Risk Ratio", "Mean Difference", "Hazard Ratio"))
reg.label$label <- as.character(reg.label$label)
reg.label$label[2] <- sub("2.2,","2.2,\n ", reg.label$label[2])
reg.label$label[2] <- sub("2.6,","2.6,\n ", reg.label$label[2])
reg.label$label[2] <- sub(", 3.1,",", 3.1,\n ", reg.label$label[2])
reg.label$label[2] <- sub("4.4, 4.4,", "4.4, 4.4,\n", reg.label$label[2])
reg.label$label[2] <- sub(", 9.5,", ", 9.5,\n", reg.label$label[2])

reg.label$label[1] <- sub(", 3.3,", ", 3.3,\n", reg.label$label[1])
reg.label$label[1] <- sub(", 7.6,", ", 7.6,\n", reg.label$label[1])

effect.diff %>% ggplot(aes(x = est.fixef - est.reg, fill = outcome.measure.new)) + 
	facet_wrap(outcome.type~., scales = "free", labeller = labeller(outcome.type = NULL)) +
	geom_histogram() + xlim(-2, 2) + xlab("Effect change") + 
	theme_bw() + theme(legend.position = "bottom", strip.text = element_blank(), legend.title = element_blank())  +
	geom_text(data = reg.label, 
						mapping = aes(x = 0, y = c(330, 50, 10), label = label, vjust = "center"), position = "dodge") #+ ylim(0, 250)


################################################################################################
################################################################################################
#CHANGE OF EFFECT SIZE AFTER ADJUSTMENT
################################################################################################
################################################################################################

#Standardized mean differences (hedges g and cohens d):
number.outside.smd <- meta.f %>% mutate(copas = est.d.copas, regression = est.d.reg) %>% 
  select(meta.id, est.d.fixef, est.d.ranef, copas, regression) %>% 
  gather(key = "method", value = "smd", est.d.fixef:regression) %>% filter(smd > 4)

meta.f %>% mutate(copas = est.d.copas, regression = est.d.reg) %>% 
  select(est.d.fixef, est.d.ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.smd", copas:regression) %>% 
  ggplot(aes(y = abs(est.d.fixef), x = abs(adjusted.smd))) + geom_point(alpha = 0.4, size = 0.6) + 
  facet_wrap(~method, scales = "free") + theme_bw() + 
  xlab("Adjusted SMD") + ylab("Fixed Effects SMD") + xlim(c(0,4)) + ylim(c(0,4)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Fixed effects vs. Adjusted Std. Mean Difference", subtitle = "with linear regression fit")

meta.f %>%   mutate(copas = est.d.copas, regression = est.d.reg) %>% 
  select(est.d.fixef, est.d.ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.smd", copas:regression) %>% 
  ggplot(aes(y = abs(est.d.ranef), x = abs(adjusted.smd))) + 
  geom_point(alpha = 0.4, size = 0.6) + 
  facet_wrap(~method, scales = "free") + theme_bw() + 
  xlab("Adjusted SMD") + ylab("Random Effects SMD") + xlim(c(0,4)) + ylim(c(0,4)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Random effects vs. Adjusted Std. Mean Difference", subtitle = "with linear regression fit")
#--------------------------------------------------------------------------------------------------------------------#

#Z-scores:
number.outside.z <- meta.f %>% mutate(copas = est.z.copas, regression = est.z.reg) %>% 
  select(meta.id, est.z.fixef, est.z.ranef, copas, regression) %>% 
  gather(key = "method", value = "smd", est.z.fixef:regression) %>% filter(smd > 4)

meta.f %>% mutate(copas = est.z.copas, regression = est.z.reg) %>% 
  select(est.z.fixef, est.z.ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.z.score", copas:regression) %>% 
  ggplot(aes(y = abs(est.z.fixef), x = abs(adjusted.z.score))) + geom_point(alpha = 0.4, size = 0.6) + 
  facet_wrap(~method, scales = "free") + theme_bw() + 
  xlab("Adjusted z-score") + ylab("Fixed Effects z-score") + xlim(c(0,1.5)) + ylim(c(0,1.5)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Fixed effects vs. Adjusted z-score", subtitle = "with linear regression fit")


meta.f %>% mutate(copas = est.z.copas, regression = est.z.reg) %>% 
  select(est.z.fixef, est.z.ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.z.score", copas:regression) %>% 
  ggplot(aes(y = abs(est.z.ranef), x = abs(adjusted.z.score))) + 
  geom_point(alpha = 0.4, size = 0.6) + 
  facet_wrap(~method, scales = "free") + theme_bw() + 
  xlab("Adjusted z-score") + ylab("Random Effects z-score") + xlim(c(0,1.5)) + ylim(c(0,1.5)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Random effects vs. Adjusted z-score", subtitle = "with linear regression fit")
#--------------------------------------------------------------------------------------------------------------------#

#log Hazard Ratios:
meta.surv %>% mutate(copas = est.copas, regression = est.reg) %>% 
  select(est.fixef, est.ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.log.hazard.ratio", copas:regression) %>% 
  ggplot(aes(y = abs(est.fixef), x = abs(adjusted.log.hazard.ratio))) + 
  geom_point(alpha = 0.4, size = 0.6) + 
  facet_wrap(~method, scales = "free") + theme_bw() + 
  xlab("Adjusted log hazard ratio") + ylab("Fixed Effects log hazard ratio") + 
  xlim(c(0.7,1.3)) + ylim(c(0.7,1.3)) + 
  geom_smooth(method = "lm", se = FALSE, ) + 
  ggtitle("Fixed effects vs. Adjusted log hazard ratio", subtitle = "with linear regression fit")

meta.surv %>% mutate(copas = est.copas, regression = est.reg) %>% 
  select(est.fixef, est.ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.log.hazard.ratio", copas:regression)  %>% 
  ggplot(aes(y = abs(est.ranef), x = abs(adjusted.log.hazard.ratio))) + 
  geom_point(alpha = 0.4, size = 0.6) + 
  facet_wrap(~method, scales = "free") + theme_bw() + 
  xlab("Adjusted log hazard ratio") + ylab("Random Effects log hazard ratio") + 
  xlim(c(0.7,1.3)) + ylim(c(0.7,1.3)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  ggtitle("Random effects vs. Adjusted log hazard ratio", subtitle = "with linear regression fit")

################################################################################################
################################################################################################
#CHANGE OF SIGNIFICANT TREATMENT EFFECT AFTER ADJUSTMENT
################################################################################################
################################################################################################

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


################################################################################################
################################################################################################
#PUBLICATION BIAS TEST RESULTS ANALYSIS
################################################################################################
################################################################################################

#Agreement proportions of publication bias tests:

#Binary:
meta.bin <- meta.bin %>%ungroup() %>%  mutate(tes.d.schwarzer = ifelse(tes.d.test == schwarzer.test, "agree", "disagree"),
                                              tes.d.peter = ifelse(tes.d.test == peter.test, "agree", "disagree"),
                                              tes.d.rucker = ifelse(tes.d.test == rucker.test, "agree", "disagree"),
                                              tes.d.harbord = ifelse(tes.d.test == harbord.test, "agree", "disagree"),
                                              schwarzer.peter = ifelse(schwarzer.test == peter.test, "agree", "disagree"),
                                              schwarzer.rucker = ifelse(schwarzer.test == rucker.test, "agree", "disagree"),
                                              schwarzer.harbord = ifelse(schwarzer.test == harbord.test, "agree", "disagree"),
                                              rucker.peter = ifelse(rucker.test == peter.test, "agree", "disagree"),
                                              rucker.harbord = ifelse(rucker.test == harbord.test, "agree", "disagree"),
                                              harbord.peter = ifelse(harbord.test == peter.test, "agree", "disagree"))
agreement.bin <- meta.bin %>% ungroup() %>% summarise(tes.d.schwarzer = sum(tes.d.schwarzer == "agree")/n(),
                                                      tes.d.peter = sum(tes.d.peter == "agree")/n(),
                                                      tes.d.rucker = sum(tes.d.rucker == "agree")/n(),
                                                      tes.d.harbord = sum(tes.d.harbord == "agree")/n(),
                                                      schwarzer.peter = sum(schwarzer.peter == "agree")/n(),
                                                      schwarzer.rucker = sum(schwarzer.rucker == "agree")/n(),
                                                      schwarzer.harbord = sum(schwarzer.harbord == "agree")/n(),
                                                      rucker.peter = sum(rucker.peter == "agree")/n(),
                                                      harbord.peter = sum(harbord.peter == "agree")/n())
correlation.bin <- meta.bin %>%ungroup() %>%  summarise(tes.d.schwarzer = cor(stat.d.tes, stat.schwarzer),
                                                        tes.d.peter = cor(stat.d.tes, stat.peter),
                                                        tes.d.rucker = cor(stat.d.tes, stat.rucker),
                                                        tes.d.harbord = cor(stat.d.tes, stat.harbord),
                                                        schwarzer.peter = cor(stat.schwarzer, stat.peter),
                                                        schwarzer.rucker = cor(stat.schwarzer, stat.rucker),
                                                        schwarzer.harbord = cor(stat.schwarzer, stat.harbord),
                                                        rucker.peter = cor(stat.rucker, stat.peter),
                                                        harbord.peter = cor(stat.harbord, stat.peter))
# rsquared.bin <- meta.bin %>%ungroup() %>%  summarise(tes.d.schwarzer = summary(lm(stat.d.tes~ stat.schwarzer))$r.squared,
#                                                         tes.d.peter = summary(lm(stat.d.tes~ stat.peter))$r.squared,
#                                                         tes.d.rucker = summary(lm(stat.d.tes~ stat.rucker))$r.squared,
#                                                         tes.d.harbord = summary(lm(stat.d.tes~ stat.harbord))$r.squared,
#                                                         schwarzer.peter = summary(lm(stat.schwarzer~ stat.peter))$r.squared,
#                                                         schwarzer.rucker = summary(lm(stat.schwarzer~ stat.rucker))$r.squared,
#                                                         schwarzer.harbord = summary(lm(stat.schwarzer~ stat.harbord))$r.squared,
#                                                         rucker.peter = summary(lm(stat.rucker~ stat.peter))$r.squared,
#                                                         harbord.peter = summary(lm(stat.harbord~ stat.peter))$r.squared)
meta.bin <- meta.bin %>%ungroup() %>%  mutate(tes.d.schwarzer = ifelse(tes.d.test + schwarzer.test > 1, "agree", "disagree"),
                                              tes.d.peter = ifelse(tes.d.test + peter.test > 1, "agree", "disagree"),
                                              tes.d.rucker = ifelse(tes.d.test + rucker.test > 1, "agree", "disagree"),
                                              tes.d.harbord = ifelse(tes.d.test + harbord.test > 1, "agree", "disagree"),
                                              schwarzer.peter = ifelse(schwarzer.test + peter.test > 1, "agree", "disagree"),
                                              schwarzer.rucker = ifelse(schwarzer.test + rucker.test > 1, "agree", "disagree"),
                                              schwarzer.harbord = ifelse(schwarzer.test + harbord.test > 1, "agree", "disagree"),
                                              rucker.peter = ifelse(rucker.test + peter.test > 1, "agree", "disagree"),
                                              rucker.harbord = ifelse(rucker.test + harbord.test > 1, "agree", "disagree"),
                                              harbord.peter = ifelse(harbord.test + peter.test > 1, "agree", "disagree"))
agreement.bin <- meta.bin %>% ungroup() %>% summarise(tes.d.schwarzer = sum(tes.d.schwarzer == "agree")/sum(schwarzer.test),
                                                      tes.d.peter = sum(tes.d.peter == "agree")/sum(peter.test),
                                                      tes.d.rucker = sum(tes.d.rucker == "agree")/n(),
                                                      tes.d.harbord = sum(tes.d.harbord == "agree")/n(),
                                                      schwarzer.peter = sum(schwarzer.peter == "agree")/n(),
                                                      schwarzer.rucker = sum(schwarzer.rucker == "agree")/n(),
                                                      schwarzer.harbord = sum(schwarzer.harbord == "agree")/n(),
                                                      rucker.peter = sum(rucker.peter == "agree")/n(),
                                                      harbord.peter = sum(harbord.peter == "agree")/n())

binary.tests.agreement <- rbind(agreement.bin, correlation.bin, rsquared.bin)
rownames(binary.tests.agreement) <- c("Test Agreement","p-value Correlation", "p-value R-squared")
colnames(binary.tests.agreement) <- c("Excess significance, Schwarzer", "Excess significance, Peter",
                                      "Excess significance, Rucker", "Excess significance, Harbord",
                                      "Peter, Schwarzer", "Schwarzer, Rucker", "Schwarzer, Harbord",
                                      "Rucker, Peter", "Harbord, Peter")
#--------------------------------------------------------------------------------------------------------------------#


#Continuous:
meta.cont <- meta.cont %>% ungroup() %>% mutate(thompson.egger = ifelse(thompson.test == egger.test, "agree", "disagree"),
                                                thompson.begg = ifelse(thompson.test == begg.test, "agree", "disagree"),
                                                egger.begg = ifelse(egger.test == begg.test, "agree", "disagree"),
                                                tes.d.egger = ifelse(tes.d.test == egger.test, "agree", "disagree"),
                                                thompson.tes.d = ifelse(thompson.test ==tes.d.test, "agree", "disagree"),
                                                tes.d.begg = ifelse(tes.d.test == begg.test, "agree", "disagree"))
agreement.cont <- meta.cont %>% ungroup() %>%  summarise(tes.d.egger = sum(tes.d.egger == "agree")/n(),
                                                         thompson.tes.d = sum(thompson.tes.d == "agree")/n(),
                                                         tes.d.begg = sum(tes.d.begg == "agree")/n(),
                                                         thompson.egger = sum(thompson.egger == "agree")/n(),
                                                         thompson.begg = sum(thompson.begg == "agree")/n(),
                                                         egger.begg = sum(egger.begg == "agree")/n())
correlation.cont <- meta.cont %>% ungroup() %>% summarise(tes.d.egger = cor(stat.egger, stat.d.tes),
                                                          thompson.tes.d = cor(stat.thompson, stat.d.tes),
                                                          tes.d.begg = cor(stat.d.tes, stat.begg),
                                                          thompson.egger = cor(stat.thompson, stat.egger),
                                                          thompson.begg = cor(stat.thompson, stat.begg),
                                                          egger.begg = cor(stat.egger, stat.begg))
rsquared.cont <- meta.cont %>% ungroup() %>% summarise(tes.d.egger = summary(lm(stat.egger~ stat.d.tes))$r.squared,
                                                       thompson.tes.d = summary(lm(stat.thompson~ stat.d.tes))$r.squared,
                                                       tes.d.begg = summary(lm(stat.d.tes~ stat.begg))$r.squared,
                                                       thompson.egger = summary(lm(stat.thompson~ stat.egger))$r.squared,
                                                       thompson.begg = summary(lm(stat.thompson~ stat.begg))$r.squared,
                                                       egger.begg = summary(lm(stat.egger~ stat.begg))$r.squared)
agreement.sig.cont <- meta.cont %>% ungroup() %>%  summarise(tes.d.egger = sum(tes.d.egger == "agree")/sum(egger.test),
                                                             thompson.tes.d = sum(thompson.tes.d == "agree")/sum(thompsom.test),
                                                             tes.d.begg = sum(tes.d.begg == "agree")/sum(begg.test),
                                                             thompson.egger = sum(thompson.egger == "agree")/sum(egger.test),
                                                             thompson.begg = sum(thompson.begg == "agree")/sum(thompson.test),
                                                             egger.begg = sum(egger.begg == "agree")/sum(egger.test))



cont.tests.agreement <- rbind(agreement.cont, correlation.cont, rsquared.cont)
rownames(cont.tests.agreement) <- c("Test Agreement","p-value Correlation", "p-value R-squared")
colnames(cont.tests.agreement) <- c( "Excess significance, Egger", "Excess significance, Thompson", 
                                     "Excess significance, Begg","Thompson, Egger", "Thompson, Begg", "Egger, Begg")

#--------------------------------------------------------------------------------------------------------------------#


#Survival:
meta.surv <- meta.surv %>% ungroup() %>% mutate(thompson.egger = ifelse(thompson.test == egger.test, "agree", "disagree"),
                                                thompson.begg = ifelse(thompson.test == begg.test, "agree", "disagree"),
                                                egger.begg = ifelse(egger.test == begg.test, "agree", "disagree"),
                                                tes.d.egger = ifelse(tes.d.test == egger.test, "agree", "disagree"),
                                                thompson.tes.d = ifelse(thompson.test ==tes.d.test, "agree", "disagree"),
                                                tes.d.begg = ifelse(tes.d.test == begg.test, "agree", "disagree"))
agreement.surv <- meta.surv %>% ungroup() %>%  summarise(tes.d.egger = sum(tes.d.egger == "agree")/n(),
                                                         thompson.tes.d = sum(thompson.tes.d == "agree")/n(),
                                                         tes.d.begg = sum(tes.d.begg == "agree")/n(),
                                                         thompson.egger = sum(thompson.egger == "agree")/n(),
                                                         thompson.begg = sum(thompson.begg == "agree")/n(),
                                                         egger.begg = sum(egger.begg == "agree")/n())
correlation.surv <- meta.surv %>% ungroup() %>% summarise(tes.d.egger = cor(stat.egger, stat.d.tes),
                                                          thompson.tes.d = cor(stat.thompson, stat.d.tes),
                                                          tes.d.begg = cor(stat.d.tes, stat.begg),
                                                          thompson.egger = cor(stat.thompson, stat.egger),
                                                          thompson.begg = cor(stat.thompson, stat.begg),
                                                          egger.begg = cor(stat.egger, stat.begg))
# rsquared.surv <- meta.surv %>% ungroup() %>% summarise(tes.d.egger = summary(lm(stat.egger~ stat.d.tes))$r.squared,
#                                                           thompson.tes.d = summary(lm(stat.thompson~ stat.d.tes))$r.squared,
#                                                           tes.d.begg = summary(lm(stat.d.tes~ stat.begg))$r.squared,
#                                                           thompson.egger = summary(lm(stat.thompson~ stat.egger))$r.squared,
#                                                           thompson.begg = summary(lm(stat.thompson~ stat.begg))$r.squared,
#                                                           egger.begg = summary(lm(stat.egger~ stat.begg))$r.squared)
meta.surv <- meta.surv %>% ungroup() %>% mutate(thompson.egger = ifelse(thompson.test + egger.test > 1, "agree", "disagree"),
                                                thompson.begg = ifelse(thompson.test + begg.test > 1, "agree", "disagree"),
                                                egger.begg = ifelse(egger.test + begg.test > 1, "agree", "disagree"),
                                                tes.d.egger = ifelse(tes.d.test + egger.test > 1, "agree", "disagree"),
                                                thompson.tes.d = ifelse(thompson.test +tes.d.test > 1, "agree", "disagree"),
                                                tes.d.begg = ifelse(tes.d.test + begg.test > 1, "agree", "disagree"))
agreement.sig.surv <- meta.surv %>% ungroup() %>%  summarise(tes.d.egger = sum(tes.d.egger == "agree")/sum(tes.d.test),
                                                             thompson.tes.d = sum(thompson.tes.d == "agree")/sum(tes.d.test),
                                                             tes.d.begg = sum(tes.d.begg == "agree")/sum(tes.d.test),
                                                             thompson.egger = sum(thompson.egger == "agree")/sum(egger.test),
                                                             thompson.begg = sum(thompson.begg == "agree")/sum(thompson.test),
                                                             egger.begg = sum(egger.begg == "agree")/sum(egger.test))

surv.tests.agreement <- rbind(agreement.surv, correlation.surv, rsquared.surv)
rownames(surv.tests.agreement) <- c("Test Agreement","p-value Correlation", "p-value R-squared")
colnames(surv.tests.agreement) <- c("Excess significance, Egger (survival)", "Excess significance, Thompson (survival)", 
                                    "Excess significance, Begg (survival)", "Thompson, Egger (survival)", 
                                    "Thompson, Begg (survival)", "Egger, Begg (survival)")

#--------------------------------------------------------------------------------------------------------------------#

#Merging:
test.agreement <- rbind(t(binary.tests.agreement), t(cont.tests.agreement), t(surv.tests.agreement))

#--------------------------------------------------------------------------------------------------------------------#

#Number of missing studies, according to copas:
hist(meta.f$missing.copas / meta.f$n)

#--------------------------------------------------------------------------------------------------------------------#

#Correlations btw. teststatistics:
cor.data <- meta.f %>% select(stat.rucker, stat.harbord, stat.peter, stat.schwarzer, stat.d.tes, 
                              stat.egger, stat.thompson, stat.begg)
colnames(cor.data) <- c("Rucker", "Harbord", "Peter", "Schwarzer", "Excess s.",  "Egger", "Thompson",
                        "Begg")

panel.pts.new <- function (x, y, corr = NULL, col.regions, cor.method, ...) 
{
  
  plot.xy(xy.coords(x, y), type = "p", ...)
  box(col = "lightgray")
  rect(xleft = -1.96, ybottom = -1.96, xright = 1.96, ytop = 1.96, density = NULL, angle = 45,
       col = NA, border = "white", lty = 1, lwd = 1.2, 
       ...)
  if (!is.null(corr)) 
    return()
  
}

corrgram(cor.data, 
         upper.panel=panel.pts.new, lower.panel=panel.cor, 
         cex.labels = 1, cex = .5, cex.cor = 2,
         text.panel=panel.txt, cor.method = "spearman")














# #Comparison of m.a. with copas correction method applied:
# sig.zcor <- meta.f %>% filter(!is.na(est.z.copas)) %>% 
# 	mutate(z.fixef = est.z.fixef/se.est.z.fixef,
# 				 z.ranef = est.z.ranef/se.est.z.ranef,
# 				 z.reg = est.z.reg/se.est.z.reg,
# 				 z.copas = est.z.copas/se.est.z.copas) %>%  
# 	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
# 	group_by(method) %>% 
# 	summarise(significant = length(which(abs(fisher.z) > 1.96)),
# 						p.significant = significant/length(fisher.z))
# sig.zcor <- data.frame(method = sig.zcor$method,
# 											 label = paste("n = ", sig.zcor$significant, ", ", round(sig.zcor$p.significant,3)*100, "% significant", sep = ""))
# 
# meta.f %>% filter(!is.na(est.z.copas)) %>% 
# 	mutate(z.fixef = est.z.fixef/se.est.z.fixef,
# 				 z.ranef = est.z.ranef/se.est.z.ranef,
# 				 z.reg = est.z.reg/se.est.z.reg,
# 				 z.copas = est.z.copas/se.est.z.copas) %>%  
# 	select(z.fixef, z.ranef, z.reg, z.copas) %>% gather(key = "method", value = "fisher.z") %>% 
# 	# filter(abs(fisher.z) < 10) %>% 
# 	ggplot(aes(x = abs(fisher.z))) + geom_histogram() + theme_bw() + facet_wrap(~method, ncol = 2) + 
# 	geom_vline(xintercept = 1.96, color = "red") + 
# 	geom_text(data = sig.zcor, aes(x = 10, y = 125, label = label), position = "dodge")
# 

