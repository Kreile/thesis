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
    sig.reg = ifelse(pval.reg > sig.level, 0, 1),
    sig.copas = ifelse(pval.copas > sig.level, 0, 1),
    
    sig.change.ranef.reg = sig.change(sig.before = sig.ranef, sig.after =  sig.reg), #Function to register if significance of estimate has changed after correction
    sig.change.ranef.copas = sig.change(sig.before = sig.ranef, sig.after =  sig.copas),
    sig.change.fixef.reg = sig.change(sig.before = sig.fixef, sig.after =  sig.reg),
    sig.change.fixef.copas = sig.change(sig.before = sig.fixef, sig.after =  sig.copas))


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


meta.bin <- meta.bin %>%ungroup() %>%  mutate(egger.schwarzer = ifelse(egger.test == schwarzer.test, "agree", "disagree"),
                                              egger.peter = ifelse(egger.test == peter.test, "agree", "disagree"),
                                              egger.rucker = ifelse(egger.test == rucker.test, "agree", "disagree"),
                                              egger.harbord = ifelse(egger.test == harbord.test, "agree", "disagree"),
                                              schwarzer.peter = ifelse(schwarzer.test == peter.test, "agree", "disagree"),
                                              schwarzer.rucker = ifelse(schwarzer.test == rucker.test, "agree", "disagree"),
                                              schwarzer.harbord = ifelse(schwarzer.test == harbord.test, "agree", "disagree"),
                                              rucker.peter = ifelse(rucker.test == peter.test, "agree", "disagree"),
                                              rucker.harbord = ifelse(rucker.test == harbord.test, "agree", "disagree"),
                                              harbord.peter = ifelse(harbord.test == peter.test, "agree", "disagree"))

meta.cont <- meta.cont %>% ungroup() %>% mutate(thompson.egger = ifelse(thompson.test == egger.test, "agree", "disagree"),
                                                thompson.begg = ifelse(thompson.test == begg.test, "agree", "disagree"),
                                                egger.begg = ifelse(egger.test == begg.test, "agree", "disagree"))

test <- meta.f %>% group_by(outcome.type) %>% summarize(
  schwarzer.test = mean(schwarzer.test),
  rucker.test = mean(rucker.test),
  harbord.test = mean(harbord.test),
  peter.test = mean(peter.test),
  egger.test = mean(egger.test),
  begg.test = mean(begg.test),
  thompson.test = mean(thompson.test))
test <- test %>% gather(key = "test.type", value = "mean", schwarzer.test:thompson.test) 
# test$mean[which(!is.na(test$mean))] <- round(test$mean[which(!is.na(test$mean))], 2)

#Compare proportions of total concordance vs. proportion of at least single positive test results among all
meta.f %>% ungroup() %>% 
  select(outcome.type, schwarzer.test, rucker.test, harbord.test, peter.test, egger.test, begg.test, thompson.test) %>% 
  gather(key = "test.type", value = "null.hypothesis", schwarzer.test:thompson.test) %>%  
  mutate(null.hypothesis = factor(ifelse(null.hypothesis == 1, "significant", "not significant"))) %>% 
  ggplot(aes(x = test.type, fill = null.hypothesis)) + geom_bar() + coord_flip() + 
  theme_bw() + xlab(label = NULL) + facet_wrap(~outcome.type, scales = "free") + theme(legend.position = "top") +
  guides(fill=guide_legend(title=NULL)) 

# + annotate("text", x = test$test.type, y = 10,
#            label = paste(round(test$mean, 2,)*100, "% rejected"),
#            color = "white")

meta.f %>% 
  select(outcome.type, schwarzer.test, rucker.test, harbord.test, peter.test, egger.test, begg.test, thompson.test) %>% 
  rowwise() %>% 
  mutate(concordance.sum = 
           sum(schwarzer.test, rucker.test, harbord.test, peter.test, egger.test, begg.test, thompson.test, na.rm = TRUE), 
         null.hypothesis = factor(ifelse(concordance.sum > 0, "significant", "not significant"))) %>% 
  ggplot(aes(outcome.type)) + 
  geom_bar(aes(fill = null.hypothesis), position= position_stack(reverse = T)) + coord_flip() + 
  theme_bw() + xlab(label = NULL)  + theme(legend.position = "top") +
  guides(fill=guide_legend(title=NULL)) 

meta.f %>% 
  select(outcome.type, schwarzer.test, rucker.test, harbord.test, peter.test, egger.test, begg.test, thompson.test) %>% 
  rowwise() %>% 
  mutate(concordance.sum = 
           sum(schwarzer.test, rucker.test, harbord.test, peter.test, egger.test, begg.test, thompson.test, na.rm = TRUE), 
         null.hypothesis = (ifelse(concordance.sum > 0, 1, 0))) %>%
  group_by(outcome.type) %>% summarise(null.hypothesis = mean(null.hypothesis, na.rm = T)) %>% 
  ggplot(aes(outcome.type)) + geom_bar(aes(weight = null.hypothesis)) + coord_flip() + 
  theme_bw() + xlab(label = NULL)  + theme(legend.position = "top") +
  guides(fill=guide_legend(title=NULL)) 
