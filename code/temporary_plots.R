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





#Comparison of p-values between meta-analysis and adjustment methods:

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









# ################################################################################################
# ################################################################################################
# #PROPORTION OF SIGNIFICANT PUBLICATION BIAS TESTS
# ################################################################################################
# ################################################################################################
# 
# sig.level <- 0.1
# meta.f <- meta.f %>% 
#   mutate(egger.test = ifelse(pval.egger < sig.level, 1, 0),
#          thompson.test = ifelse(pval.thompson < sig.level, 1, 0),
#          begg.test = ifelse(pval.begg < sig.level, 1, 0),
#          
#          schwarzer.test = ifelse(pval.schwarzer < sig.level, 1, 0),
#          rucker.test = ifelse(pval.rucker < sig.level, 1, 0),
#          harbord.test = ifelse(pval.harbord < sig.level, 1, 0),
#          peter.test = ifelse(pval.peter < sig.level, 1, 0))
# 
# meta.bin <- meta.f %>% filter(outcome.type == "bin")
# meta.cont <- meta.f %>% filter(outcome.type == "cont")
# meta.surv <- meta.f %>% filter(outcome.type == "surv")




number.outside.smd <- meta.f %>% mutate(copas = est.d.copas, regression = est.d.reg) %>%
  select(meta.id, est.d.fixef, est.d.ranef, copas, regression) %>%
  gather(key = "method", value = "smd", est.d.fixef:regression) %>% filter(smd > 4)

meta.f %>%   mutate(copas = est.d.copas, regression = est.d.reg) %>%
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

# number.outside.surv <- meta.f %>% mutate(copas = est.z.copas, regression = est.z.reg) %>%
#   select(meta.id, est.z.fixef, est.z.ranef, copas, regression) %>%
#   gather(key = "method", value = "smd", est.z.fixef:regression) %>% filter(smd > 4)

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
################################################################################################
################################################################################################
################################################################################################


meta.f %>% 
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
  xlab(expression(paste("Difference between fixed effects and adjusted ", z, "-score")))


meta.f %>% 
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
  xlab(expression(paste("Difference between random effects and adjusted ", z, "-score")))



meta.f %>% 
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
  xlab(expression(paste("Difference between fixed effects and adjusted Hedges ", g)))
  
  
meta.f  %>% 
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
  xlab(expression(paste("Difference between random effects and adjusted Hedges ", g)))




meta.surv %>% 
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
  facet_wrap(~method) + theme_bw() + 
  xlab("Difference between fixed effects and adjusted log hazard ratio")

meta.surv %>% 
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
  facet_wrap(~method) + theme_bw() + 
  xlab("Difference between random effects and adjusted log hazard ratio")









