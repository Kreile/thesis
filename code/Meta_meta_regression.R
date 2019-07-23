#Load data:
rm(list = ls())
PATH_HOME = path.expand("~") # user home
PATH = file.path(PATH_HOME, 'Data/PubBias')
PATH2 = file.path(PATH_HOME, 'PubBias')
FILE = 'cochrane_2019-07-04.csv'
PATH_DATA = file.path(PATH, 'data')
PATH_CODE = file.path(PATH2, 'code')
PATH_RESULTS = file.path(PATH2, 'results_new')
PATH_FIGURES = file.path(PATH_RESULTS, 'figures')

require(nlme)

source(file.path(PATH_CODE, 'PubBias_functions.R'))

load(file.path(PATH_RESULTS, "data2.RData"))
# data.ext2 <- pb.process3(data)
load(file.path(PATH_RESULTS, "data_used_for_analysis.RData"))

load(file.path(PATH_RESULTS, "meta_analyses_summary_complete.RData"))

load(file.path(PATH_RESULTS, "meta_id_vector.RData"))
#----------------------------------------------------------------------------------------------#

##-------Prepare data--------#
reg.dat <- data.ext2 %>% filter(meta.id %in% meta.id.vector) #Filter meta-analysis data

reg.dat <- reg.dat %>%  mutate(smd = case_when(outcome.flag == "DICH" ~ smd.ordl,
                                         outcome.flag == "CONT" ~ smd),
                         var.smd = case_when(outcome.flag == "DICH" ~ var.smd.ordl,
                                             outcome.flag == "CONT" ~ var.smd),
                         smd = case_when(outcome.measure.new == "Std. Mean Difference" ~ effect, 
                                         TRUE ~ smd),
                         var.smd = case_when(outcome.measure.new == "Std. Mean Difference" ~ se^2,
                                             TRUE ~ var.smd)) #put smd's all together

#Mirror all effects on one positive side:
reg.dat <- reg.dat %>% group_by(meta.id) %>% 
  mutate(bias.side = bias.side.fct2(outcome = outcome.flag, lrr = lrr, var.lrr = var.lrr, 
                                    smd = smd, var.smd = var.smd, effect = effect, se = se),
  smd = smd * bias.side, #Mirror effect sizes on one side
  se.smd = sqrt(var.smd)) 

#Get tau^2 based on smd's and introduce new study se:
reg.dat <- reg.dat %>% group_by(meta.id) %>% 
  mutate(tau.squared = metagen(TE = smd, seTE = se.smd)$tau^2,
         new.se.smd = sqrt(se.smd^2 + tau.squared))

reg.dat <- reg.dat %>% filter(se.smd != 0) %>% filter(!is.na(smd)) %>% 
  filter(!is.na(se.smd)) %>% filter(!is.na(study.year))

#Get centered study years (centered with respect to meta-analysis study years and all study years (ov))
reg.dat <- reg.dat %>% group_by(meta.id) %>% 
  mutate(study.year.center = scale(scale = F, center = T, x = study.year)) %>% 
  ungroup() %>% mutate(study.year.center.ov = scale(scale = T, center = T, x = study.year))

#----------------------------------------------------------------------------------------------#

##-------Plots:-------#
# par(mfrow = c(1,1))
# plot(y = reg.dat$smd, x = reg.dat$se.smd, cex = .1)
# plot(y = reg.dat$smd, x = reg.dat$study.year.center, cex = .1)

#Check how largest m.a. looks like..
# largest.id <- meta.f$meta.id[2957]
# reg.dat <- reg.dat %>%  filter(meta.id == largest.id)
# plot(y = reg.dat$smd, x = reg.dat$se.smd, cex = .1, xlim = c(0, 1.3), ylim = c(-2.4, 6))
# plot(y = reg.dat$smd, x = reg.dat$study.year.center, cex = .1, xlim = c(-45, 38), ylim = c(-2.4, 6))


#----------------------------------------------------------------------------------------------#

##-------Checks:-------#

# # Compare lm() to metabias() and limitmeta():

# reg.dat <- reg.dat %>% filter(meta.id %in% meta.id.vector[c(1:2)])
# print(round(c(meta.f$est.reg, meta.f$est.d.reg, meta.f$est.se, meta.f$est.d.se)),1)
# Check if it works:
ex.id <- 79
ex.data <- reg.dat %>% filter(meta.id == meta.f$meta.id[ex.id])
# #funnel.id(meta.f$meta.id[ex.id])
# meta.f[ex.id, c("pval.rucker", "est.d.fixef", "est.d.reg", "tau")]
# lm2.fit <- lm(smd ~ new.se.smd, data = ex.data, weights = 1/(new.se.smd)^2)
# summary(lm2.fit)
# metabias(metagen(TE = smd, seTE = se.smd, data = ex.data), method.bias = "mm")
# # 
# tau.2 <- metagen(TE = smd, seTE = se.smd, data = ex.data, method.tau = "DL")$tau^2
# ex.data$new.se.smd <- sqrt(ex.data$se.smd^2 + tau.2)
# ex.data$new.var.smd <- ex.data$se.smd^2 + tau.2
# lm2.fit <- lm(smd ~ se.smd, data = ex.data, weights = 1/new.se.smd^2)
# summary(lm2.fit)
# metabias(metagen(TE = smd, seTE = se.smd, data = ex.data, method.tau = "DL"), method.bias = "mm")
# # lm2.fit <- lm(smd ~ se.smd, data = ex.data, weights = 1/se.smd^2)
# # summary(lm2.fit)
# # metabias(metagen(TE = smd, seTE = se.smd, data = ex.data, method.tau = "DL"), method.bias = "linreg")
# coef(lm2.fit)[1] + coef(lm2.fit)[2]*sqrt(tau.2)
# limitmeta(x = metagen(TE = smd, seTE = se.smd, data = ex.data))
# # lm.fit <- lm(formula = smd ~ se.smd, weights = 1/se.smd, data = reg.dat)
# # plot(y = reg.dat$smd, x = reg.dat$se.smd, cex = .1)
# # abline(lm.fit, col = 2)

# # Compare lme() fit to lm() fit:

lm2.fit <- lm(smd ~ se.smd, data = ex.data, weights = 1/(new.se.smd^2))
small.study.fit <- lme(fixed = smd ~ se.smd, random = list(~ 1 | id, ~ 1 | meta.id),
                       weights = ~ new.se.smd^2,
                       data = ex.data, method = "ML")

plot(y = ex.data$smd, x = ex.data$se.smd)
abline(small.study.fit)
abline(lm2.fit)
# 
# coef(lm2.fit)
# coef(small.study.fit)
#----------------------------------------------------------------------------------------------#

#Modelling
small.study.fit <- lme(fixed = smd ~ se.smd, random = list(~ 1 | id, ~ 1 | meta.id),
                       weights = ~ new.se.smd^2, method = "ML", data = reg.dat)

year.fit <- lme(fixed = smd ~ study.year.center, random = list(~ 1 | id, ~ 1 | meta.id),
                weights = ~ new.se.smd^2, method = "ML", data = reg.dat)

small.study.year.ov.fit <- lme(fixed = smd ~ se.smd + study.year.center.ov, random = list(~ 1 | id, ~ 1 | meta.id),
                weights = ~ new.se.smd^2, method = "ML", data = reg.dat)
summary(small.study.year.fit)
anova(small.study.fit, small.study.year.fit) 
small.study.year.ov.int.fit <- lme(fixed = smd ~ se.smd * study.year.center.ov, random = list(~ 1 | id, ~ 1 | meta.id),
                                   weights = ~ new.se.smd^2, method = "ML", data = reg.dat)

anova(small.study.year.ov.fit, small.study.year.ov.int.fit) 
summary(small.study.year.ov.int.fit)

reg.dat$prediction <- predict(object = small.study.year.ov.fit, newdata = reg.dat)

plot(y = reg.dat$prediction, reg.dat$study.year.center.ov, cex = .1)
plot(y = reg.dat$smd, reg.dat$study.year.center.ov, cex = .1)
abline(small.study.year.ov.fit)

plot(y = reg.dat$smd, x = reg.dat$se.smd, cex = .1)
small.study.fit.coef <- summary(small.study.fit)$coefficients$fixed
plot(small.study.fit)
abline(a = small.study.fit.coef, col = "red", lty = 2)
plot(small.study.year.ov.fit$fitted[,1], x=  reg.dat$study.year.center.ov, cex = .1)




time.fit <- glmer(formula = smd ~ study.year.center + (1 | meta.id) + (1 | id), weights = 1/se.smd,
                         family = gaussian, data = reg.dat, nAGQ = 10)
time.small.study.fit <- glmer(formula = smd ~ se.smd + study.year.center + (1 | meta.id) + (1 | id), weights = 1/se.smd,
                  family = gaussian, data = reg.dat, nAGQ = 10)
summary(time.fit)
anova(small.study.fit, time.small.study.fit)


time.interaction.fit <- glmer(formula = smd ~  study.year.center + (1 | meta.id) + (1 | id), weights = 1/se.smd,
                                             family = gaussian, data = reg.dat, nAGQ = 10)


#----------------------------------------------------------------------------------------------#

set.seed(894)
data.ext2 %>% group_by(meta.id) %>% 
  filter(!is.na(study.year)) %>% 
  filter(!is.na(pval.single)) %>% 
  filter(n > 2 ) %>% 
  mutate(timerank = rank(x = study.year, na.last = NA, ties.method = "random")) %>% 
  filter(timerank < 10) %>% 
  filter(abs(effect/se) < 20) %>% 
  ggplot(aes(x = timerank, y = abs(effect/se))) + geom_jitter(size = .1, alpha = .1, height = 0) +
  geom_smooth(method = "lm")

data.ext2 %>% group_by(meta.id) %>% filter(n > 9) %>% 
  mutate(participants = total1 + total2,
         sizerank = rank(x = participants, na.last = NA, ties.method = "random"),
         ) %>% filter(sizerank < 10) %>% 
  ggplot(aes(x = sizerank, y = log(pval.single))) + geom_jitter(height = 0)
  

