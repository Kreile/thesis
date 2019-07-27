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

load(file.path(PATH_DATA, "PubBias_2019-07-19.RData"))
#----------------------------------------------------------------------------------------------#

##-------Prepare data--------#
reg.dat <- data.ext2 %>% filter(meta.id %in% meta.id.vector) #Filter meta-analysis data

#Mirror all effects on one positive side:
reg.dat <- reg.dat %>% group_by(meta.id) %>% 
  mutate(bias.side = bias.side.fct2(outcome = outcome.flag, outcome.measure.merged = outcome.measure.merged,
                                    lrr = lrr, var.lrr = var.lrr, smd = cohensd, var.smd = var.cohensd, 
                                    effect = effect, se = se),
  smd.pool = smd.pool * bias.side) #Mirror effect sizes on one side

#Get tau^2 based on smd's and introduce new study se:
reg.dat <- reg.dat %>% group_by(meta.id) %>% 
  mutate(tau.squared = metagen(TE = smd.pool, seTE = se.smd.pool)$tau^2,
         new.se.smd = sqrt(se.smd.pool^2 + tau.squared))

reg.dat <- reg.dat %>% filter(se.smd.pool != 0) %>% filter(!is.na(smd.pool)) %>% 
  filter(!is.na(se.smd.pool))

#Get centered study years (centered with respect to meta-analysis study years and all study years (ov))
reg.dat <- reg.dat %>% group_by(meta.id) %>% 
  mutate(study.year.center = scale(scale = F, center = T, x = study.year)) %>% 
  ungroup() %>% mutate(study.year.center.ov = scale(scale = F, center = T, x = study.year)[,1]) 

reg.dat <- merge(reg.dat, data.review[,c("id", "year")], by = "id")

reg.dat <- reg.dat %>% filter(se.smd.pool < 3)
#----------------------------------------------------------------------------------------------#

##-------Modelling--------#
small.study.fit <- lme(fixed = smd.pool ~ se.smd.pool, random = list(~ 1 | id, ~ 1 | meta.id),
                       weights = ~ new.se.smd^2, method = "ML", data = reg.dat)
anova(small.study.fit)

#----------------------------------------------------------------------------------------------#

##-------Illustration--------#
lmer.table <- function(model){
  sum <- summary(model)
  estimate <- round(sum$coefficients$fixed)
  ci <- intervals(model)
  tb <- data.frame(estimate = round(ci$fixed[,2], 3), CI = round(ci$fixed[,c(1,3)], 3))
  colnames(tb)[c(2,3)] <- c("2.5%CL", "97.5%CL")
  return(tb)
}

require(GGally)
#topairsplot.pb.stats.bin <- abs(topairsplot.pb.stats.bin)
notplot.dat <- reg.dat %>%  filter(se.smd.pool >= 2 | smd.pool >= 5)
plot.dat <- reg.dat %>% select(smd.pool, se.smd.pool, new.se.smd, study.year.center, 
                               study.year.center.ov, id, meta.id, outcome.measure.merged) %>% 
  filter(se.smd.pool < 1.53) %>% filter(smd.pool < 5)

# pairs(plot.dat[,c(1,2,5)], lower.panel = NULL, cex = .1)
ggpairs(plot.dat, columns = c(1,2,5), 
        upper = list(continuous = wrap("points", size = 0.15, alpha = .1, )),
        axisLabels = "show")
ggpairs(topairsplot.pb.stats.bin, columns = 1:5, aes(alpha = 0.5))

ggpairs(topairsplot.pb.stats.cont, columns = 1:3, aes(colour = factor(egger.test)))
ggpairs(topairsplot.pb.stats.cont, columns = 1:3, aes(alpha = 0.7))


anova.xtable <- xtable(anova(small.study.fit, small.study.year.ov.fit, small.study.year.ov.int.fit))

small.study.fit.table <- lmer.table(small.study.fit)
small.study.year.ov.table <- lmer.table(small.study.year.ov.fit)
# small.study.year.ov.int.table <- lmer.table(small.study.year.ov.int.fit)




# print(xtable(kreatinin.tba, align = "lcccc"), floating = F, size = "scriptsize")





year.fit <- lme(fixed = smd.pool ~ study.year.center, random = list(~ 1 | id, ~ 1 | meta.id),
                weights = ~ new.se.smd^2, method = "ML", data = reg.dat)




summary(small.study.year.fit)
anova(small.study.fit, small.study.year.fit) 



summary(small.study.year.ov.int.fit)

reg.dat$prediction <- predict(object = small.study.year.ov.fit, newdata = reg.dat)

plot(y = reg.dat$prediction, reg.dat$study.year.center.ov, cex = .1)
plot(y = reg.dat$smd.pool, reg.dat$study.year.center.ov, cex = .1)
abline(small.study.year.ov.fit)

plot(y = reg.dat$smd.pool, x = reg.dat$se.smd.pool, cex = .1)
small.study.fit.coef <- summary(small.study.fit)$coefficients$fixed
plot(small.study.fit)
abline(a = small.study.fit.coef, col = "red", lty = 2)
plot(small.study.year.ov.fit$fitted[,1], x=  reg.dat$study.year.center.ov, cex = .1)




time.fit <- glmer(formula = smd.pool ~ study.year.center + (1 | meta.id) + (1 | id), weights = 1/se.smd.pool,
                  family = gaussian, data = reg.dat, nAGQ = 10)
time.small.study.fit <- glmer(formula = smd.pool ~ se.smd.pool + study.year.center + (1 | meta.id) + (1 | id), weights = 1/se.smd.pool,
                              family = gaussian, data = reg.dat, nAGQ = 10)
summary(time.fit)
anova(small.study.fit, time.small.study.fit)


time.interaction.fit <- glmer(formula = smd.pool ~  study.year.center + (1 | meta.id) + (1 | id), weights = 1/se.smd.pool,
                              family = gaussian, data = reg.dat, nAGQ = 10)


#----------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------#

##-------Plots:-------#
# par(mfrow = c(1,1))
# plot(y = reg.dat$smd.pool, x = reg.dat$se.smd.pool, cex = .1)
# plot(y = reg.dat$smd.pool, x = reg.dat$study.year.center, cex = .1)

#Check how largest m.a. looks like..
# largest.id <- meta.f$meta.id[2957]
# reg.dat <- reg.dat %>%  filter(meta.id == largest.id)
# plot(y = reg.dat$smd.pool, x = reg.dat$se.smd.pool, cex = .1, xlim = c(0, 1.3), ylim = c(-2.4, 6))
# plot(y = reg.dat$smd.pool, x = reg.dat$study.year.center, cex = .1, xlim = c(-45, 38), ylim = c(-2.4, 6))


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
# lm2.fit <- lm(smd.pool ~ new.se.smd, data = ex.data, weights = 1/(new.se.smd)^2)
# summary(lm2.fit)
# metabias(metagen(TE = smd.pool, seTE = se.smd.pool, data = ex.data), method.bias = "mm")
# # 
# tau.2 <- metagen(TE = smd.pool, seTE = se.smd.pool, data = ex.data, method.tau = "DL")$tau^2
# ex.data$new.se.smd <- sqrt(ex.data$se.smd.pool^2 + tau.2)
# ex.data$new.var.smd.pool <- ex.data$se.smd.pool^2 + tau.2
# lm2.fit <- lm(smd.pool ~ se.smd.pool, data = ex.data, weights = 1/new.se.smd^2)
# summary(lm2.fit)
# metabias(metagen(TE = smd.pool, seTE = se.smd.pool, data = ex.data, method.tau = "DL"), method.bias = "mm")
# # lm2.fit <- lm(smd.pool ~ se.smd.pool, data = ex.data, weights = 1/se.smd.pool^2)
# # summary(lm2.fit)
# # metabias(metagen(TE = smd.pool, seTE = se.smd.pool, data = ex.data, method.tau = "DL"), method.bias = "linreg")
# coef(lm2.fit)[1] + coef(lm2.fit)[2]*sqrt(tau.2)
# limitmeta(x = metagen(TE = smd.pool, seTE = se.smd.pool, data = ex.data))
# # lm.fit <- lm(formula = smd.pool ~ se.smd.pool, weights = 1/se.smd.pool, data = reg.dat)
# # plot(y = reg.dat$smd.pool, x = reg.dat$se.smd.pool, cex = .1)
# # abline(lm.fit, col = 2)

# # Compare lme() fit to lm() fit:

lm2.fit <- lm(smd.pool ~ se.smd.pool, data = ex.data, weights = 1/(new.se.smd^2))
small.study.fit <- lme(fixed = smd.pool ~ se.smd.pool, random = list(~ 1 | id, ~ 1 | meta.id),
                       weights = ~ new.se.smd^2,
                       data = ex.data, method = "ML")

plot(y = ex.data$smd.pool, x = ex.data$se.smd.pool)
abline(small.study.fit)
abline(lm2.fit)
# 
# coef(lm2.fit)
# coef(small.study.fit)
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
  


# small.study.year.ov.fit <- lme(fixed = smd.pool ~ se.smd.pool + study.year.center.ov, random = list(~ 1 | id, ~ 1 | meta.id),
#                                weights = ~ new.se.smd^2, method = "ML", data = reg.dat)
# small.study.year.ov.int.fit <- lme(fixed = smd.pool ~ se.smd.pool * study.year.center.ov, random = list(~ 1 | id, ~ 1 | meta.id),
#                                    weights = ~ new.se.smd^2, method = "ML", data = reg.dat)
# 
# anova(small.study.fit, small.study.year.ov.fit, small.study.year.ov.int.fit) 

