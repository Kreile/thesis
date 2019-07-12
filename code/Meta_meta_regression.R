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

require(lme4)

source(file.path(PATH_CODE, 'PubBias_functions.R'))

load(file.path(PATH_RESULTS, "data2.RData"))
# data.ext2 <- pb.process3(data)
load(file.path(PATH_RESULTS, "data_used_for_analysis.RData"))

load(file.path(PATH_RESULTS, "meta_id_vector.RData"))

reg.dat <- data.ext2 %>% filter(meta.id %in% meta.id.vector) #Filter meta-analysis data

reg.dat %>% ungroup() %>% transmute(smd = case_when(outcome.flag == "DICH" ~ smd.ordl,
                                                 outcome.flag == "CONT" ~ smd),
                                    var.smd = case_when(outcome.flag == "DICH" ~ var.smd.ordl,
                                                        outcome.flag == "CONT" ~ var.smd),
                                    smd = case_when(outcome.measure.new == "Std. Mean Difference" ~ effect, 
                                                    TRUE ~ smd),
                                    var.smd = case_when(outcome.measure.new == "Std. Mean Difference" ~ se^2,
                                                        TRUE ~ smd))

reg.dat <- reg.dat %>%  mutate(smd = case_when(outcome.flag == "DICH" ~ smd.ordl,
                                         outcome.flag == "CONT" ~ smd),
                         var.smd = case_when(outcome.flag == "DICH" ~ var.smd.ordl,
                                             outcome.flag == "CONT" ~ var.smd),
                         smd = case_when(outcome.measure.new == "Std. Mean Difference" ~ effect, 
                                         TRUE ~ smd),
                         var.smd = case_when(outcome.measure.new == "Std. Mean Difference" ~ se^2,
                                             TRUE ~ var.smd)) #put smd's all together


reg.dat <- reg.dat %>% group_by(meta.id) %>% 
  mutate(bias.side = bias.side.fct2(outcome = outcome.flag, lrr = lrr, var.lrr = var.lrr, 
                                    smd = smd, var.smd = var.smd, effect = effect, se = se),
  smd = smd * bias.side, #Mirror effect sizes on one side
  se.smd = sqrt(var.smd)) 



reg.dat <- reg.dat %>% filter(se.smd != 0) %>% filter(!is.na(smd)) %>% 
  filter(!is.na(se.smd))

reg.dat <- reg.dat %>% mutate(study.year.center = scale(scale = F, center = T, x = study.year),
                              study.year.scale = scale(scale = T, center = T, x = study.year))

plot(y = reg.dat$smd, x = reg.dat$se.smd)

reg.dat <- reg.dat %>% filter(meta.id %in% meta.id.vector[c(1:2)])
print(round(c(meta.f$est.reg, meta.f$est.d.reg, meta.f$est.se, meta.f$est.d.se)),1)


lm.fit <- lm(formula = smd ~ se.smd, weights = 1/se.smd, data = reg.dat)

coef(lm.fit)[1] + coef(lm.fit)[2]



small.study.fit <- glmer(formula = smd ~ se.smd + (1 | meta.id) + (1 | id), weights = 1/se.smd,
                       family = gaussian, data = reg.dat, nAGQ = 10)

coef(small.study.fit)[1] + coef(small.study.fit)[2]





time.fit <- glmer(formula = smd ~ study.year.center + (1 | meta.id) + (1 | id), weights = 1/se.smd,
                         family = gaussian, data = reg.dat, nAGQ = 10)
time.small.study.fit <- glmer(formula = smd ~ se.smd + study.year.center + (1 | meta.id) + (1 | id), weights = 1/se.smd,
                  family = gaussian, data = reg.dat, nAGQ = 10)

anova(small.study.fit, time.small.study.fit)


time.interaction.fit <- glmer(formula = smd ~  study.year.center + (1 | meta.id) + (1 | id), weights = 1/se.smd,
                                             family = gaussian, data = reg.dat, nAGQ = 10)





