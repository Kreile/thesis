rm(list = ls())
require(biostatUZH)
require(tidyverse)
require(meta)
require(psych)

##########################################################
# read the generated (single) CSV file into a data frame:
# PATH_HOME = path.expand("~") # user home
# PATH = file.path(PATH_HOME, 'thesis/cochrane')
# source(file.path(PATH, 'prepare.R'))
# madata <- read.csv("~/thesis/cochrane/cochrane_2018-06-09.csv",
#                    colClasses=c("integer", "character", "character",  # file       :  nr, name, doi
#                                 "integer", "integer",                 # file       :  index, version
#                                 "integer", "character",               # comparison :  nr, name
#                                 "integer", "character", "character",  # outcome    :  nr, name, measure
#                                 "integer", "character",               # subgroup   :  nr, name
#                                 "character", "integer",               # study      :  name, year
#                                 "numeric", "numeric",                 # effect, std.err.
#                                 "integer", "numeric",                 # group 1    :  events, total
#                                 "numeric", "numeric",                 # group 1    :  mean, se
#                                 "integer", "numeric",                 # group 2    :  events, total
#                                 "numeric", "numeric",                 # group 2    :  mean, se
#                                 "integer"))                           # total N
# madata$total1 <- as.integer(madata$total1)
# madata$total2 <- as.integer(madata$total2)
##########################################################
# read the generated (single) CSV file into a data frame:

yodata <- madata %>% filter(!is.na(outcome.measure))

#Convert identical but differently spelled outcome names (~ 50).

yodata[yodata$outcome.measure == "RR", ] <- yodata %>% filter(outcome.measure == "RR") %>% mutate(outcome.measure = "Risk Ratio")
yodata[yodata$outcome.measure == "Risk ratio", ] <- yodata %>% filter(outcome.measure == "Risk ratio") %>% mutate(outcome.measure = "Risk Ratio")
#yodata[yodata$outcome.measure == "Incidence Risk Ratio", ] <- yodata %>% filter(outcome.measure == "Incidence Risk Ratio") %>% mutate(outcome.measure = "Risk Ratio")
#yodata[yodata$outcome.measure == "Relatvie risks", ] <- yodata %>% filter(outcome.measure == "Relatvie risks") %>% mutate(outcome.measure = "Risk Ratio")

yodata[yodata$outcome.measure == "Mean difference", ] <- yodata %>% filter(outcome.measure == "Mean difference") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "mean difference", ] <- yodata %>% filter(outcome.measure == "mean difference") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "mean diff (L/min)", ] <- yodata %>% filter(outcome.measure == "mean diff (L/min)") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "mean difference (L)", ] <- yodata %>% filter(outcome.measure == "mean difference (L)") %>% mutate(outcome.measure = "Mean Difference")
#yodata[yodata$outcome.measure == "Mean Change Difference", ] <- yodata %>% filter(outcome.measure == "Mean Change Difference") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "MD", ] <- yodata %>% filter(outcome.measure == "MD") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "difference in means", ] <- yodata %>% filter(outcome.measure == "difference in means") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "MD in SBP", ] <- yodata %>% filter(outcome.measure == "MD in SBP") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "MD in DBP", ] <- yodata %>% filter(outcome.measure == "MD in DBP") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "Diff in mean score", ] <- yodata %>% filter(outcome.measure == "Diff in mean score") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "Mean Difference (SDs)", ] <- yodata %>% filter(outcome.measure == "Mean Difference (SDs)") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "Mean diff in DBP", ] <- yodata %>% filter(outcome.measure == "Mean diff in DBP") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "Mean diff in SBP", ] <- yodata %>% filter(outcome.measure == "Mean diff in SBP") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "Mean Difference [%]", ] <- yodata %>% filter(outcome.measure == "Mean Difference [%]") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "MD or Difference-in-Differences (SDs)", ] <- yodata %>% filter(outcome.measure == "MD or Difference-in-Differences (SDs)") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "Diff in mean score", ] <- yodata %>% filter(outcome.measure == "Diff in mean score") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "MD in serum Ca", ] <- yodata %>% filter(outcome.measure == "MD in serum Ca") %>% mutate(outcome.measure = "Mean Difference")
yodata[yodata$outcome.measure == "Mean % change", ] <- yodata %>% filter(outcome.measure == "Mean % change") %>% mutate(outcome.measure = "Mean Difference")
# yodata[yodata$outcome.measure == "% change from baseline", ] <- yodata %>% filter(outcome.measure == "% change from baseline") %>% mutate(outcome.measure = "Mean Difference")

yodata[yodata$outcome.measure == "SMD", ] <- yodata %>% filter(outcome.measure == "SMD") %>% mutate(outcome.measure = "Std. Mean Difference")

yodata[yodata$outcome.measure == "OR", ] <- yodata %>% filter(outcome.measure == "OR") %>% mutate(outcome.measure = "Odds Ratio")
#yodata[yodata$outcome.measure == "Odds Ratio (Non-event)", ] <- yodata %>% filter(outcome.measure == "Odds Ratio (Non-event)") %>% mutate(outcome.measure = "Odds Ratio")
yodata[yodata$outcome.measure == "odds ratio", ] <- yodata %>% filter(outcome.measure == "odds ratio") %>% mutate(outcome.measure = "Odds Ratio")
yodata[yodata$outcome.measure == "Odds Ratios", ] <- yodata %>% filter(outcome.measure == "Odds Ratios") %>% mutate(outcome.measure = "Odds Ratio")
yodata[yodata$outcome.measure == "Odds ratio", ] <- yodata %>% filter(outcome.measure == "Odds ratio") %>% mutate(outcome.measure = "Odds Ratio")
yodata[yodata$outcome.measure == "odds ratios", ] <- yodata %>% filter(outcome.measure == "odds ratios") %>% mutate(outcome.measure = "Odds Ratio")

#yodata[yodata$outcome.measure == "Peto Odds Ratio (Non-event)", ] <- yodata %>% filter(outcome.measure == "Peto Odds Ratio (Non-event)") %>% mutate(outcome.measure = "Peto Odds Ratio")

yodata[yodata$outcome.measure == "Hazard ratio", ] <- yodata %>% filter(outcome.measure == "Hazard ratio") %>% mutate(outcome.measure = "Hazard Ratio")
yodata[yodata$outcome.measure == "Hazards ratio", ] <- yodata %>% filter(outcome.measure == "Hazards ratio") %>% mutate(outcome.measure = "Hazard Ratio")
yodata[yodata$outcome.measure == "hazards ratio", ] <- yodata %>% filter(outcome.measure == "hazards ratio") %>% mutate(outcome.measure = "Hazard Ratio")
yodata[yodata$outcome.measure == "hazard ratio", ] <- yodata %>% filter(outcome.measure == "hazard ratio") %>% mutate(outcome.measure = "Hazard Ratio")
yodata[yodata$outcome.measure == "Survival HR", ] <- yodata %>% filter(outcome.measure == "Survival HR") %>% mutate(outcome.measure = "Hazard Ratio")
#yodata[yodata$outcome.measure == "RRorHR", ] <- yodata %>% filter(outcome.measure == "RRorHR") %>% mutate(outcome.measure = "Hazard Ratio")
#yodata[yodata$outcome.measure == "HR and variance", ] <- yodata %>% filter(outcome.measure == "HR and variance") %>% mutate(outcome.measure = "Hazard Ratio")

yodata[yodata$outcome.measure == "Incidence Rate Ratio", ] <- yodata %>% filter(outcome.measure == "Incidence Rate Ratio") %>% mutate(outcome.measure = "Rate Ratio")
yodata[yodata$outcome.measure == "Rate ratio", ] <- yodata %>% filter(outcome.measure == "Rate ratio") %>% mutate(outcome.measure = "Rate Ratio")
yodata[yodata$outcome.measure == "rate ratio", ] <- yodata %>% filter(outcome.measure == "rate ratio") %>% mutate(outcome.measure = "Rate Ratio")
yodata[yodata$outcome.measure == "Incidence rate ratio", ] <- yodata %>% filter(outcome.measure == "Incidence rate ratio") %>% mutate(outcome.measure = "Rate Ratio")
yodata[yodata$outcome.measure == "incidence rate ratio", ] <- yodata %>% filter(outcome.measure == "incidence rate ratio") %>% mutate(outcome.measure = "Rate Ratio")
#yodata[yodata$outcome.measure == "Relative rate", ] <- yodata %>% filter(outcome.measure == "Relative rate") %>% mutate(outcome.measure = "Rate Ratio")


#yodata[yodata$outcome.measure == "Annualized risk difference (%)", ] <- yodata %>% filter(outcome.measure == "Annualized risk difference (%)") %>% mutate(outcome.measure = "Risk Difference")
yodata[yodata$outcome.measure == "Risk difference (RD)", ] <- yodata %>% filter(outcome.measure == "Risk difference (RD)") %>% mutate(outcome.measure = "Risk Difference")

yodata[yodata$outcome.measure == "rate difference", ] <- yodata %>% filter(outcome.measure == "rate difference") %>% mutate(outcome.measure = "Rate difference")
yodata[yodata$outcome.measure == "Prevented fraction", ] <- yodata %>% filter(outcome.measure == "Prevented fraction") %>% mutate(outcome.measure = "Prevented Fraction")
yodata[yodata$outcome.measure == "Hedges &#180;g", ] <- yodata %>% filter(outcome.measure == "Hedges &#180;g") %>% mutate(outcome.measure = "Hedges` g")
yodata[yodata$outcome.measure == "Hedges&#180; g", ] <- yodata %>% filter(outcome.measure == "Hedges&#180; g") %>% mutate(outcome.measure = "Hedges` g")
yodata[yodata$outcome.measure == "Hedges ` g", ] <- yodata %>% filter(outcome.measure == "Hedges ` g") %>% mutate(outcome.measure = "Hedges` g")
yodata[yodata$outcome.measure == "Hedges`g", ] <- yodata %>% filter(outcome.measure == "Hedges`g") %>% mutate(outcome.measure = "Hedges` g")

#yodata %>% group_by(outcome.measure) %>% count() %>% arrange(desc(n))



#Transform Binary outcomes to odds ratios -> standardized mean difference -> correlation -> fisher scala correlation
# yodata <- yodata %>% 
#   filter(events1 > 0 & events2 > 0) %>% 
#   mutate(odds.ratio = (events1/(total1 - events1)) / (events2/(total2 - events2))) %>% 
#   mutate(std.mean.d = log(odds.ratio) * (sqrt(3)/pi)) %>% 
#   mutate(correlation = std.mean.d/sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2*total1))) %>% 
#   mutate(fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), variance = 1/(total1 + total2 - 3))
# 
# #Transform mean difference to standardized mean difference -> correlation -> fisher scala correlation
# yodata[yodata$outcome.measure == "Mean Difference", ] <- yodata %>% filter(outcome.measure == "Mean Difference") %>% 
#   mutate(std.mean.d = effect/se) %>% 
#   mutate(correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1))) %>% 
#   mutate(fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), variance = 1/(total1 + total2 - 3))
# 
# yodata[yodata$outcome.measure == "Std. Mean Difference", ] <- yodata %>% filter(outcome.measure == "Std. Mean Difference") %>% 
#   mutate(std.mean.d = effect) %>% 
#   mutate(correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1))) %>% 
#   mutate(fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), variance = 1/(total1 + total2 - 3))
# 
# yodata[yodata$outcome.measure == "Mean difference", ] <- yodata %>% filter(outcome.measure == "Mean difference") %>% 
#   mutate(std.mean.d = effect) %>% 
#   mutate(correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1))) %>% 
#   mutate(fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), variance = 1/(total1 + total2 - 3))
# 
# yodata[yodata$outcome.measure == "mean difference", ] <- yodata %>% filter(outcome.measure == "mean difference") %>% 
#   mutate(std.mean.d = effect) %>% 
#   mutate(correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1))) %>% 
#   mutate(fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), variance = 1/(total1 + total2 - 3))

#Filter out fishersz abs(correlation) > 1 and subgroups with fewer than 2 reproductions
# yodata <- yodata %>% filter(fishersz < 1 & fishersz > -1) %>% filter(!is.na(fishersz))
# yodata <- yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
#   filter(fishersz < 1 & fishersz > -1) %>% 
#   mutate(counts = n()) %>% filter(counts > 1)

save(yodata, file = "yodata.RData")







dd <- data_frame(x = rnorm(10), g = c(rep(c("a", "b", "c"), 3), "a"))
dd <- dd %>% group_by(g) %>% mutate(var = var(x))
dd
dd %>% group_by(g) %>% mutate(id  = row_number(), s = sum(x/var) - x/var)
dd %>% group_by(g) %>% mutate(id  = row_number(), s = sum(x) - x[id])


loo.mean <- function(x){
    result <- c()
    for ( i in seq_along(x) ) {
      result[i] <- mean(x[-i], na.rm = TRUE) 
    }
  return(result)
}

dd %>% mutate(l.m = loo.mean(x))
dd %>% group_by(g) %>% mutate(l.m = loo.mean(x))


loo.mean <- function(df){
  n <- count(df)
  return(n)
}

dd %>% loo.mean()
dd %>% group_by(g) %>% loo.mean()

#dd <- data_frame(x = rnorm(10), y = sample(c(1:10), size = 10, replace = T), g = c(rep(c("a", "b", "c"), 3), "a"))
dd %>% group_by(g) %>% mutate( id  = row_number())
dd %>% group_by(g) %>% arrange(y) %>% mutate(id  = row_number())
dd %>% group_by(g) %>% summarise(max(x[-1]))
dd %>% group_by(g) %>% top_n(1)




dd <- dd %>% mutate(jo = "nai")

dd[dd$g == "b",] <- dd %>% filter(g == "b") %>% mutate(jo = "jo")























