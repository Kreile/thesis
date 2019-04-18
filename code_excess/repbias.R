PATH_HOME = path.expand("~") # user home
PATH = file.path(PATH_HOME, 'thesis/code')
source(file.path(PATH, 'prepare.R'))

require(tidyverse)
require(biostatUZH)
require(meta)

#Do a Peterts Test for Publication bias for all reproduction trial groups with Odds ratios or Risk Ratios and 
#more than 10 trials per group. 

#Problem child subgroup: Metabias regression test fails because solve.default(t(X) %*% W %*% X) doesnt work.
print(yodata %>% filter(file.nr == 3014 & outcome.nr == 3 & subgroup.nr == 1) %>% select(events1, events2, total1, total2))
#excluded from analysis.

yodata <- yodata %>% filter(file.nr != 3014)

yodata %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% #filter(file.nr < 503) %>% 
  filter(events1 > 0 | events2 > 0) %>% 
  filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
  group_by(file.nr, outcome.nr, subgroup.nr) %>% 
  mutate(n = n()) %>% filter(n > 9) %>% 
  summarize(pval = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "peters")$p.val) %>% 
  ggplot(aes(x = pval)) + geom_histogram(col = "gray15", fill = "dodgerblue", bins = 20) +
  theme_bw() + labs(title = "Peters Reporting Bias Test P-values for Binary Outcome Meta-Analyses") + xlab("P-value") + ylab("Frequency")

#Do a Eggers and Thompson Sharp and ?Begg and Mazumdar test for continuous outcomes.
yodata %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% #filter(file.nr < 503) %>% 
  filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
  filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
  group_by(file.nr, outcome.nr, subgroup.nr) %>% 
  mutate(n = n()) %>% filter(n > 9) %>% 
  summarize(pval = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                            method = "linreg")$p.val) %>% 
  ggplot(aes(x = pval)) + geom_histogram(col = "gray15", fill = "dodgerblue", bins = 20) +
  theme_bw() + labs(title = "Eggers Reporting Bias Test P-values for Continuous Outcome Meta-Analyses") + xlab("P-value") + ylab("Frequency")

yodata %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% #filter(file.nr < 503) %>% 
  filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
  filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
  group_by(file.nr, outcome.nr, subgroup.nr) %>% 
  mutate(n = n()) %>% filter(n > 9) %>% 
  summarize(pval = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                            method = "mm")$p.val) %>% 
  ggplot(aes(x = pval)) + geom_histogram(col = "gray15", fill = "dodgerblue", bins = 20) +
  theme_bw() + labs(title = "Thomson Sharp Reporting Bias Test P-values for Continuous Outcome Meta-Analyses") + xlab("P-value") + ylab("Frequency")

yodata %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% #filter(file.nr < 503) %>% 
  filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
  filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
  group_by(file.nr, outcome.nr, subgroup.nr) %>% 
  mutate(n = n()) %>% filter(n > 9) %>% 
  summarize(pval = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                            method = "rank")$p.val) %>% 
  ggplot(aes(x = pval)) + geom_histogram(col = "gray15", fill = "dodgerblue", bins = 20) +
  theme_bw() + labs(title = "Begg and Mazumdar Reporting Bias Test P-values for Continuous Outcome Meta-Analyses") + xlab("P-value") + ylab("Frequency")


#Trim and fill - fraction of trimmed comparisons of all comparisons in meta-analysis.
#Continuous outcomes
yodata %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% #filter(file.nr < 50) %>% 
  filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
  filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% 
  group_by(file.nr, outcome.nr, subgroup.nr) %>% 
  mutate(n = n()) %>% filter(n > 9) %>% 
  summarize(trim = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2))$k0 / n()) %>% 
  ggplot(aes(x = trim)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Fraction of Trimmed Comparisons") + xlab("Fraction") + ylab("Frequency")

#Binary outcomes
yodata %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% #filter(file.nr < 503) %>% 
  filter(events1 > 0 | events2 > 0) %>% 
  filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
  group_by(file.nr, outcome.nr, subgroup.nr) %>% 
  mutate(n = n()) %>% filter(n > 9) %>% 
  summarize(trim = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2))$k0 / n()) %>% 
  ggplot(aes(x = trim)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Fraction of Trimmed Comparisons") + xlab("Fraction") + ylab("Frequency")




tp <- yodata %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% filter(file.nr < 50) %>% 
  filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% 
  filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2))
  
ms <- metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, data=  tp)
trimfill(ms)$k0






















  
# print(yodata %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% filter(file.nr == 3014) %>% 
#         filter(events1 > 0 | events2 > 0) %>% 
#         filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
#   group_by(file.nr, outcome.nr, subgroup.nr) %>% mutate(n = n()) %>% filter(n > 9) %>% 
#     select(events1, events2, total1, total2), n = 250)
# 
# tp <- yodata %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% filter(file.nr == 3014) %>% 
#   filter(outcome.nr == 1 & subgroup.nr == 2) %>% 
#   filter(events1 > 0 | events2 > 0) %>% 
#   filter(total1 - events1 > 0 | total2 - events2 > 0) %>%
#   group_by(file.nr, outcome.nr, subgroup.nr) %>% mutate(n = n()) %>% filter(n > 9) %>% 
#   select(events1, events2, total1, total2)
# 
# tp.m <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR", data = tp)
# metabias(tp.m, method = "peters")$p.value
# yodata %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% filter(file.nr == 12) %>% 
#   filter(events1 > 0 | events2 > 0) %>% 
#   group_by(file.nr, outcome.nr, subgroup.nr) %>% mutate(n = n()) %>% filter(n > 9) %>% 
#   summarize(pval = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "peters")$p.val)