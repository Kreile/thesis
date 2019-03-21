##################################################################################################################################
#Scaled OR over time with initial OR > 1
madata %>% filter(outcome.measure == "Odds Ratio") %>% filter(!is.na(study.year)) %>% 
  group_by(file.nr, comparison.nr, subgroup.nr) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% 
  mutate(first.effect = effect[which.min(study.year)]) %>% 
  filter(first.effect > 1) %>% 
  mutate(scaled.effect = scale(effect), scaled.time = scale(study.year)) %>% 
  ggplot(aes(x = scaled.time, y = scaled.effect)) + geom_point(size = 1, alpha = 0.7)  + 
  geom_smooth(method = lm) + theme_bw()  + 
  labs(title = "Scaled and centered OR over time", subtitle = "Initial finding > 1") + ylab("Odds Ratio")  + ylab("Time")

#Scaled OR over time with initial OR < 1
madata %>% filter(outcome.measure == "Odds Ratio") %>% filter(!is.na(study.year)) %>% 
  group_by(file.nr, comparison.nr, subgroup.nr) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% 
  mutate(first.effect = effect[which.min(study.year)]) %>% 
  filter(first.effect < 1) %>% 
  mutate(scaled.effect = scale(effect), scaled.time = scale(study.year)) %>% 
  ggplot(aes(x = scaled.time, y = scaled.effect)) + geom_point(size = 1.2, alpha = 0.2)  + 
  geom_smooth(method = lm) + theme_bw() + 
  labs(title = "Scaled and centered OR over time", subtitle = "Initial finding < 1") + ylab("Odds Ratio")  + ylab("Time")

#Scaled RR over time with initial RR > 1
madata %>% filter(outcome.measure == "Risk Ratio") %>% filter(!is.na(study.year)) %>% 
  group_by(file.nr, comparison.nr, subgroup.nr) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% 
  mutate(first.effect = effect[which.min(study.year)]) %>% 
  filter(first.effect > 1) %>% 
  mutate(scaled.effect = scale(effect), scaled.time = scale(study.year)) %>% 
  ggplot(aes(x = scaled.time, y = scaled.effect)) + geom_point(size = .5, alpha = 0.2)  + 
  geom_smooth(method = lm) + theme_bw() + 
  labs(title = "Scaled and centered RR over time", subtitle = "Initial finding > 1") + ylab("Risk Ratio")  + ylab("Time")

#Scaled RR over time with initial RR > 1
madata %>% filter(outcome.measure == "Risk Ratio") %>% filter(!is.na(study.year)) %>% 
  group_by(file.nr, comparison.nr, subgroup.nr) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% 
  mutate(first.effect = effect[which.min(study.year)]) %>% 
  filter(first.effect < 1) %>% 
  mutate(scaled.effect = scale(effect), scaled.time = scale(study.year)) %>% 
  ggplot(aes(x = scaled.time, y = scaled.effect)) + geom_point(size = .5, alpha = 0.1)  + 
  geom_smooth(method = lm) + theme_bw() + 
  labs(title = "Scaled and centered RR over time", subtitle = "Initial finding < 1") + ylab("Risk Ratio")  + ylab("Time")

#Regression to test if an interaction initial effect estimate >/< 1 and the slope is significant
regression.temp <- madata %>% filter(outcome.measure == "Risk Ratio") %>% filter(!is.na(study.year)) %>% 
  group_by(file.nr, comparison.nr, subgroup.nr) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% 
  mutate(first.effect = effect[which.min(study.year)]) %>% 
  mutate(effect.sign = ifelse(first.effect > 1, "pos", "neg"), scaled.effect = scale(effect), scaled.time = scale(study.year))
regression.temp$effect.sign <- as.factor(regression.temp$effect.sign)


m1 <- lm(scaled.effect ~ scaled.time, data = regression.temp)
m2 <- lm(scaled.effect ~ scaled.time * effect.sign, data = regression.temp)
anova(m1, m2)


#Histogram of scaled effect sizes
madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% mutate(scaled.effect = scale(effect), scaled.time = scale(study.year)) %>% 
  ggplot(aes(x = scaled.effect)) + geom_histogram(col = "gray15", fill = "dodgerblue", bins = 40) + theme_bw() + labs(title = "Scaled effect sizes")

#Histogram of effect sizes
madata %>%  filter(effect > -10 & effect < 10) %>% ggplot(aes(x = effect)) + geom_histogram(col = "gray15", fill = "dodgerblue", bins = 100) + 
  theme_bw() + labs(title = "Effect sizes")


# #Checks:
# madata %>% filter(file.nr == 21) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
#   mutate(counts = n()) %>% filter(counts > 1) %>% mutate(scaled.effect = scale(effect), scaled.time = scale(study.year)) %>% 
#   select(scaled.effect, effect, scaled.time, study.year, study.name)
# madata %>% filter(file.nr == 21) %>% group_by(file.nr, outcome.name, comparison.name, subgroup.nr) %>% 
#   mutate(counts = n()) %>% filter(counts > 1) %>% mutate(scaled.effect = scale(effect), scaled.time = scale(study.year)) %>% 
#   ggplot(aes(x = scaled.time, y = scaled.effect, group = outcome.nr, colour = outcome.nr)) + geom_line() 

##################################################################################################################################
#Heterogeneity tests for reviews:
#Transformation to log ORs and RRs - calculation of the according se()
madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% mutate(log.effect = log(effect)) 

#Transformation to log ORs- calculation of the according se()
madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  filter(outcome.measure == "Odds Ratio") %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>%
  mutate(log.OR = log(effect), 
         se.log.OR = sqrt( (1/events1) + ifelse(total1 - events1 > 0, (1/(total1 - events1)), 0) + 
                             (1/events2) + ifelse(total2 - events2 > 0, (1/(total2 - events2)), 0))) %>% 
  filter(!is.na(se.log.OR)) %>% mutate(summary.effect = sum( log.OR/se.log.OR^2 ) / sum( 1/se.log.OR^2 )) %>% 
  summarise(Q = sum(((log.OR - summary.effect)^2)/se.log.OR^2), n = n(), OR =  exp(min(summary.effect))) %>% 
  mutate(p.value = 1 - pchisq(Q, n - 1)) %>% 
  ggplot(aes(x = p.value)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Cochrane Heterogeneity Test P-values of ORs")

madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  filter(outcome.measure == "Odds Ratio") %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>%
  summarise(Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, method = "Inverse", sm = "OR")$Q,
            p.value = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, method = "Inverse", sm = "OR")$pval.Q) %>% 
  ggplot(aes(x = p.value)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Cochrane Heterogeneity Test P-values of ORs")


#Check if identical: (almost ..)
file.number.to.compare <- 23
ds <- madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% filter(file.nr == file.number.to.compare) %>% 
  filter(outcome.measure == "Odds Ratio") %>% filter(events1 > 0 & events2 > 0) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>%
  summarise(
    Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, method = "Inverse", sm = "OR")$Q,
    p.value = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, 
                      method = "Inverse", sm =  "OR")$pval.Q)

ds2 <- madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% filter(file.nr == file.number.to.compare) %>% 
  filter(outcome.measure == "Odds Ratio") %>% filter(events1 > 0 & events2 > 0) %>%
  mutate(counts = n()) %>% filter(counts > 1) %>%
  mutate(log.OR = log(effect), 
         se.log.OR = sqrt( (1/events1) + ifelse(total1 - events1 > 0, (1/(total1 - events1)), 0) + 
                             (1/events2) + ifelse(total2 - events2 > 0, (1/(total2 - events2)), 0))) %>% 
  filter(!is.na(se.log.OR)) %>% mutate(summary.effect = sum( log.OR/se.log.OR^2 ) / sum( 1/se.log.OR^2 )) %>% 
  summarise(Q = sum(((log.OR - summary.effect)^2)/se.log.OR^2), n = n(), OR =  exp(min(summary.effect))) %>% 
  mutate(p.value = 1 - pchisq(Q, df = n - 1))

identical(ds$Q, ds2$Q)
print(data_frame(ds$Q, ds2$Q, ds$p.value, ds2$p.value), n = 20)

# #Check: 
# madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
#   filter(outcome.measure == "Odds Ratio") %>% 
#   filter(file.nr == 11) %>% 
#   mutate(counts = n()) %>% filter(counts > 1) %>%
#   mutate(log.OR = log(effect), 
#          se.log.OR = sqrt( (1/events1) + (1/(total1 - events1)) + (1/events2) + (1/(total2 - events2)))) %>% 
#   filter(!is.na(se.log.OR)) %>% mutate(summary.effect = sum( log.OR/se.log.OR^2 ) / sum( 1/se.log.OR^2 )) %>% 
#   summarise(Q = sum(((log.OR - summary.effect)^2)/se.log.OR^2), n = n(), OR =  exp(min(summary.effect))) %>% 
#   mutate(p.value = 1 - pchisq(Q, n - 1))
# 
# 
# madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% filter(outcome.measure == "Odds Ratio") %>% 
#   mutate(counts = n()) %>% filter(counts > 1 & !is.na(effect) & events1 > 0 & events2 > 0 & !is.na(events1) & !is.na(events2)) %>% 
#   mutate(log.OR = log(effect), se.log.OR = sqrt(1/events1 + 1/(total1 - events1) + 1/events2 + 1/(total1 - events2)))
# 
# madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
#   filter(outcome.measure == "Odds Ratio") %>% 
#   mutate(counts = n()) %>% filter(counts > 1) %>% filter(!is.na(effect)) %>% 
#   summarise(isna = sum(is.na(effect)), isna1 = sum(is.na(events1)), isna2 = sum(is.na(events2))) %>%
#   ungroup %>% count(isna, isna1, isna2)


#Q heterogeneity statistic calculation (manually)
yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% filter(file.nr < 100) %>% 
  mutate(summary.effect = sum(fishersz/variance ) / sum( 1/variance )) %>% 
  summarise(Q = sum(((fishersz - summary.effect)^2)/variance ), n = n()) %>%
  mutate(p.value = 1 - pchisq(Q, n - 1)) %>% 
  ggplot(aes(x = p.value)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Cochrane Heterogeneity Test P-values of ORs")

yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% filter(file.nr < 100) %>% 
  summarise(Q = metacor(cor = fishersz, n = total1 + total2, sm = "ZCOR")$Q,
            p.value = metacor(cor = fishersz, n = total1 + total2, sm = "ZCOR")$pval.Q) %>% 
  ggplot(aes(x = p.value)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Cochrane Heterogeneity Test P-values of ORs")


##################################################################################################################################

#Frequencies of trials with the same research subject per review
yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
  filter(n < 100) %>%
  full_join( yodata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
               filter(n > 99) %>% ungroup %>% summarise(n = 100, nn = sum(nn))) %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Trials Comparing the Same Subject per Review")

yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
  ungroup %>% arrange(desc(nn)) %>% mutate(trials = cumsum(n * nn)) %>%
  ggplot(aes(x = nn, y = trials)) + geom_line(col = "gray15") +
  theme_bw() + labs(title = "Trials Comparing the Same Subject per Review")

##################################################################################################################################

#Reproduction probability 
print(yodata %>% filter(file.nr == 21) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
        distinct(study.year, .keep_all = T) %>% 
        mutate(counts = n()) %>% filter(counts > 1) %>% mutate(time.rank = rank(study.year)) %>% 
        mutate(pval = 2*(1-pnorm(abs(fishersz), mean = 0, sd = sqrt(variance)))) %>% 
        select(fishersz, study.year, time.rank, pval) %>% filter(time.rank == 1 | time.rank == 2), n = 40) 


temp.pval <- yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
  distinct(study.year, .keep_all = T) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% mutate(time.rank = rank(study.year)) %>% 
  mutate(pval = 2*(1-pnorm(abs(fishersz), mean = 0, sd = sqrt(variance))))  

temp.pval %>% filter(time.rank < 10) %>% ggplot(aes(x = time.rank, y = pval)) + 
  geom_jitter(col = "gray15", size = 0.3, alpha = .25) + 
  theme_bw() + xlab("Time rank") + ylab("P-value") + labs(title = "P-values over (ranked) time")

temp.pval %>% filter(time.rank < 10) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  mutate(diff = pval - pval[time.rank - 1]) %>% select(diff, pval, time.rank) %>% filter(time.rank != 1) %>% 
  ggplot(aes(x = time.rank, y = diff)) + 
  geom_jitter(col = "gray15", size = 0.3, alpha = .25) + 
  theme_bw() + xlab("Significance of first trial") + ylab("Significance of second trial")

temp.pval %>% summarise(first.p = pval[time.rank == 1], second.p = pval[time.rank == 2]) %>% 
  ggplot(aes(x = first.p, y = second.p)) + geom_point(col = "gray15", size = 0.1) + 
  theme_bw() + xlab("Significance of first trial") + ylab("Significance of second trial")

temp.pval %>% summarise(first.p = pval[time.rank == 1], second.p = pval[time.rank == 2]) %>% 
  mutate(repr = ifelse(first.p > second.p, "rep.yes", "rep.no")) %>%
  ggplot(aes(x = repr)) + geom_histogram(col = "gray15", fill = "dodgerblue", stat = "count") + 
  theme_bw() + xlab("Reproduction success")

yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
  distinct(study.year, .keep_all = T) %>% mutate(time.rank = rank(study.year)) %>% 
  mutate(pval = 2*(1-pnorm(abs(fishersz), mean = 0, sd = sqrt(variance)))) %>% 
  ungroup() %>% count(time.rank == 2)

yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
  distinct(study.year, .keep_all = T) %>% mutate(time.rank = rank(study.year)) %>% 
  mutate(pval = 2*(1-pnorm(abs(fishersz), mean = 0, sd = sqrt(variance)))) %>% 
  summarise(first.p = pval[time.rank == 1], second.p = pval[time.rank == 2]) %>% 
  filter(!is.na(first.p))
ggplot(aes(x = first.p, y = second.p)) + geom_point(col = "gray15") + 
  theme_bw() + xlab("Significance of first trial") + ylab("Significance of second trial")

yodata %>% filter(file.nr == 21) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
  distinct(study.year, .keep_all = T) %>% mutate(time.rank = rank(study.year)) %>% 
  mutate(pval = 2*(1-pnorm(abs(fishersz), mean = 0, sd = sqrt(variance)))) %>% 
  summarise(first.p = pval[time.rank == 1], second.p = pval[time.rank == 2]) %>% 
  ggplot(aes(x = first.p, y = second.p)) + geom_point(col = "gray15") + 
  theme_bw()

yodata %>% filter(file.nr == 21) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
  distinct(study.year, .keep_all = T) %>% mutate(time.rank = rank(study.year)) %>% 
  mutate(pval = 2*(1-pnorm(abs(fishersz), mean = 0, sd = sqrt(variance)))) %>% select(pval[time.rank == 1])

##################################################################################################################################

#Leave one out fixed effect prediction (manually):
yodata.temp.fixef <- yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  filter(file.nr < 1000) %>% 
  filter(!is.na(fishersz) & variance > 0) %>% mutate(counts = n()) %>% filter(counts > 2) %>% 
  mutate(summary.effect = (sum( fishersz/variance) - fishersz/variance) / (sum( 1/variance ) - 1/variance), 
         summary.variance =  (sum(variance ) - variance))
levels <- seq(0.5, 0.99, by = 0.01)
e.cs <- c()

for(level in levels){
  empirical.coverage <- yodata.temp.fixef %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  
    mutate(l.pred = summary.effect - qnorm( (1+level)/2 ) * sqrt(summary.variance),
           u.pred = summary.effect + qnorm( (1+level)/2 ) * sqrt(summary.variance)) %>% 
    arrange(study.year) %>%
    mutate(covered = ifelse(fishersz < u.pred & fishersz > l.pred, 1, 0), id = row_number()) %>%
    filter(row_number() > 1) %>% 
    ungroup() %>% summarise(sum(covered)/length(covered))
  e.cs <- c(e.cs, empirical.coverage)
}

plot(levels, as.numeric(e.cs), ylab = "nominal coverage", xlab = "empirical coverage", ylim = c(.5, 1), xlim = c(.5,1))
lines(c(0,1),c(0,1), col = 11)

##################################################################################################################################

#Leave one out random effects prediction:
yodata.temp.ranef <- yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  filter(file.nr < 1000) %>% 
  filter(!is.na(fishersz) & variance > 0) %>% mutate(counts = n()) %>% filter(counts > 2) %>% 
  mutate(summary.effect = (sum( fishersz/variance) - fishersz/variance) / (sum( 1/variance ) - 1/variance), 
         summary.variance =  (sum(variance ) - variance),
         Q = sum(((fishersz - summary.effect)^2)/variance) - ((fishersz - summary.effect)^2)/variance, n = n(), 
         tau  = max(c(0, Q - (n - 2) / (sum(1/variance) - 1/variance - (sum(variance^2) - variance^2)/(sum(variance) - variance)))),
         summary.ranef = (sum(fishersz/(variance + tau)) - fishersz/(variance + tau))/ (sum(1/(variance + tau)) - 1/(variance-tau)),
         variance.ranef = sum(variance + tau) - (variance + tau)) 

# yodata.temp.ranef <- yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  filter(file.nr < 1000) %>% 
#   filter(!is.na(fishersz) & variance > 0) %>% mutate(counts = n()) %>% filter(counts > 2) %>% 
#   mutate(summary.effect = (sum( fishersz/variance) - fishersz/variance) / (sum( 1/variance ) - 1/variance), 
#          summary.variance =  (sum(variance ) - variance),
#          Q = sum(((fishersz - summary.effect)^2)/variance) - ((fishersz - summary.effect)^2)/variance, n = n(), 
#          tau  = max(c(0, Q - (n - 2) / (sum(1/variance) - 1/variance - (sum(1/(variance^2)) - 1/(variance^2))/(sum(1/variance) - 1/variance)))),
#          summary.ranef = (sum(fishersz/(variance + tau)) - fishersz/(variance + tau))/ (sum(1/(variance + tau)) - 1/(variance-tau)),
#          variance.ranef = sum(variance + tau) - (variance + tau))

# #Check1 (no leave one out, compare calculations). Tau is different, but possibly due to different methods.
# yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  
#   filter(file.nr == 22, outcome.nr == 2, comparison.nr == 1) %>% select(fishersz, variance)
# 
# tp.s <- yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  
#   filter(file.nr == 22, outcome.nr == 2, comparison.nr == 1) %>% 
#   slice(c(2:4)) 
# 
# 
# tp.s %>% summarise(summary.ranef = metacor(cor = fishersz, n = total1 + total2, sm = "COR")$TE.random, 
#                    se.ranef = metacor(cor = fishersz, n = total1 + total2, sm = "COR")$seTE.random, 
#                    tau = metacor(cor = fishersz, n = total1 + total2, sm = "COR")$tau, 
#                    Q = metacor(cor = fishersz, n = total1 + total2, sm = "COR")$Q)
# 
# tp.s %>% mutate( summary.effect = (sum( fishersz/variance)) / (sum( 1/variance )), 
#                  summary.variance =  (sum(variance)),
#                  Q = sum(((fishersz - summary.effect)^2)/variance), 
#                  n = n(),  
#                  tau  = max(c(0, (Q - (n - 1)) / (sum(1/variance) - (sum(1/variance^2))/(sum(1/variance))))),
#                  summary.ranef = (sum(fishersz/(variance + tau)))/ (sum(1/(variance + tau))),
#                  se.ranef = 1/sqrt(sum(1/(variance + tau)))) %>% 
#   select(summary.ranef, se.ranef, tau, Q) 
# 
# yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
#   filter(outcome.measure == "Odds Ratio") %>% 
#   mutate(counts = n()) %>% filter(counts > 1) %>%
#   mutate(log.OR = log(effect), 
#          se.log.OR = sqrt( (1/events1) + (1/(total1 - events1)) + (1/events2) + (1/(total2 - events2)))) %>% 
#   filter(!is.na(se.log.OR)) %>% mutate(summary.effect = sum( log.OR/se.log.OR^2 ) / sum( 1/se.log.OR^2 )) %>% 
#   summarise(Q = sum(((log.OR - summary.effect)^2)/se.log.OR^2), n = n(), OR =  exp(min(summary.effect)))
# 
# 
# levels <- seq(0.5, 0.99, by = 0.01)
# e.cs <- c()
# 
# for(level in levels){
#   empirical.coverage <- yodata.temp.ranef %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  
#     mutate(l.pred = summary.effect - qnorm( (1+level)/2 ) * sqrt(summary.variance),
#            u.pred = summary.effect + qnorm( (1+level)/2 ) * sqrt(summary.variance)) %>% 
#     arrange(study.year) %>%
#     mutate(covered = ifelse(fishersz < u.pred & fishersz > l.pred, 1, 0), id = row_number()) %>%
#     filter(row_number() > 1) %>% 
#     ungroup() %>% summarise(sum(covered)/length(covered))
#   e.cs <- c(e.cs, empirical.coverage)
# }
# 
# plot(levels, as.numeric(e.cs), ylab = "nominal coverage", xlab = "empirical coverage", ylim = c(.5, 1), xlim = c(.5,1))
# lines(c(0,1),c(0,1), col = 11)
# 
# 
# 
# yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  filter(file.nr < 100) %>% 
#   filter(!is.na(fishersz)) %>% mutate(counts = n()) %>% filter(counts > 2) %>% 
#   mutate(summary.effect = (sum( fishersz/variance) - fishersz/variance) / (sum( 1/variance ) - 1/variance), 
#          summary.variance =  (sum(variance ) - variance),
#          l.pred = summary.effect - qnorm( (1+level)/2 ) * sqrt(summary.variance),
#          u.pred = summary.effect + qnorm( (1+level)/2 ) * sqrt(summary.variance)) %>% 
#   mutate(covered = ifelse(fishersz < u.pred & fishersz > l.pred, 1, 0)) %>% select(variance, fishersz, l.pred, u.pred, summary.effect, summary.variance, covered)
# 
# yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%  filter(file.nr < 100) %>% 
#   filter(!is.na(fishersz)) %>% mutate(counts = n()) %>% filter(counts > 2) %>% 
#   mutate(summary.effect = (sum( fishersz/variance) - fishersz/variance) / (sum( 1/variance ) - 1/variance), 
#          summary.variance =  (sum(variance ) - variance),
#          l.pred = summary.effect - qnorm( (1+level)/2 ) * sqrt(summary.variance),
#          u.pred = summary.effect + qnorm( (1+level)/2 ) * sqrt(summary.variance)) %>% 
#   mutate(covered = ifelse(fishersz < u.pred & fishersz > l.pred, 1, 0)) %>%
#   ungroup() %>% summarise(sum(covered)/length(covered))
# 
# 
# yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
#   filter(!is.na(fishersz)) %>% mutate(counts = n()) %>% filter(counts > 2) %>%
#   summarise(Q = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, method = "Inverse", sm = "OR")$Q,
#             p.value = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, studlab = study.name, method = "Inverse", sm = "OR")$pval.Q) %>% 
#   ggplot(aes(x = p.value)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
#   theme_bw() + labs(title = "Cochrane Heterogeneity Test P-values of ORs")

##################################################################################################################################

#Checks: Is number of fishersz outcomes appr. equal to RR + OR + Meandifferences
# (n <- length(madata$study.name))
# 
# (n.z <- sum(!is.na(madata$fishersz)))
# 
# (n.b <- sum(madata$events1 > 0)) #binary data number
# (n.d <- sum(madata$outcome.measure == "Mean Difference")) #mean diff. number
# 
# n.z - n.b - n.d #27000 
# 
# total <- madata$total1 + madata$total2
# hist(total[total < 500], 100)

# sum(!is.na(madata$std.mean.d))
# sum(!is.na(madata$correlation))
# sum(!is.na(madata$std.mean.d)) - sum(!is.na(madata$correlation))
# #41971 std.mean.d are lost..
# sum(is.na(madata$total1[!is.na(madata$std.mean.d)] * madata$total2[!is.na(madata$std.mean.d)])) #198 occurences due to too large product -> 41773
# nas <- ifelse(madata$total1[!is.na(madata$std.mean.d)] > 0, 0, 1) + ifelse(madata$total2[!is.na(madata$std.mean.d)] > 0, 0, 1) - 
#   ifelse(madata$std.mean.d[!is.na(madata$std.mean.d)] != 0, 2, 0)
# sum(nas > 0)                     #No occurences where denominator is = 0
# sum(!is.na(madata$std.mean.d[!is.na(madata$std.mean.d)]/sqrt((madata$std.mean.d[!is.na(madata$std.mean.d)]^2))))
# sum(!is.na(1/sqrt((madata$std.mean.d[!is.na(madata$std.mean.d)]^2))))
# sum(is.na((madata$total2*madata$total1)))
# sum(!is.na(madata$std.mean.d))
# sum(!is.na(madata$correlation))
# sum(!is.na(madata$std.mean.d)) - sum(!is.na(madata$correlation))
# sum(madata$std.mean.d[!is.na(madata$std.mean.d)] > 0)


# #Check: Expected number of fishersz outcomes  != 0
# sum(madata$events1 > 0) - sum(is.na(madata$total1 * madata$total2))
# sum(!is.na(madata$odds.ratio))
# sum(!is.na(madata$std.mean.d))
# sum(!is.na(madata$correlation))
# sum(!is.na(madata$fishersz))

##################################################################################################################################

##################################################################################################################################

##################################################################################################################################

##################################################################################################################################