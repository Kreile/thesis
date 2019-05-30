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






################################################################################################################
################################################################################################################

#Look at Migraine sympton relieve review
arrange(filter(data, file.name=="MYdwnld_CD005220StatsDataOnly_Version2.csv") %>% 
          select(file.nr, comparison.nr, comparison.name, outcome.nr, subgroup.nr, study.name), study.name)

arrange(filter(data, file.nr == 21) %>% 
          select(file.nr, comparison.nr, outcome.nr, subgroup.nr, study.name), study.name)

arrange(filter(data, file.nr == 21) %>% 
          select(file.nr, comparison.name, outcome.name, subgroup.name, study.name), study.name)


arrange(filter(data, file.nr == 761) %>% 
          select(study.name, comparison.name, outcome.name), study.name)

arrange(filter(data, file.nr == 2897) %>% 
          select(study.name, comparison.name, outcome.name), study.name)


#How many studies are in the dataset?
data %>% group_by(file.nr) %>% distinct(study.name) %>% ungroup %>% count()

#Distribution of number of distinct studies per review (1: >=30):
data %>% group_by(file.nr) %>% distinct(study.name) %>% count() %>% 
  filter(n < 30) %>% 
  full_join( data %>% group_by(file.nr) %>% distinct(study.name) %>% count %>% 
               filter(n > 29) %>% mutate(n = 30)) %>% 
  group_by(n) %>% count() %>% 
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Histogram of Study Number per Review") + xlab("Study number per review")
#Finding review with one study only:
#data %>% group_by(file.nr) %>% distinct(study.name) %>% count() %>% filter(n < 2) -> File.nr 38 for example
arrange(filter(data, file.nr == 38) %>% 
          select(doi, study.name, comparison.name, outcome.name), study.name)

#Number of different research subjects per review:
data %>% group_by(file.nr) %>% distinct(comparison.name) %>% count() %>%
  filter(n < 50) %>%
  full_join( data %>% group_by(file.nr) %>% distinct(comparison.name) %>% count %>%
               filter(n > 49) %>% mutate(n = 50)) %>%
  group_by(n) %>% count() %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Research Subject Number per Review") + xlab("Research subject number") + ylab("Number of reviews")

#Number of different outcomes for different research subjects summed up per review:
data %>% group_by(file.nr, comparison.nr) %>% distinct(outcome.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% filter(n < 50) %>%
  full_join( data %>% group_by(file.nr, comparison.nr) %>% distinct(outcome.name) %>% ungroup() %>% group_by(file.nr) %>% 
               count() %>% filter(n > 49) %>% mutate(n = 50)) %>%
  group_by(n) %>% count() %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Outcome Number per Review") + xlab("Outcome number") + ylab("Number of reviews")

#Mean and median number of different outcomes for different research subjects summed up per review:
data %>% group_by(file.nr, comparison.nr) %>% distinct(outcome.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% ungroup %>% summarise(mean = mean(n))

data %>% group_by(file.nr, comparison.nr) %>% distinct(outcome.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% ungroup %>% summarise(mean = median(n))

#Number of different subgroups given that they share research subject and outcome:
data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% filter(n < 100) %>%
  full_join( data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
               ungroup() %>% group_by(file.nr) %>% 
               count() %>% filter(n > 99) %>% mutate(n = 100)) %>%
  group_by(n) %>% count() %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Subgroup Number per Review") + xlab("Subgroup number") + ylab("Number of reviews")

data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% 
  ggplot(aes(x = n)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Subgroup Number per Review") + xlab("Subgroup number") + ylab("Number of reviews")

data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% filter(n < 100) %>%
  full_join( data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
               ungroup() %>% group_by(file.nr) %>% 
               count() %>% filter(n > 99) %>% mutate(n = 100)) %>%
  ggplot(aes(x = n)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Subgroup Number per Review") + xlab("Subgroup number") + ylab("Number of reviews")

mean.diff.subg <- data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% ungroup() %>% summarise(mean = mean(n))

median.diff.subg <- data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% ungroup() %>% summarise(median = median(n))

#Barbiturate Review:
arrange(filter(data, file.nr == 21) %>% 
select(file.nr, comparison.name,  subgroup.name, outcome.name, study.name), study.name)

#Frequencies of outcome type among study results:
data %>% group_by(outcome.measure.new) %>% count() %>% arrange(desc(n)) %>% ungroup() %>% filter(row_number() < 9) %>% 
  ggplot(aes(x = outcome.measure, y = n)) +
  geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Outcome measure frequencies") + xlab("Frequency") + ylab("Outcome measure")

#Study year dsitribution:
data %>% filter(study.year < 2019 & study.year > 1920) %>% ggplot(aes(x = study.year)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Study publication year frequency") + xlab("Year") + ylab("Frequency")

#The frequency of missing or erronous variables is also important. Especially, missing outcome measure indication, sample size indication, 
#standard deviation, mean and event indication for the cases where the latter are meaningful.
missing.mean <- data %>% filter(!is.na(sd1) & sd1 > 0) %>% filter(!is.na(sd2) & sd2 > 0) %>% 
  summarise(total = length(sd1), na  = sum(is.na(mean1) | is.na(mean1))) %>% mutate(f = na/total) %>% select(f)

missing.sd <- data %>% filter(!is.na(mean1) & !is.na(mean2)) %>% filter(mean1 != 0 & mean2 != 0) %>%
  summarise(total = length(mean1), na.sd  = sum(is.na(sd1) | is.na(sd1)), null.sd = length(sd1[sd1 <= 0 | sd2 <= 0])) %>% 
  mutate(f = (na.sd + null.sd)/total) %>% select(f)

missing.ssize <- data %>% summarize(total = length(total1), null.ss = length(total1[total1 <= 0 | total2 <= 0])) %>% 
  mutate(f = null.ss / total) %>% select(f)

missing.events <- data %>% filter(is.na(mean1) | mean1 == 0) %>% filter(is.na(mean2) | mean2 == 0) %>% 
  filter(is.na(sd1) | sd1 <= 0) %>% filter(is.na(sd2) | sd2 <= 0) %>% 
  summarise(null.ev = length(total1[events1 <= 0 | events2 <= 0]), total = length(events1)) %>% mutate(f = null.ev/total) %>% select(f)

missing.mean <- data %>% filter(!is.na(sd1) & sd1 > 0) %>% filter(!is.na(sd2) & sd2 > 0) %>% 
  summarise(total = length(sd1), na  = sum(is.na(mean1) | is.na(mean1))) %>% mutate(f = na/total) %>% select(f)

missing.sd <- data %>% filter(!is.na(mean1) & !is.na(mean2)) %>% filter(mean1 != 0 & mean2 != 0) %>%
  summarise(total = length(mean1), na.sd  = sum(is.na(sd1) | is.na(sd1)), null.sd = length(sd1[sd1 <= 0 | sd2 <= 0])) %>% 
  mutate(f = (na.sd + null.sd)/total) %>% select(f)

missing.ssize <- data %>% summarize(total = length(total1), null.ss = length(total1[total1 <= 0 | total2 <= 0])) %>% 
  mutate(f = null.ss / total) %>% select(f)

null.events <- data %>% filter(is.na(mean1) | mean1 == 0) %>% filter(is.na(mean2) | mean2 == 0) %>% 
  filter(is.na(sd1) | sd1 <= 0) %>% filter(is.na(sd2) | sd2 <= 0) %>% 
  summarise(null.ev = length(total1[events1 <= 0 & events2 <= 0]), total = length(events1)) %>% 
  mutate(f = null.ev/total) %>% select(f)

null.events.single <- data %>% filter(is.na(mean1) | mean1 == 0) %>% filter(is.na(mean2) | mean2 == 0) %>% 
  filter(is.na(sd1) | sd1 <= 0) %>% filter(is.na(sd2) | sd2 <= 0) %>% 
  summarise(null.ev = length(total1[events1 <= 0 | events2 <= 0]), total = length(events1)) %>% 
  mutate(f = null.ev/total) %>% select(f)

missing.events <- data %>% filter(is.na(mean1) | mean1 == 0) %>% filter(is.na(mean2) | mean2 == 0) %>% 
  filter(is.na(sd1) | sd1 <= 0) %>% filter(is.na(sd2) | sd2 <= 0) %>% 
  summarise(na.ev = sum(is.na(events1) | is.na(events2)), total = length(events1)) %>% 
  mutate(f = na.ev/total) %>% select(f)
  

#Number of groups with n repr trials:
data %>% group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>%
  filter(n < 20) %>%
  full_join( data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
               filter(n > 19) %>% ungroup %>% summarise(n = 20, nn = sum(nn))) %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Number of groups with number of reproduction trials = n") + xlab("n") + ylab("Number of groups")

#Cumulative number of groups with n repr trials:
data %>% group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>%
  filter(n < 20) %>%
  full_join( data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
               filter(n > 19) %>% ungroup %>% summarise(n = 20, nn = sum(nn))) %>% 
  ungroup() %>% arrange(desc(n)) %>% mutate(csum  = cumsum(nn)) %>% 
  ggplot(aes(x = n, y = csum)) + geom_point(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Number of groups with number of reproduction trials >= n") + xlab("n") + ylab("Number of groups")


#Frequencies of unique combinations of review, compared research subject and outcome measure 
# -> Frequencies of trials analyzing the same subject in a review (with the same outcome)
data %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count 

data %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count %>% filter(n < 15) %>% 
  full_join( data %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 14) %>% ungroup %>% summarise(n = 15, nn = sum(nn)))

data %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count %>% 
  filter(n < 15) %>% 
  full_join( data %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count %>% 
    filter(n > 14) %>% ungroup %>% summarise(n = 15, nn = sum(nn))) %>%  
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Trials comparing the same subject per Review")

#Frequencies of unique combinations of review, compared research subject, outcome measure and subgroup 
# -> Frequencies of trials analyzing the same subject with the same subgroups in a review
data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count 

data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% filter(n < 15) %>% 
  full_join( data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 14) %>% ungroup %>% summarise(n = 15, nn = sum(nn)))

data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% filter(n < 15) %>% 
  full_join( data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 14) %>% ungroup %>% summarise(n = 15, nn = sum(nn))) %>%  
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Trials comparing the same subject in the same subgroups per Review")

#Frequencies of different outcomes per review
data %>% group_by(file.nr) %>% summarise(diff.out = length(unique(outcome.name))) %>% count(diff.out)

data %>% group_by(file.nr) %>% summarise(diff.out = length(unique(outcome.name))) %>% count(diff.out) %>%  
  ggplot(aes(x = diff.out, y = n)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs("Unique Comparisons per Review")

#Study year distribution:
data %>% filter(study.year > 1900 & study.year < 2019) %>% ggplot(aes(x = study.year)) + geom_histogram()

#Study year distribution of distinct study years per review (i.e. duplicates are not considered)
data %>% group_by(file.nr) %>% filter(study.year > 1900 & study.year < 2019) %>% distinct(study.year) %>% 
  ggplot(aes(x = study.year)) + 
  geom_histogram(col = "gray15", fill = "dodgerblue") + theme_bw() + labs(title = "Study years of analyses")

#Scaled effect sizes over time
data %>% filter(file.nr < 5016) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% mutate(scaled.effect = scale(effect), scaled.time = scale(study.year)) %>% 
  ggplot(aes(x = scaled.time, y = scaled.effect)) + geom_point(size = .005, alpha = 0.15) 

data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  filter(study.year > 1950 | study.year < 2019) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% 
  filter(max(diff(study.year)) < 50) %>% 
  mutate(abs.scaled.effect = abs(scale(effect)), scaled.time = scale(study.year)) %>% 
  ggplot(aes(x = scaled.time, y = abs.scaled.effect)) + geom_point(size = .005, alpha = 0.15) + geom_smooth()

#Pvalues over time
temp.pval <- data.ext %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
  distinct(study.year, .keep_all = T) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% mutate(time.rank = rank(study.year))

temp.pval %>% filter(time.rank < 10) %>% ggplot(aes(x = time.rank, y = pval.single)) + 
  geom_jitter(col = "gray15", size = 0.3, alpha = .25) + 
  theme_bw() + xlab("Time rank") + ylab("P-value") + labs(title = "P-values over (ranked) Time") +geom_smooth(method = "lm")

########################################################################################################################################
########################################################################################################################################
# Cumulative number of total trials with
data %>% filter(!is.na(fishersz)) %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
  ungroup %>% arrange(desc(n)) %>% mutate(trials = cumsum(n * nn)) %>% 
  filter(n < 50) %>% 
  full_join( data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 49) %>% ungroup %>% summarise(n = 50, nn = sum(nn))) %>%  
  ggplot(aes(x = n, y = trials)) + geom_point(col = "gray15") + 
  theme_bw() + labs(title = "Sum of Trials with Number of Replicates >= n (ignoring subgroups") + ylab("Number of trials") + xlab("n, number of replicates")

data %>% filter(!is.na(fishersz)) %>% group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>%
  ungroup %>% arrange(desc(n)) %>% mutate(trials = cumsum(n * nn)) %>% 
  filter(n < 50) %>% 
  full_join( data %>% group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 49) %>% ungroup %>% summarise(n = 50, nn = sum(nn))) %>%  
  ggplot(aes(x = n, y = trials)) + geom_point(col = "gray15") + 
  theme_bw() + labs(title = "Sum of Trials with Number of Replicates >= n") + ylab("Number of trials") + xlab("n, number of replicates")

data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
  filter(n < 100) %>%
  full_join( data %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
               filter(n > 99) %>% ungroup %>% summarise(n = 100, nn = sum(nn))) %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Number of groups with number of reproduction trials = n") + xlab("n") + ylab("Number of groups")

data %>% filter(!is.na(fishersz)) %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
  ungroup %>% arrange(desc(nn)) %>% mutate(trials = cumsum(n * nn)) %>%
  ggplot(aes(x = nn, y = trials)) + geom_line(col = "gray15") +
  theme_bw() + labs(title = "Trials Comparing the Same Subject per Review")

########################################################################################################################################
########################################################################################################################################
#Look for sample size effect on effect size
########################################################################################################################################
########################################################################################################################################

data %>% filter(total1 + total2 > 15 & total1 + total2 < 500) %>% 
  mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>%
  group_by(sample.size) %>% 
  summarize(mean.effect = mean(scaled.effect, na.rm = T)) %>% 
  ggplot(aes(x = sample.size, y = mean.effect)) + geom_point() + theme_bw() + xlab("sample size") + ylab("mean absolute normalized effect size")

data %>% filter(total1 + total2 > 10 & total1 + total2 < 100) %>% 
  mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>%
  group_by(sample.size) %>% 
  summarize(median.effect = median(scaled.effect, na.rm = T)) %>% 
  ggplot(aes(x = sample.size, y = median.effect)) + geom_point() + theme_bw() + xlab("sample size") + ylab("median absolute normalized effect size")

#Plot separately for outcome measures (median effect size)
data %>% filter(total1 + total2 > 10 & total1 + total2 < 100) %>% 
  filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio" | outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% 
  group_by(outcome.measure) %>% 
  mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>%
  group_by(sample.size, outcome.measure) %>% 
  summarize(median.effect = median(scaled.effect, na.rm = T)) %>% 
  ggplot(aes(x = sample.size, y = median.effect, group = outcome.measure)) + geom_point() + facet_wrap(~outcome.measure) + 
  theme_bw() + xlab("sample size") + ylab("median absolute normalized effect size")

#Plot separately for outcome measurres (mean effect size)
data %>% filter(total1 + total2 > 10 & total1 + total2 < 100) %>% 
  filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio" | outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% 
  group_by(outcome.measure) %>% 
  mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>%
  group_by(sample.size, outcome.measure) %>% 
  summarize(mean.effect = mean(scaled.effect, na.rm = T)) %>% 
  ggplot(aes(x = sample.size, y = mean.effect, group = outcome.measure)) + geom_point() + facet_wrap(~outcome.measure) + 
  theme_bw() + xlab("sample size") + ylab("mean absolute normalized effect size") + geom_smooth(se = F)

#Plot all data
data %>% filter(total1 + total2 > 10 & total1 + total2 < 100) %>% 
  filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio" | outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% 
  group_by(outcome.measure) %>% 
  mutate(sample.size = total1 + total2, scaled.effect = abs(scale(effect, center = T, scale = T))) %>%
  ggplot(aes(x = sample.size, y = scaled.effect, group = outcome.measure)) + geom_point() + facet_wrap(~outcome.measure) + 
  theme_bw() + xlab("sample size") + ylab("mean absolute normalized effect size") + geom_smooth(se = F)

#All data histogramms
data %>% filter(total1 + total2 > 10 & total1 + total2 < 100) %>% 
  filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio" | outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>% 
  group_by(outcome.measure) %>% 
  mutate(sample.size = total1 + total2, scaled.effect = scale(effect, center = T, scale = T)) %>%
  ggplot(aes(x = scaled.effect, group = outcome.measure)) + geom_histogram() + facet_wrap(~outcome.measure) + 
  theme_bw() + xlab("sample size") + ylab("mean absolute normalized effect size")

#Find out more about the outliers of std.mean differences as seen in plot before (around 100 outliers, likely missspecifications from pooling effect sizes)
data %>% filter(outcome.measure == "Std. Mean Difference") %>% 
  mutate(seffect = as.vector(scale(effect, center = T, scale = T)), size = ifelse(total1 + total2 < 27, 1, 0)) %>% filter(seffect > 10) %>% 
  ggplot(aes(x = seffect, fill = factor(size))) + geom_histogram()

#Check Risk ratios
data %>% filter(outcome.measure == "Risk Ratio") %>% 
  mutate(seffect = as.vector(scale(effect, center = T, scale = T)), size = ifelse(total1 + total2 < 27, 1, 0)) %>% #filter(seffect > 10) %>% 
  ggplot(aes(x = seffect, fill = factor(size))) + geom_histogram()

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

#Find a study with different comparisons:
data %>% group_by(file.nr, study.name) %>% distinct(comparison.name) %>% count %>% ggplot(aes(x = n)) + geom_histogram()

print(arrange(filter(data, file.nr == 21) %>% 
          select(comparison.name, study.name), study.name))


# #Check out metabin and metabias function with one meta-analysis.
# print(data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
#   filter(outcome.measure == "Odds Ratio" | outcome.measure == "Risk Ratio") %>% filter(events1 > 0 & events2 > 0) %>% 
#   mutate(counts = n()) %>% filter(counts > 9) %>% ungroup() %>% distinct(file.nr) %>% select(file.nr), n = 200)
# 
# 
# ds <- data %>% filter(file.nr == 11) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
#   filter(outcome.measure == "Odds Ratio" | outcome.measure == "Risk Ratio") %>% filter(events1 > 0 & events2 > 0) %>% 
#   mutate(counts = n()) %>% filter(counts == 18)
# print(ds %>% select(effect, events1, events2, study.name, outcome.name, counts), n = 200)
# 
# meta.ds <- metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, data = ds, method = "Inverse")
# trimfill.ds <- trimfill(meta.ds)
# funnel(trimfill.ds)
# 
# metabias(meta.ds, method = "peters")
# metabias(meta.ds, method = "peters")$statistic
# 
# cp <- copas(meta.ds)
# limitmeta(meta.ds)
# 
# 

# #Tasks:
# #Heterogeneity of outcomes (outcome.nr), possibly in relation with the total number of studies (study name)
# data %>% group_by(file.nr) %>% distinct(outcome.name) %>% count() %>% 
# 	filter(n < 50) %>% 
# 	full_join( data %>% group_by(file.nr) %>% distinct(outcome.name) %>% count %>% 
# 						 	filter(n > 49) %>% mutate(n = 50)) %>% 
# 	ggplot(aes(x = n)) + geom_histogram(col = "gray15", fill = "dodgerblue", bins = 50) +
# 	theme_bw() + labs(title = "Histogram of Research Subjects per Review") + xlab("Research subject number per review")
