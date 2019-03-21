load(file = "cochrane.RData")
load(file = "yodata.RData")

require(tidyverse)
require(biostatUZH)

#Look at Migraine sympton relieve review
arrange(filter(madata, file.name=="MYdwnld_CD005220StatsDataOnly_Version2.csv") %>% 
          select(file.nr, comparison.nr, comparison.name, outcome.nr, subgroup.nr, study.name), study.name)

arrange(filter(madata, file.nr == 21) %>% 
          select(file.nr, comparison.nr, outcome.nr, subgroup.nr, study.name), study.name)

arrange(filter(madata, file.nr == 21) %>% 
          select(file.nr, comparison.name, outcome.name, subgroup.name, study.name), study.name)


arrange(filter(madata, file.nr == 761) %>% 
          select(study.name, comparison.name, outcome.name), study.name)

arrange(filter(madata, file.nr == 2897) %>% 
          select(study.name, comparison.name, outcome.name), study.name)


#How many studies are in the dataset?
madata %>% group_by(file.nr) %>% distinct(study.name) %>% ungroup %>% count()

#Distribution of number of distinct studies per review (1: >=30):
madata %>% group_by(file.nr) %>% distinct(study.name) %>% count() %>% 
  filter(n < 30) %>% 
  full_join( madata %>% group_by(file.nr) %>% distinct(study.name) %>% count %>% 
               filter(n > 29) %>% mutate(n = 30)) %>% 
  group_by(n) %>% count() %>% 
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Histogram of Study Number per Review") + xlab("Study number per review")

#Number of different research subjects per review:
madata %>% group_by(file.nr) %>% distinct(comparison.name) %>% count() %>%
  filter(n < 50) %>%
  full_join( madata %>% group_by(file.nr) %>% distinct(comparison.name) %>% count %>%
               filter(n > 49) %>% mutate(n = 50)) %>%
  group_by(n) %>% count() %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Research Subject Number per Review") + xlab("Research subject number") + ylab("Number of reviews")

#Number of different outcomes for different research subjects summed up per review:
madata %>% group_by(file.nr, comparison.nr) %>% distinct(outcome.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% filter(n < 50) %>%
  full_join( madata %>% group_by(file.nr, comparison.nr) %>% distinct(outcome.name) %>% ungroup() %>% group_by(file.nr) %>% 
               count() %>% filter(n > 49) %>% mutate(n = 50)) %>%
  group_by(n) %>% count() %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Outcome Number per Review") + xlab("Outcome number") + ylab("Number of reviews")

#Mean and median number of different outcomes for different research subjects summed up per review:
madata %>% group_by(file.nr, comparison.nr) %>% distinct(outcome.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% ungroup %>% summarise(mean = mean(n))

madata %>% group_by(file.nr, comparison.nr) %>% distinct(outcome.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% ungroup %>% summarise(mean = median(n))

#Number of different subgroups given that they share research subject and outcome:
madata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% filter(n < 100) %>%
  full_join( madata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
               ungroup() %>% group_by(file.nr) %>% 
               count() %>% filter(n > 99) %>% mutate(n = 100)) %>%
  group_by(n) %>% count() %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Subgroup Number per Review") + xlab("Subgroup number") + ylab("Number of reviews")

madata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% 
  ggplot(aes(x = n)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Subgroup Number per Review") + xlab("Subgroup number") + ylab("Number of reviews")

madata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% filter(n < 100) %>%
  full_join( madata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
               ungroup() %>% group_by(file.nr) %>% 
               count() %>% filter(n > 99) %>% mutate(n = 100)) %>%
  ggplot(aes(x = n)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Subgroup Number per Review") + xlab("Subgroup number") + ylab("Number of reviews")

mean.diff.subg <- madata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% ungroup() %>% summarise(mean = mean(n))

median.diff.subg <- madata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% distinct(subgroup.name) %>% 
  ungroup() %>% group_by(file.nr) %>% 
  count() %>% ungroup() %>% summarise(median = median(n))

#Barbiturate Review:
arrange(filter(madata, file.nr == 21) %>% 
select(file.nr, comparison.name,  subgroup.name, outcome.name, study.name), study.name)

#Frequencies of outcome type among study results:
yodata %>% group_by(outcome.measure) %>% count() %>% arrange(desc(n)) %>% ungroup() %>% filter(row_number() < 9) %>% 
  ggplot(aes(x = outcome.measure, y = n)) +
  geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Outcome measure frequencies") + xlab("Frequency") + ylab("Outcome measure")

#Study year dsitribution:
madata %>% filter(study.year < 2019 & study.year > 1920) %>% ggplot(aes(x = study.year)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Study publication year frequency") + xlab("Year") + ylab("Frequency")

#The frequency of missing or erronous variables is also important. Especially, missing outcome measure indication, sample size indication, 
#standard deviation, mean and event indication for the cases where the latter are meaningful.
missing.mean <- madata %>% filter(!is.na(sd1) & sd1 > 0) %>% filter(!is.na(sd2) & sd2 > 0) %>% 
  summarise(total = length(sd1), na  = sum(is.na(mean1) | is.na(mean1))) %>% mutate(f = na/total) %>% select(f)

missing.sd <- madata %>% filter(!is.na(mean1) & !is.na(mean2)) %>% filter(mean1 != 0 & mean2 != 0) %>%
  summarise(total = length(mean1), na.sd  = sum(is.na(sd1) | is.na(sd1)), null.sd = length(sd1[sd1 <= 0 | sd2 <= 0])) %>% 
  mutate(f = (na.sd + null.sd)/total) %>% select(f)

missing.ssize <- madata %>% summarize(total = length(total1), null.ss = length(total1[total1 <= 0 | total2 <= 0])) %>% 
  mutate(f = null.ss / total) %>% select(f)

missing.events <- madata %>% filter(is.na(mean1) | mean1 == 0) %>% filter(is.na(mean2) | mean2 == 0) %>% 
  filter(is.na(sd1) | sd1 <= 0) %>% filter(is.na(sd2) | sd2 <= 0) %>% 
  summarise(null.ev = length(total1[events1 <= 0 | events2 <= 0]), total = length(events1)) %>% mutate(f = null.ev/total) %>% select(f)
  


#Tasks:
#Heterogeneity of outcomes (outcome.nr), possibly in relation with the total number of studies (study name)
madata %>% group_by(file.nr) %>% distinct(outcome.name) %>% count() %>% 
  filter(n < 50) %>% 
  full_join( madata %>% group_by(file.nr) %>% distinct(outcome.name) %>% count %>% 
               filter(n > 49) %>% mutate(n = 50)) %>% 
  ggplot(aes(x = n)) + geom_histogram(col = "gray15", fill = "dodgerblue", bins = 50) +
  theme_bw() + labs(title = "Histogram of Research Subjects per Review") + xlab("Research subject number per review")


#Frequencies of outcomes used
madata %>% group_by(outcome.measure) %>% 
  count %>% arrange(desc(n))

#Entries per review
madata %>% group_by(file.nr) %>% count

madata %>% group_by(file.nr) %>% count %>%  ggplot(aes(x = n)) + geom_histogram(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Entries per Review")

#Frequencies of unique combinations of review, compared research subject and outcome measure 
# -> Frequencies of trials analyzing the same subject in a review (with the same outcome)
madata %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count 

madata %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count %>% filter(n < 15) %>% 
  full_join( madata %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 14) %>% ungroup %>% summarise(n = 15, nn = sum(nn)))

madata %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count %>% 
  filter(n < 15) %>% 
  full_join( madata %>% group_by(file.nr, outcome.nr, comparison.nr) %>% count %>% group_by(n) %>% count %>% 
    filter(n > 14) %>% ungroup %>% summarise(n = 15, nn = sum(nn))) %>%  
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Trials comparing the same subject per Review")

#Frequencies of unique combinations of review, compared research subject, outcome measure and subgroup 
# -> Frequencies of trials analyzing the same subject with the same subgroups in a review
madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count 

madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% filter(n < 15) %>% 
  full_join( madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 14) %>% ungroup %>% summarise(n = 15, nn = sum(nn)))

madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% filter(n < 15) %>% 
  full_join( madata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 14) %>% ungroup %>% summarise(n = 15, nn = sum(nn))) %>%  
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Trials comparing the same subject in the same subgroups per Review")

#Frequencies of different outcomes per review
madata %>% group_by(file.nr) %>% summarise(diff.out = length(unique(outcome.name))) %>% count(diff.out)

madata %>% group_by(file.nr) %>% summarise(diff.out = length(unique(outcome.name))) %>% count(diff.out) %>%  
  ggplot(aes(x = diff.out, y = n)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs("Unique Comparisons per Review")

#Study year distribution:
madata %>% filter(study.year > 1900 & study.year < 2019) %>% ggplot(aes(x = study.year)) + geom_histogram()

#Study year distribution of distinct study years per review (i.e. duplicates are not considered)
madata %>% group_by(file.nr) %>% filter(study.year > 1900 & study.year < 2019) %>% distinct(study.year) %>% 
  ggplot(aes(x = study.year)) + 
  geom_histogram(col = "gray15", fill = "dodgerblue") + theme_bw() + labs(title = "Study years of analyses")

#Scaled effect sizes over time
madata %>% filter(file.nr < 5016) %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% mutate(scaled.effect = scale(effect), scaled.time = scale(study.year)) %>% 
  ggplot(aes(x = scaled.time, y = scaled.effect)) + geom_point(size = .005, alpha = 0.15) 

#Pvalues over time
temp.pval <- yodata %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
  distinct(study.year, .keep_all = T) %>% 
  mutate(counts = n()) %>% filter(counts > 1) %>% mutate(time.rank = rank(study.year)) %>% 
  mutate(pval = 2*(1-pnorm(abs(fishersz), mean = 0, sd = sqrt(variance))))  

temp.pval %>% filter(time.rank < 10) %>% ggplot(aes(x = time.rank, y = pval)) + 
  geom_jitter(col = "gray15", size = 0.3, alpha = .25) + 
  theme_bw() + xlab("Time rank") + ylab("P-value") + labs(title = "P-values over (ranked) Time")

# Cumulative number of total trials with
yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
  ungroup %>% arrange(desc(n)) %>% mutate(trials = cumsum(n * nn)) %>% 
  filter(n < 50) %>% 
  full_join( yodata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 49) %>% ungroup %>% summarise(n = 50, nn = sum(nn))) %>%  
  ggplot(aes(x = n, y = trials)) + geom_point(col = "gray15") + 
  theme_bw() + labs(title = "Sum of Trials with Number of Replicates >= n (ignoring subgroups") + ylab("Number of trials") + xlab("n, number of replicates")

yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>%
  ungroup %>% arrange(desc(n)) %>% mutate(trials = cumsum(n * nn)) %>% 
  filter(n < 50) %>% 
  full_join( yodata %>% group_by(file.nr, comparison.nr, outcome.nr, subgroup.nr) %>% count %>% group_by(n) %>% count %>% 
               filter(n > 49) %>% ungroup %>% summarise(n = 50, nn = sum(nn))) %>%  
  ggplot(aes(x = n, y = trials)) + geom_point(col = "gray15") + 
  theme_bw() + labs(title = "Sum of Trials with Number of Replicates >= n") + ylab("Number of trials") + xlab("n, number of replicates")

yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
  filter(n < 100) %>%
  full_join( yodata %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
               filter(n > 99) %>% ungroup %>% summarise(n = 100, nn = sum(nn))) %>%
  ggplot(aes(x = n, y = nn)) + geom_col(col = "gray15", fill = "dodgerblue") +
  theme_bw() + labs(title = "Sum of Reviews Comparing the Same Subject per Review")

yodata %>% filter(!is.na(fishersz)) %>% group_by(file.nr, comparison.nr, outcome.nr) %>% count %>% group_by(n) %>% count %>%
  ungroup %>% arrange(desc(nn)) %>% mutate(trials = cumsum(n * nn)) %>%
  ggplot(aes(x = nn, y = trials)) + geom_line(col = "gray15") +
  theme_bw() + labs(title = "Trials Comparing the Same Subject per Review")









  
  





  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
