# #Transform Binary outcomes to odds ratios -> standardized mean difference -> correlation -> fisher scala correlation


data <- data %>% mutate(odds.ratio = NA,
               std.mean.d = NA,
               correlation = NA,
               fishersz = NA)

data[data$events1 > 0 & data$events2 > 0,] <- data[data$events1 > 0 & data$events2 > 0,] %>%
  mutate(odds.ratio = (events1/(total1 - events1)) / (events2/(total2 - events2)),
         std.mean.d = log(odds.ratio) * (sqrt(3)/pi),
         correlation = std.mean.d/sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2*total1)),
         fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), variance = 1/(total1 + total2 - 3))
#
# #Transform mean difference to standardized mean difference -> correlation -> fisher scala correlation
data[data$outcome.measure == "Mean Difference", ] <- data %>% 
  filter(outcome.measure == "Mean Difference") %>%
  mutate(std.mean.d = effect/se,
         correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
         fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
         variance = 1/(total1 + total2 - 3))
#
data[data$outcome.measure == "Std. Mean Difference", ] <- data %>% 
  filter(outcome.measure == "Std. Mean Difference") %>%
  mutate(std.mean.d = effect,
         correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
         fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
         variance = 1/(total1 + total2 - 3))
#
data[data$outcome.measure == "Mean difference", ] <- data %>% 
  filter(outcome.measure == "Mean difference") %>%
  mutate(std.mean.d = effect, 
         correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
         fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
         variance = 1/(total1 + total2 - 3))
#
data[data$outcome.measure == "mean difference", ] <- data %>% 
  filter(outcome.measure == "mean difference") %>%
  mutate(std.mean.d = effect,
         correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
         fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), variance = 1/(total1 + total2 - 3))

# #Filter out fishersz abs(correlation) > 1 and subgroups with fewer than 2 reproductions
data <- data %>% filter(fishersz < 1 & fishersz > -1) %>% filter(!is.na(fishersz))
data <- data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
  filter(fishersz < 1 & fishersz > -1) %>%
  mutate(counts = n()) %>% filter(counts > 1)

mly.fishersz <- function(data){
  data <- data %>% mutate(odds.ratio = NA,
                          std.mean.d = NA,
                          correlation = NA,
                          fishersz = NA)
  
  data[data$events1 > 0 & data$events2 > 0,] <- data[data$events1 > 0 & data$events2 > 0,] %>%
    mutate(odds.ratio = (events1/(total1 - events1)) / (events2/(total2 - events2)),
           std.mean.d = log(odds.ratio) * (sqrt(3)/pi),
           correlation = std.mean.d/sqrt((std.mean.d^2) + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), variance = 1/(total1 + total2 - 3))
  #
  # #Transform mean difference to standardized mean difference -> correlation -> fisher scala correlation
  data[data$outcome.measure == "Mean Difference", ] <- data %>% 
    filter(outcome.measure == "Mean Difference") %>%
    mutate(std.mean.d = effect/se,
           correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
           variance = 1/(total1 + total2 - 3))
  #
  data[data$outcome.measure == "Std. Mean Difference", ] <- data %>% 
    filter(outcome.measure == "Std. Mean Difference") %>%
    mutate(std.mean.d = effect,
           correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
           variance = 1/(total1 + total2 - 3))
  #
  data[data$outcome.measure == "Mean difference", ] <- data %>% 
    filter(outcome.measure == "Mean difference") %>%
    mutate(std.mean.d = effect, 
           correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), 
           variance = 1/(total1 + total2 - 3))
  #
  data[data$outcome.measure == "mean difference", ] <- data %>% 
    filter(outcome.measure == "mean difference") %>%
    mutate(std.mean.d = effect,
           correlation = std.mean.d/sqrt(std.mean.d^2 + ((total1 + total2)^2)/(total2*total1)),
           fishersz = 0.5 * log( (1 + correlation)/(1 - correlation) ), variance = 1/(total1 + total2 - 3))
  
  # #Filter out fishersz abs(correlation) > 1 and subgroups with fewer than 2 reproductions
  data <- data %>% filter(fishersz < 1 & fishersz > -1) %>% filter(!is.na(fishersz))
  data <- data %>% group_by(file.nr, outcome.nr, comparison.nr, subgroup.nr) %>%
    filter(fishersz < 1 & fishersz > -1) %>%
    mutate(counts = n()) %>% filter(counts > 1)
  return(data)
}