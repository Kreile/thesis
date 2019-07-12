# by Simon Schwab, 2019

library(testit)
require(tidyverse)
require(meta)
require(metasens)
require(biostatUZH)
require(gridExtra)
require(metafor)
require(xtable)

# Reads the data
pb.readData = function(path, file) {
  
  data = read.csv(file.path(path, file),
                  colClasses=c("integer", "character", "character",  # file       :  nr, name, doi
                               "integer", "integer",                 # file       :  index, version
                               "integer", "character",               # comparison :  nr, name
                               "integer", "character", "character",  # outcome    :  nr, name, measure
                               "integer", "character",               # subgroup   :  nr, name
                               "character", "integer",               # study      :  name, year
                               "numeric", "numeric",                 # effect, std.err.
                               "numeric", "numeric",                 # group 1    :  events, total
                               "numeric", "numeric",                 # group 1    :  mean, se
                               "numeric", "numeric",                 # group 2    :  events, total
                               "numeric", "numeric",                 # group 2    :  mean, se
                               "numeric"))                           # total N
  # some cleanup
  data$N = NULL # this colums has all NA's
  data$outcome.measure[grep("&#223;-agonist", data$outcome.measure)] = "beta-agonist"
  data$outcome.measure[grep("Hedge[s ]*&#180", data$outcome.measure)] = "Hedges' g"
  
  return(data)
}

# Cleans the raw data with some reg expressions
pb.clean = function(data) {
  
  names = unique(data$outcome.measure)
  data$outcome.measure.new = data$outcome.measure
  
  # Define aliases, first element will be used for all aliases
  # https://handbook-5-1.cochrane.org/chapter_9/9_2_2_5_what_is_the_event.htm
  aliases=list()
  aliases[[1]] = grep("risk ratio|^RR$|RR Ratio[s]*|IRR|relat[ive]* risk", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[2]] = grep("^mean dif[ference]*$|^mean dif[ference]*\\s+.[^%]|^MD", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[3]] = grep("change in", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[4]] = grep("Std. mean|SMD", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[5]] = grep("hedges", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[6]] = sort(grep("% change|% increase|% rate|risk difference \\(%\\)|mean difference \\[%\\]|changes in MVC \\[%]|changes \\[%]",
                           names, perl = TRUE, ignore.case = TRUE, value = TRUE))
  aliases[[7]] = grep("^odds ratio[s]*|[^peto] odds ratio", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[8]] = grep("peto", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[9]] = grep("hazard|^HR$|Survival HR|HR and variance", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[10]] = grep("rate ratio", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[11]] = grep("risk dif[ference]*$|risk dif[ference]*\\s+.[^%]", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[12]] = grep("^rate dif", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[13]] = grep("prevented fraction", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  
  
  
  for (i in 1:length(aliases)) {
    idx =  data$outcome.measure.new %in% aliases[[i]] # find all aliases
    data$outcome.measure.new[idx] = aliases[[i]][1] # replace consistently with 1st element
  }
  
  return(list(data, aliases))
}

# Adds pool.nr. We need to know what outcomes and subgroups can be combined in a pooled analysis.
# Searches for duplicates in {file.nr, comparison.nr, outcome.nr, subgroup.nr}.
pb.pool = function(data) {
  
  # create table to check for duplicates, these get the same pool nr
  tab = cbind(data$file.nr, data$comparison.nr, data$outcome.nr, data$subgroup.nr)
  idx = duplicated(tab)
  
  data$pool.nr = rep(NA, nrow(data))
  from = which(!idx)
  to   = c(from[2:length(from)] - 1, length(idx))
  
  # iterate through list of duplicates and assign nr
  c = 1
  for (i in 1:length(from)) {
    data$pool.nr[from[i]:to[i]] = c
    c = c + 1
  }
  return(data)
}

## Creates a database of reviews (each line one review).
## Fetches the title, year for each review
## Adds additional variables
library(roadoi)
pb.createReviews = function(data) {ls
  file_fetch = 'oadoi_fetch.RData'
  table = data.frame(file.nr = data$file.nr, doi=data$doi,
                     rev.title = NA, rev.year = NA)
  table$doi = tolower(as.character(table$doi))
  table = table[!duplicated(table),]
  
  # do not fetch if file already exists, takes ~30 min to get 5,000 titles
  if (file.exists(file.path(PATH_RESULTS, file_fetch))) {
    load(file.path(PATH_RESULTS, file_fetch))
  } else {
    fetch = oadoi_fetch(dois = table$doi, email = "simon.schwab@uzh.ch")
    save(fetch, file = file.path(PATH_RESULTS, file.fetch))
  }
  
  # merge fetched data with database
  idx = match(fetch$doi, table$doi)
  table$rev.title = rep(NA, nrow(table))
  table$rev.year = rep(NA, nrow(table))
  table$rev.title[idx] = fetch$title
  table$rev.year[idx] = fetch$year
  
  # populate review database with important variables
  
  # add number of studies per review and pool variables
  table$nr.studies = rep(NA, nrow(table))
  table$pool.nr1 = rep(NA, nrow(table)) 
  table$pool.outName1 = rep(NA, nrow(table))
  table$pool.count1 = rep(NA, nrow(table)) # nr of outcomes that can be pooled
  table$pool.nr2 = rep(NA, nrow(table))
  table$pool.outName2 = rep(NA, nrow(table)) 
  table$pool.count2 = rep(NA, nrow(table))
  
  for (i in 1:nrow(table)) {
    d = subset(data, data$file.nr == table$file.nr[i])
    table$nr.studies[i] = length(unique(d$study.name))
    s = sort(table(d$pool.nr), decreasing = TRUE) # count and sort pool.nr
    
    table$pool.nr1[i] = names(s[1])
    table$pool.nr2[i] = names(s[2])
    
    table$pool.count1[i] = s[1]
    table$pool.count2[i] = s[2]
    
    # get autcome name
    on = d$outcome.name[d$pool.nr == as.numeric(names(s[1]))]
    assert(length(unique(on)) == 1) # check all the same outcomes
    table$pool.outName1[i] = on[1]
    
    on = d$outcome.name[d$pool.nr == as.numeric(names(s[2]))]
    assert(length(unique(on)) == 1) # check all the same outcomes
    table$pool.outName2[i] = on[2]
  }
  
  return(table)
}

## Search keyword in database and returns all rows that match.
pb.search = function(keyword, data) {
  if ("rev.title" %in% names(data)) {
    idx1 = grepl(keyword, data$rev.title, ignore.case = T)
    return(data[idx1,c("file.nr", "doi", "rev.title", 
                       "rev.year","nr.studies","pool.count1","pool.count2")])
    
  } else {
    idx1 = grepl(keyword, data$outcome.name, ignore.case = T)
    idx2 = grepl(keyword, data$comparison.name, ignore.case = T)
    idx3 = grepl(keyword, data$study.name, ignore.case = T)
    idx4 = grepl(keyword, data$subgroup.name, ignore.case = T)
    return(data[idx1|idx2|idx3|idx4,c("file.nr","comparison.name","outcome.name","outcome.measure", 
                                 "subgroup.name","study.name","total1","total2")])
  }
  
}

########################################################################################################################
#Meta-analysis Functions

#Copas selection model automatic estimate and std. error and N.unpubl extraction:
#The estimate with smallest N.unpubl and a p-value larger than sig.level + 0.05 is chosen.
auto.copas <- function(meta.obj, sig.level){
  sig.level <- sig.level
  gamma0 <- -1.7 #analog to P(select|small trial w. sd = 0.4) = 0.1 and P(select|large trial w. sd  = 0.05) = 0.9
  gamma1 <- 0.16 #from limitmeta paper (Rücker 2011): "small range" procedure - if no nonsignificance - "broad range"
  copas <- copas(meta.obj, gamma0.range = c(gamma0, 2), gamma1.range = c(0, gamma1))
  pval.rsb <- copas$pval.rsb
  N.unpubl <- copas$N.unpubl
  if(all(pval.rsb < sig.level)){
    copas <- copas(meta.obj, , gamma0.range = c(2*gamma0 - 2, 2), gamma1.range = c(0, 2*gamma1))
    pval.rsb <- copas$pval.rsb
    N.unpubl <- copas$N.unpubl
    if(all(pval.rsb < sig.level)){
      corr.est <- NA
      se.corr.est <- NA
      N.unpubl <- NA
    } else{ 
      ind.nonsig <- which(pval.rsb > sig.level)
      ind.estimate <- which.min(N.unpubl[ind.nonsig])
      corr.est <- copas$TE.slope[ind.nonsig[ind.estimate]]
      se.corr.est <- copas$seTE.slope[ind.nonsig[ind.estimate]]
      N.unpubl <- copas$N.unpubl[ind.nonsig[ind.estimate]]
    }
  } else{
    if(all(pval.rsb > sig.level)){
      corr.est <- NA
      se.corr.est <- NA
      N.unpubl <- NA
    } else{
      ind.nonsig <- which(pval.rsb > sig.level)
      ind.estimate <- which.min(N.unpubl[ind.nonsig])
      corr.est <- copas$TE.slope[ind.nonsig[ind.estimate]]
      se.corr.est <- copas$seTE.slope[ind.nonsig[ind.estimate]]
      N.unpubl <- copas$N.unpubl[ind.nonsig[ind.estimate]]
    }
    return(c(corr.est, se.corr.est, N.unpubl))
  }
}


dupl.max.finder <- function(duplicate.index, study.number, metas){
  results <- rep(NA, length(duplicate.index))
  
  for(u in seq_along((duplicate.index))){
    
    duplicate.indices <- which(duplicate.index %in% duplicate.index[u])
    
    if(duplicate.index[u] != 0){
      
      meta.ids <- metas[duplicate.indices]
      
      max.meta.id <- meta.ids[which.max(study.number[duplicate.indices])]
      max.meta.indices <- which(metas %in% max.meta.id)
      if(length(unique(max.meta.id)) < 2){
        results[max.meta.indices] <- 0
      } else print("error")
      
    } else results[duplicate.indices] <- 0
    
  }
  
  results[is.na(results)] <- 1
  
  return(results)
}


dupl.finder <- function(effects, names, metas){
  
  results <- rep(NA, times = length(effects))
  meta.double.marker <- 1
  
  for(u in seq_along(effects)){
    
    if(is.na(results[u])){
      
      double.indices <- c(u, which(effects[u] == effects & names[u] == names))
      
      if(length(double.indices) > 1){
        
        meta.ids <- metas[double.indices]
        meta.indices <- which(metas %in% meta.ids)
        results[meta.indices] <- meta.double.marker
        meta.double.marker <- meta.double.marker + 1
        
      } else results[u] <- 0
      
    }}
  
  return(results)
}

tes.z <- function(data){
  
  if(all(data$outcome.type == "surv")){
    return(rep(NA, times = 7))
    
  } else{
    yi = data$z
    sei = sqrt(data$var.z)
  }
  

  alpha = 0.025
  to.ommit <- which(is.na(yi) | is.na(sei))
  
  if(length(to.ommit) > 0){
    yi <-  yi[-to.ommit]
    sei <- sei[-to.ommit]
  }
  
  side = ifelse(sum(pnorm(yi/sei) < .05) > sum(pnorm(yi/sei, lower.tail=FALSE) < .05), 
                "left", "right")
  
  ### FE meta-analysis for statistical power analysis
  est.fe <- rma(yi = yi, sei = sei, method = "FE")$b[1] 
  
  ### Compute statistical power and determine the number of observed 
  # statistically significant results
  if (side == "right") { 
    pow <- pnorm(qnorm(alpha, lower.tail = FALSE, sd = sei), mean = est.fe,
                 sd = sei, lower.tail = FALSE) #Probability to 
    O <- sum(pnorm(yi/sei, lower.tail = FALSE) < alpha)
  } else if (side == "left") {
    pow <- pnorm(qnorm(alpha, sd = sei), mean = est.fe, sd = sei)
    O <- sum(pnorm(yi/sei) < alpha)
  }
  
  E <- sum(pow) # Expected number of statistically significant result  
  n <- length(yi) # Number of studies in meta-analysis
  A <- (O - E)^2/E + (O - E)^2/(n - E) # Compute chi-square statistic
  pval.chi <- pchisq(A, 1, lower.tail = FALSE) # Compute p-value
  pval.chi <- ifelse(pval.chi < 0.5, pval.chi*2, (1-pval.chi)*2)
  
  est.p <- O/n
  A.wald <- ((est.p) - (E/n)) / sqrt((est.p)*(1-est.p)/n)
  pval.bin <- 1 - pnorm(A.wald)
  
  return(c(A = A, A.wald = A.wald, pval.chi = pval.chi, pval.bin = pval.bin, O = O, E = E, n = n))
}


tes.d <- function(data){
  if(all(data$outcome.type == "bin")){
    yi = data$smd.ordl
    sei = sqrt(data$var.smd.ordl)
  } else{
    if(all(data$outcome.type == "surv")){
      yi = data$effect
      sei = sqrt(data$se)
    } else{
      yi = data$smd
      sei = sqrt(data$var.smd)
    }
  }

  alpha = 0.05
  
  to.ommit <- which(is.na(yi) | is.na(sei))
  
  if(length(to.ommit) > 0){
    yi <-  yi[-to.ommit]
    sei <- sei[-to.ommit]
  }
  
  to.ommit <- which(sei == 0)
  
  if(length(to.ommit) > 0){
    yi <-  yi[-to.ommit]
    sei <- sei[-to.ommit]
  }
  
  ### FE meta-analysis for statistical power analysis
  est.fe <- rma(yi = yi, sei = sei, method = "FE")$b[1] 
  
  side = ifelse(sum(pnorm(yi/sei) < .05) > sum(pnorm(yi/sei, lower.tail=FALSE) < .05), 
                "left", "right")
  
  ### Compute statistical power and determine the number of observed 
  # statistically significant results
  if (side == "right") { 
    pow <- pnorm(qnorm(alpha, lower.tail = FALSE, sd = sei), mean = est.fe,
                 sd = sei, lower.tail = FALSE) #Probability to 
    O <- sum(pnorm(yi/sei, lower.tail = FALSE) < alpha)
  } else if (side == "left") {
    pow <- pnorm(qnorm(alpha, sd = sei), mean = est.fe, sd = sei)
    O <- sum(pnorm(yi/sei) < alpha)
  }
  
  E <- sum(pow) # Expected number of statistically significant result  
  n <- length(yi) # Number of studies in meta-analysis
  A <- (O - E)^2/E + (O - E)^2/(n - E) # Compute chi-square statistic
  pval.chi <- pchisq(A, 1, lower.tail = FALSE) # Compute p-value
  pval.chi <- ifelse(pval.chi < 0.5, pval.chi*2, (1-pval.chi)*2)
  
  est.p <- O/n
  A.wald <- ((est.p) - (E/n)) / sqrt((est.p)*(1-est.p)/n)
  pval.bin <- 1 - pnorm(A.wald)
  
  return(c(A = A, A.wald = A.wald, pval.chi = pval.chi, pval.bin = pval.bin, O = O, E = E, n = n))
}




pb.process2 <- function(data){
	data <- data %>% mutate(meta.id = group_indices(., file.nr, comparison.nr, outcome.nr, subgroup.nr)) %>%
		group_by(meta.id) %>% mutate(study.id = row_number()) %>% ungroup()
	
	data <- data %>% group_by(meta.id) %>%
	  mutate(n = n())  %>% ungroup() %>% group_by(file.nr) %>% 
	  mutate(dupl.id = dupl.finder(effects = effect, names = study.name, metas = meta.id), 
	         dupl.remove = dupl.max.finder(duplicate.index = dupl.id, study.number = n, metas = meta.id)) %>% ungroup()
	
	data <- data %>% mutate(study.year = ifelse(study.year < 2019, study.year, NA)) %>% 
		mutate(study.year = ifelse(study.year > 1920, study.year, NA)) 
	
	data <- data %>% mutate(
		outcome.type = ifelse(outcome.measure.new == "Hazard Ratio", "surv", 
													ifelse(outcome.measure.new == "Odds Ratio" | outcome.measure.new == "Risk Ratio" |
																 	outcome.measure.new == "Peto Odds Ratio" | outcome.measure.new == "Risk Difference", "bin", "cont")),
		outcome.type = ifelse(outcome.measure.new == "Rate Ratio" | outcome.measure.new == "Rate difference", "rate", outcome.type),
		lrr = NA,
		var.lrr = NA,
		smd = NA,
		var.smd = NA,
		smd.pbit = NA,
		var.smd.pbit = NA,
		smd.ordl = NA,
		var.smd.ordl = NA,
		cor.phi = NA,
		var.cor.phi = NA,
		pval.single = NA, 
		cor.pearson = NA,
		var.cor.pearson = NA,
		z = NA,
		var.z = NA,
		events1c = NA,
		events2c = NA,
		sig.single = NA)
	cont.ind <- which(data$outcome.type == "cont" & !is.na(data$outcome.type))
	
	data[cont.ind, "outcome.type"] <- data[cont.ind, ] %>% 
		mutate(outcome.type = ifelse(outcome.measure.new == "Std. Mean Difference"  | outcome.measure.new == "Mean Difference", 
																 outcome.type, NA)) %>% select(outcome.type)
	
	cont.ind <- which(data$outcome.type == "cont" & !is.na(data$outcome.type))
	
	bin.ind <- which(data$outcome.type == "bin" & !is.na(data$outcome.type))
	data[bin.ind,] <- escalc(data = data[bin.ind,], ai = events1, n2i = total2, ci = events2, n1i = total1, 
													 to = "only0",
													 add = 1/2,
													 measure = "RR", append = T, var.names = c("lrr", "var.lrr"))
	data[bin.ind,] <- escalc(data = data[bin.ind,], ai = events1, n2i = total2, ci = events2, n1i = total1, 
													 to = "only0",
													 add = 1/2,
													 measure = "PBIT", append = T, var.names = c("smd.pbit", "var.smd.pbit"))
	data[bin.ind,] <- escalc(data = data[bin.ind,], ai = events1, n2i = total2, ci = events2, n1i = total1, 
													 to = "only0",
													 add = 1/2,
													 measure = "OR2DL", append = T, var.names = c("smd.ordl", "var.smd.ordl"))
	data[bin.ind,] <- escalc(data = data[bin.ind,], ai = events1, n2i = total2, ci = events2, n1i = total1, 
													 to = "only0",
													 add = 1/2,
													 measure = "PHI", append = T, var.names = c("cor.phi", "var.cor.phi"))
	
	#What to do if there are 0/total cells:
	data[bin.ind,] <- data[bin.ind,] %>% mutate(events1c  = case_when(events1 == 0 ~ events1 + 0.5,
																																		events2 == 0 ~ events1 + 0.5,
																																		events1 - total1 == 0 ~ events1 - 0.5,
																																		events2 - total2 == 0 ~ events1 - 0.5,
																																		TRUE ~ events1),
																							events2c = case_when(events2 == 0 ~ events2 + 0.5,
																																	 events1 == 0 ~ events2 + 0.5,
																																	 events2 - total2 == 0 ~ events2 - 0.5,
																																	 events1 - total1 == 0 ~ events2 - 0.5,
																																	 TRUE ~ events2))
	
	#What to do if there is one zero and one total:
	data[bin.ind,] <- data[bin.ind,] %>% mutate(events1c  = case_when(events2 == 0 & events1 - total1 == 0 ~ events1,
																																		events1 == 0 & events2 - total2 == 0 ~ events1,
																																		TRUE ~ events1c),
																							events2c = case_when(events2 == 0 & events1 - total1 == 0 ~ events2,
																																	 events1 == 0 & events2 - total2 == 0 ~ events2,
																																	 TRUE ~ events2c))
	
	
	nomean1.ind <- which(data$outcome.type == "cont" & is.na(data$mean1))
	nomean2.ind <- which(data$outcome.type == "cont" & is.na(data$mean2))
	noeffects.ind <- which(data$outcome.type == "cont" & !is.na(data$effect))
	nomeans.ind <- intersect(nomean1.ind, nomean2.ind)
	nomeans.ind <- intersect(noeffects.ind, nomeans.ind) #no means, but a (std.) mean difference is provided.
	
	data[nomeans.ind,"mean1"] <- data[nomeans.ind,"effect"]
	data[nomeans.ind,"mean2"] <- 0
	
	data[cont.ind,] <- escalc(data = data[cont.ind,], m1i = mean1, m2i = mean2, 
														sd1i = sd1, sd2i = sd2, n1i = total1, n2i = total2,
														measure = "SMD" , append = T, var.names = c("smd", "var.smd"))
	
	inpute.given.ind <- which(is.na(data[cont.ind, "smd"]) & data[cont.ind, "outcome.measure.new"] == "Std. Mean Difference")
	data[inpute.given.ind, "smd"] <- data[inpute.given.ind, "effect"]
	data[inpute.given.ind, "var.smd"] <- data[inpute.given.ind, "se"]^2
	
	surv.ind <- which(data$outcome.type == "surv")
	data[surv.ind, "pval.single"] <- data[surv.ind, ] %>% mutate(pval.single = 2*(1-pnorm(abs((effect)/se)))) %>% select(pval.single)
	
	rate.ind <- which(data$outcome.type == "rate")
	data[rate.ind, "pval.single"] <- data[rate.ind, ] %>% mutate(pval.single = 2*(1-pnorm(abs((effect)/se)))) %>% select(pval.single)
	
	data[cont.ind, "pval.single"] <- data[cont.ind, ] %>% 
		mutate(t = (mean1 - mean2)/(sqrt((((total1-1)*sd1^2) + (total2-1)*sd2^2)/(total1 + total2 -2))*sqrt((1/total1)+(1/total2))), 
					 pval.single = 2*(1-pt(abs(t), df = total1 + total2 - 2))) %>% select(pval.single)
	
	sd1.0 <- which(data$outcome.type == "cont" & data$sd1 == 0) #Such that p-value here is not equal to zero
	sd2.0 <- which(data$outcome.type == "cont" & data$sd2 == 0)
	sd.0 <- union(sd1.0, sd2.0)
	data[sd.0, "pval.single"] <- NA
	
	data[bin.ind, "pval.single"] <- data[bin.ind, ] %>% 
	  mutate(pval.single = 2*(1-pnorm(abs(lrr/sqrt(var.lrr))))) %>% select(pval.single)
	
	data[cont.ind, "cor.pearson"]<- data[cont.ind, ] %>% mutate(
		a = ((total1 + total2)^2)/(total1*total2),
		cor.pearson  = smd/sqrt((smd^2) + a),
		var.cor.pearson = ((a^2)*var.smd)/(((smd^2)+a)^3)) %>% select(cor.pearson)
	
	data[cont.ind, "var.cor.pearson"]<- data[cont.ind, ] %>% mutate(
		a = ((total1 + total2)^2)/(total1*total2),
		cor.pearson  = smd/sqrt((smd^2) + a),
		var.cor.pearson = ((a^2)*var.smd)/(((smd^2)+a)^3)) %>% select(var.cor.pearson)
	
	data[bin.ind, "cor.pearson"]<- data[bin.ind, ] %>% mutate(
	  a = ((total1 + total2)^2)/(total1*total2),
	  cor.pearson  = smd.ordl/sqrt((smd.ordl^2) + a),
	  var.cor.pearson = ((a^2)*var.smd.ordl)/(((smd.ordl^2)+a)^3)) %>% select(cor.pearson)
	
	data[bin.ind, "var.cor.pearson"]<- data[bin.ind, ] %>% mutate(
	  a = ((total1 + total2)^2)/(total1*total2),
	  cor.pearson  = smd.ordl/sqrt((smd.ordl^2) + a),
	  var.cor.pearson = ((a^2)*var.smd.ordl)/(((smd.ordl^2)+a)^3)) %>% select(var.cor.pearson)
	
	
	pearson.ind <- which(!is.na(data$cor.pearson))
	
	data[pearson.ind, "z"] <- data[pearson.ind, ] %>% 
	  mutate(z = 0.5 * log( (1 + cor.pearson)/(1 - cor.pearson) , base = exp(1)), 
	         var.z= 1/(total1 + total2 - 3)) %>% select(z)
	data[pearson.ind, "var.z"] <- data[pearson.ind, ] %>% 
	  mutate(z = 0.5 * log( (1 + cor.pearson)/(1 - cor.pearson) , base = exp(1)), 
	         var.z= 1/(total1 + total2 - 3)) %>% select(var.z)
	
	
	# data[cont.ind, "z"] <- data[cont.ind, ] %>% 
	# 	mutate(z = 0.5 * log( (1 + cor.pearson)/(1 - cor.pearson) , base = exp(1)), 
	# 				 var.z= 1/(total1 + total2 - 3)) %>% 
	# 	select(z)
	# 
	# data[cont.ind, "var.z"] <- data[cont.ind, ] %>% 
	# 	mutate(z = 0.5 * log( (1 + cor.pearson)/(1 - cor.pearson) , base = exp(1)), 
	# 				 var.z= 1/(total1 + total2 - 3)) %>% 
	# 	select(var.z)
	# 
	# data[bin.ind, "z"] <- data[bin.ind, ] %>% 
	# 	mutate(z = 0.5 * log( (1 + cor.phi)/(1 - cor.phi) , base = exp(1)), 
	# 				 var.z= 1/(total1 + total2 - 3)) %>% 
	# 	select(z)
	# 
	# data[bin.ind, "var.z"] <- data[bin.ind, ] %>% 
	# 	mutate(z = 0.5 * log( (1 + cor.phi)/(1 - cor.phi) , base = exp(1)), 
	# 				 var.z= 1/(total1 + total2 - 3)) %>% 
	# 	select(var.z)
	
	data$sig.single <- ifelse(data$pval.single < 0.05, 1, 0)
	
	return(data)
}


#Function to classify how significance changed after correction
sig.change <- function(sig.before, sig.after){
  if(!is.na(sig.before & !is.na(sig.after))){
    if(sig.before == 0){
      results <- ifelse(sig.after == 0, 0, 1)
    } else{
      results <- ifelse(sig.after == 1, 2, 3)}
  } else{
    results <- NA
  }
  return(results)
}

meta.id <- function(meta.id.tolook){
  dt <- data.ext2 %>% filter(meta.id == meta.id.tolook)
  settings.meta(hakn=TRUE, method.tau="PM", method = "Inverse")
  if(any(is.na(dt$events1)) | all(dt$events1 == 0)){
    meta.ex <- metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name, data = dt)
  } else{
    meta.ex <- metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = "RR", data = dt, method = "Inverse")
  }
  return(meta.ex)
}

funnel.id <- function(meta.id.tolook){
  funnel(meta.id(meta.id.tolook))
}

trimfill.id <- function(meta.id.tolook){
  funnel(trimfill(meta.id(meta.id.tolook)))
}

limitmeta.id <- function(meta.id.tolook){
  funnel(limitmeta(meta.id(meta.id.tolook)))
}

copas.id <- function(meta.id.tolook){
  auto.copas(meta.id(meta.id.tolook), 0.1)
}

setting.id <- function(meta.id.tolook){
  print(data.ext2 %>% filter(meta.id == meta.id.tolook) %>% select(doi, comparison.name, comparison.nr, outcome.name, subgroup.name), n = 200)
}

est.id <- function(meta.id.tolook){
  return(meta %>% filter(meta.id == meta.id.tolook) %>% select(est.fixef, est.ranef, est.copas, est.trimfill.fixef, est.reg.ranef))
}

data.id <- function(meta.id.tolook){
  return(print(data.ext2 %>% filter(meta.id == meta.id.tolook) %>% 
                 select(study.name, outcome.measure, effect, lrr, smd, se, var.lrr, var.smd,  total1, total2, events1c, events2c, mean1, sd1, mean2, sd2) %>% 
                 arrange(se), n = 200))
}

data.cor.id <- function(meta.id.tolook){
  return(print(data.ext2 %>% filter(meta.id == meta.id.tolook) %>% 
                 select(outcome.measure.new, outcome.type, study.name, effect, se, total1, total2, events1, events1c, events2, events2c, 
                        mean1, mean2, cor.phi, cor.pearson, z, var.z) %>% 
                 arrange(se), n = 200))
}

tes.fct <- function(data){
  outcome.type <- unique(data$outcome.type)
  
  if(outcome.type == "bin"){
    yi <- data$lrr
    sei <- sqrt(data$var.lrr)
  }
  
  if(outcome.type == "cont"){
    yi <- data$smd
    sei <- sqrt(data$var.smd)
  }
  
  if(outcome.type == "surv"){
    yi <- data$effect
    sei <- data$se
  }
  
  
  alpha = 0.05
  
  to.ommit <- which(is.na(yi) | is.na(sei))
  
  if(length(to.ommit) > 0){
    yi <-  yi[-to.ommit]
    sei <- sei[-to.ommit]
  }
  
  to.ommit <- which(sei == 0)
  
  if(length(to.ommit) > 0){
    yi <-  yi[-to.ommit]
    sei <- sei[-to.ommit]
  }
  
  ### FE meta-analysis for statistical power analysis
  est.fe <- rma(yi = yi, sei = sei, method = "FE")$b[1] 
  
  side = ifelse(sum(pnorm(yi/sei) < .05) > sum(pnorm(yi/sei, lower.tail=FALSE) < .05), 
                "left", "right")
  
  ### Compute statistical power and determine the number of observed 
  # statistically significant results
  if (side == "right") { 
    pow <- pnorm(qnorm(alpha, lower.tail = FALSE, sd = sei), mean = est.fe,
                 sd = sei, lower.tail = FALSE) #Probability to 
    O <- sum(pnorm(yi/sei, lower.tail = FALSE) < alpha)
  } else if (side == "left") {
    pow <- pnorm(qnorm(alpha, sd = sei), mean = est.fe, sd = sei)
    O <- sum(pnorm(yi/sei) < alpha)
  }
  
  n <- length(yi) # Number of studies in meta-analysis
  E <- sum(pow) # Expected number of statistically significant result  
  
  A <- (O - E)^2/E + (O - E)^2/(n - E) # Compute chi-square statistic
  pval.chi <- pchisq(A, 1, lower.tail = FALSE) # Compute p-value
  pval.chi <- ifelse(pval.chi < 0.5, pval.chi*2, (1-pval.chi)*2) #Van Aert's test
  
  pval.bin <- pbinom(q = O-1, size = n, prob = E/n, lower.tail = F)
  
  return(c(A = A, Expected = E, pval.chi = pval.chi, pval.bin = pval.bin, O = O, E = E, n = n))
}

tes.fct.old <- function(data){
  outcome.type <- unique(data$outcome.type)
  
  if(outcome.type == "bin"){
    yi <- data$lrr
    sei <- sqrt(data$var.lrr)
  }
  
  if(outcome.type == "cont"){
    yi <- data$smd
    sei <- sqrt(data$var.smd)
  }
  
  if(outcome.type == "surv"){
    yi <- data$effect
    sei <- data$se
  }
  
  
  alpha = 0.05
  
  to.ommit <- which(is.na(yi) | is.na(sei))
  
  if(length(to.ommit) > 0){
    yi <-  yi[-to.ommit]
    sei <- sei[-to.ommit]
  }
  
  to.ommit <- which(sei == 0)
  
  if(length(to.ommit) > 0){
    yi <-  yi[-to.ommit]
    sei <- sei[-to.ommit]
  }
  
  ### FE meta-analysis for statistical power analysis
  est.fe <- rma(yi = yi, sei = sei, method = "FE")$b[1] 
  
  side = ifelse(sum(pnorm(yi/sei) < .05) > sum(pnorm(yi/sei, lower.tail=FALSE) < .05), 
                "left", "right")
  
  ### Compute statistical power and determine the number of observed 
  # statistically significant results
  if (side == "right") { 
    pow <- pnorm(qnorm(alpha, lower.tail = FALSE, sd = sei), mean = est.fe,
                 sd = sei, lower.tail = FALSE) #Probability to 
    O <- sum(pnorm(yi/sei, lower.tail = FALSE) < alpha)
  } else if (side == "left") {
    pow <- pnorm(qnorm(alpha, sd = sei), mean = est.fe, sd = sei)
    O <- sum(pnorm(yi/sei) < alpha)
  }
  
  n <- length(yi) # Number of studies in meta-analysis
  E <- sum(pow) # Expected number of statistically significant result  
  
  if(O == n){
    O <- O - 1/2
    E <- E - 1/2
  }
  
  A <- (O - E)^2/E + (O - E)^2/(n - E) # Compute chi-square statistic
  pval.chi <- pchisq(A, 1, lower.tail = FALSE) # Compute p-value
  pval.chi <- ifelse(pval.chi < 0.5, pval.chi*2, (1-pval.chi)*2)
  
  est.p <- O/n
  A.wald <- ((est.p) - (E/n)) / sqrt((est.p)*(1-est.p)/n)
  pval.bin <- 1 - pnorm(A.wald)
  
  return(c(A = A, A.wald = A.wald, pval.chi = pval.chi, pval.bin = pval.bin, O = O, E = E, n = n))
}




meta.fct.ranef <- function(data){
  outcome.type <- unique(data$outcome.type)
  
  if(outcome.type == "bin"){
    meta.result <- rma.uni(yi = lrr, sei = sqrt(var.lrr),
                           method = "PM", measure = "RR",  data = data, control = list(tau2.max = 300000))
  }
  
  if(outcome.type == "cont"){
    meta.result <- rma.uni(yi = smd, sei = sqrt(var.smd),
                           method = "PM", measure = "SMD",  data = data, control = list(tau2.max = 300000))
  }
  
  if(outcome.type == "surv"){
    meta.result <- rma.uni(yi = effect, sei = se, subset = which(se != 0),
                           method = "PM",  data = data, control = list(tau2.max = 300000))
  }
  
  return(meta.result)
}

metagen.bincont <- function(data){
  if (all(data$outcome.type == "bin")){
    meta <- metagen(TE = smd.ordl, seTE = sqrt(var.smd.ordl), studlab = study.name, data, sm = "SMD")
  } else{
    meta <- metagen(TE = smd, seTE = sqrt(var.smd), studlab = study.name, data, sm = "SMD")
  }
}


bias.side.fct <- function(outcome, lrr, var.lrr, smd, var.smd, effect, se){
  alpha = 0.05
  outcome <- unique(outcome)
  
  if(outcome == "bin"){
    yi <- lrr
    sei <- sqrt(var.lrr)
  }
  
  if(outcome == "cont"){
    yi <- smd
    sei <- sqrt(var.smd)
  }
  
  if(outcome == "surv"){
    yi <- effect
    sei <- se
  }
  
  to.ommit <- which(is.na(yi) | is.na(sei))
  
  if(length(to.ommit) > 0){
    yi <-  yi[-to.ommit]
    sei <- sei[-to.ommit]
  }
  
  to.ommit <- which(sei == 0)
  
  if(length(to.ommit) > 0){
    yi <-  yi[-to.ommit]
    sei <- sei[-to.ommit]
  }
  
  if(sum(pnorm(yi/sei) < .05) == sum(pnorm(yi/sei, lower.tail=FALSE) < .05)){
    side <- sign(est.fe <- rma(yi = yi, sei = sei, method = "FE")$b[1] )
  } else{
    side = ifelse(sum(pnorm(yi/sei) < .05) > sum(pnorm(yi/sei, lower.tail=FALSE) < .05), 
                  -1, 1)
  }
  return(side)
  
}


onesided.p <- function(stat, side, n, test.type){
  if(test.type == "reg"){
    if(side == -1){
      p <- pt(stat, df = n - 2)
    } else{
      p <- 1 - pt(stat, df = n -2)
    }
  } else if(side == -1){
    p <- pnorm(stat)
  } else{
    p <- 1 - pnorm(stat)
  }
  return(p)
}