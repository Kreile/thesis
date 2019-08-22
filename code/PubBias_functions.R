# by Simon Schwab, 2019

library(testit)
require(tidyverse)
require(meta)
require(metasens)
# require(biostatUZH)
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

# Reads the data
pb.readData2 = function(path, file) {
  
  data = read.csv(file.path(path, file),
                  colClasses=c("character",                          # id
                               "integer", "character", "character",  # comparison :  nr, name, id
                               "integer", "character", "character",  # outcome    :  nr, name, measure
                               "character", "character",             #               id, flag
                               "integer",  "character", "character", # subgroup   :  nr, name, id
                               "character","character","character",  # study      :  id, name, year
                               "character",                          #               data_source  
                               "numeric", "numeric",                 # effect, std.err.
                               "numeric", "numeric",                 # group 1    :  events, total
                               "numeric", "numeric",                 # group 1    :  mean, sd
                               "numeric", "numeric",                 # group 2    :  events, total
                               "numeric", "numeric"                  # group 2    :  mean, sd
                  ))
  
  # Fix invalid years from study.year
  idx_invalid = (!grepl("^(19|20)\\d{2}$", data$study.year, perl = TRUE))
  idx_hasyear = grepl(".*(19|20)(\\d{2}).*", data$study.year, perl = TRUE)
  fixme = data$study.year[idx_invalid & idx_hasyear]
  fixed = gsub(".*(19|20)(\\d{2}).*", "\\1\\2", fixme, perl = TRUE)
  data$study.year[idx_invalid & idx_hasyear] = fixed
  
  # Fix remaining invalid years from study.id
  idx_invalid = (!grepl("^(19|20)\\d{2}$", data$study.year, perl = TRUE))
  idx_hasyear = grepl(".*(19|20)(\\d{2}).*", data$study.id, perl = TRUE)
  fixme = data$study.id[idx_invalid & idx_hasyear]
  fixed = gsub(".*(19|20)(\\d{2}).*", "\\1\\2", fixme, perl = TRUE)
  data$study.year[idx_invalid & idx_hasyear] = fixed
  
  # set remaining to NA
  idx_invalid = (!grepl("^(19|20)\\d{2}$", data$study.year, perl = TRUE))
  data$study.year[idx_invalid] = NA
  data$study.year = as.numeric(data$study.year)
  
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
#library(roadoi)
pb.createReviews = function(data) {
  file_fetch = 'oadoi_fetch.RData'
  table = data.frame(file.nr = data$file.nr, doi=data$doi,
                     rev.title = NA, rev.year = NA)
  table$doi = as.character(table$doi)
  table = table[!duplicated(table),]
  
  # do not fetch if file already exists, takes ~30 min to get 5,000 titles
  if (file.exists(file.path(PATH_RESULTS, file_fetch))) {
    load(file.path(PATH_RESULTS, file_fetch))
  } else {
    fetch = oadoi_fetch(dois = table$doi, email = "simon.schwab@uzh.ch")
    save(fetch, file = file.path(PATH_RESULTS, file.fetch))
  }
  
  # merge fetched data with database
  idx = match(tolower(fetch$doi), tolower(table$doi))
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

library(rvest)
pb.crawlCochrane = function(dois) {
  
  mylist = list() # stores a list of data frames per review with all the studies
  
  # iterate through reviews
  for (m in 1:length(dois)) {
    
    print(paste("review nr:", m))
    url = paste0("https://www.cochranelibrary.com/cdsr/doi/", dois[m], "/references")
    src = read_html(url)
    #if (m %% 8 == 0) { Sys.sleep(40) } # or we get blacklisted for a while
    Sys.sleep(7)
    studies = html_nodes(css="div.bibliographies.references_includedStudies", x=src)
    n = length(studies)
    
    studies.df = data.frame(refShort = rep(NA, n), authors = rep(NA, n), title = rep(NA, n),
                            journal = rep(NA, n), year = rep(NA, n), doi = rep(NA, n),
                            timesCited = rep(NA, n), refLong = rep(NA, n), error = rep(NA, n))
    
    for (i in 1:n) { # iterate through studies of a meta-analysis
      
      #print(i)
      refShort = html_text(html_nodes(css="h4.title", x=studies[i]))
      refShort = sub("\\s[\\{\\+].*", "", refShort)
      # often multiple references are listed for each study called records
      records = html_nodes(css="div.bibliography-section", x=studies[i])
      
      idx_rec = 1 # take fist record if the is only one
      
      # if there are more than one
      if (length(records) > 1) {
        # select correct record and remove all special characters like Malmström -> Malmstr.m and 2018a -> 2018
        idx_rec = grep(paste0("^\\s*[the]*\\s*", gsub("[^A-Za-z0-9]+|[a-z]$", ".*", refShort)),
                       html_text(records), ignore.case = TRUE)
      }
      
      # if none has been found
      if (length(idx_rec) == 0) {  # relax on year only (due Country codes like Chile, Europe as First Author name)
        idx_rec = grep(sub("(.*)(19[0-9]{2}|20[0-9]{2})(.*)", "\\2", refShort, perl = TRUE), # match year 19** or 20**
                       html_text(records), ignore.case = TRUE)
      }
      
      # if there are still more than 1...
      if (length(idx_rec) > 1) { # if more than one record, take the one that is cited more often, or has more links to original work
        moreCites = rep(NA, length(idx_rec))
        moreLinks = rep(NA, length(idx_rec))
        for (j in 1:length(idx_rec)) { # for each duplicate with same shortRef author/year get citations and links
          links = html_nodes("a", x = records[idx_rec[j]])
          moreLinks[j] = length(links)
          idx = grep("Web of Science", html_text(links), ignore.case = TRUE) # determine link times cited
          timesCited = s(html_text(html_nodes(css="span", x=links[idx])))
          moreCites[j] = as.numeric(sub(".*:", "", timesCited))
        }
        if (all(is.na(moreCites))) {
          idx_rec = idx_rec[which.max(moreLinks)]
        } else {
          idx_rec = idx_rec[which.max(moreCites)]
        }
      }
      
      # if idx_rec is still not 1 then something went wrong.
      if (length(idx_rec) != 1) {
        authors = title = journal = year = doi = timesCited = refLong = NA
        error = 1
        
      } else { # success: fetch and fill data
        
        refLong  = html_text(html_nodes(css="div", x=records[idx_rec]))
        authors  = sub("\\..*", "", refLong)
        title    = s(html_text(html_nodes(css="span.citation-title", x=records[idx_rec])))
        journal  = gsub("^\\s|\\s$", "", s(html_text(html_nodes(css="span.citation", x=records[idx_rec]))))
        year     = s(html_text(html_nodes(css="span.pubYear", x=records[idx_rec])))
        
        links    = html_nodes("a", x = records[idx_rec]) # select
        idx = grep("Link to article", html_text(links), ignore.case = TRUE) # determine link doi
        doi = s(html_attr(name="href", x=links[idx]))
        
        idx = grep("Web of Science", html_text(links), ignore.case = TRUE) # determine link times cited
        timesCited = s(html_text(html_nodes(css="span", x=links[idx])))
        timesCited = sub(".*:", "", timesCited)
        error = 0
        
      }
      # n x 9 data frame
      studies.df[i,] = c(refShort, authors, title, journal, year, doi, timesCited, refLong, error)
    }
    
    mylist[[m]] = studies.df
  }
  
  return(mylist)
}




# Returns NA if string has length 0
s <- function(object) {
  if (length(object) == 0) {
    return(NA)
  } else {
    return(object)
  }
}


#Random effects meta-analysis function:
meta.fct.ranef.mom <- function(data){
  outcome.flag <- unique(data$outcome.flag)
  outcome.measure.merged <- unique(data$outcome.measure.merged)
  
  if(outcome.flag == "DICH"){
    meta.result <- rma.uni(yi = lrr, sei = sqrt(var.lrr),
                           method = "DL", measure = "RR",  data = data)
  }
  
  if(outcome.flag == "CONT"){
    if(outcome.measure.merged == "SMD"){
      meta.result <- rma.uni(yi = cohensd, sei = sqrt(var.cohensd),
                             method = "DL", measure = "SMD",  data = data)
    }
    else{
    data <- data %>% filter(se != 0)
    meta.result <- rma.uni(yi = effect, sei = se,
                             method = "DL", data = data)
    }
  }
  
  if(outcome.flag == "IV"){
    if(outcome.measure.merged == "SMD" | outcome.measure.merged == "MD"){
      meta.result <- rma.uni(yi = effect, sei = se,
                             method = "DL", measure = outcome.measure.merged,  
                             data = data[data$se != 0,])
    } else {
      if(outcome.measure.merged == "RR" | outcome.measure.merged == "OR"){
        meta.result <- rma.uni(yi = log(effect), sei = se,
                               method = "DL", measure = outcome.measure.merged,  
                               data = data[data$effect != 0 & data$se != 0,])
      } else{
        if(outcome.measure.merged == "Hazard Ratio"){
          meta.result <- rma.uni(yi = log(effect), sei = se,
                                 method = "DL",  
                                 data = data[data$effect != 0 & data$se != 0,])
        } else{
          if(outcome.measure.merged == "Rate Ratio"){
            meta.result <- rma.uni(yi = log(effect), sei = se,
                                   method = "DL", measure = "IRR",  
                                   data = data[data$effect != 0 & data$se != 0,])
          } else{
            yi = NA
            sei = NA
          }
        }
      }
    }
  }
  
  return(meta.result)
}

#Meta-analysis for IV outcomes:
metagen.iv <- function(data){
  if(all(data$outcome.measure.merged == "SMD") | all(data$outcome.measure.merged == "MD")){
    meta <- metagen(TE = effect, seTE = se, studlab = study.name, sm = unique(data$outcome.measure.merged),
                    data = data)
  } else {
    if(all(data$outcome.measure.merged == "RR") | all(data$outcome.measure.merged == "OR")){
      meta <- metagen(TE = log(effect), seTE = se, studlab = study.name, sm = unique(data$outcome.measure.merged),
                      data = data)
    } else{
      if(all(data$outcome.measure.merged == "Hazard Ratio")){
        meta <- metagen(TE = log(effect), seTE = se, studlab = study.name, sm = "HR",
                        data = data)
      } else{
        if(all(data$outcome.measure.merged == "Rate Ratio")){
          meta <- metagen(TE = log(effect), seTE = se, studlab = study.name, 
                          data = data)
        } else{
          meta <- NA
        }
      }
    }
  }
  return(meta)
}

#Meta-analysis for standardized mean difference
metagen.bincont2 <- function(data){
  if (all(data$outcome.flag == "DICH")){
    meta <- metagen(TE = smd.ordl, seTE = sqrt(var.smd.ordl), studlab = study.name, data, sm = "SMD")
  } else{
    if(all(data$outcome.flag == "CONT")){
      meta <- metagen(TE = cohensd, seTE = sqrt(var.cohensd), studlab = study.name, data, sm = "SMD")
    } else{
      meta <- metagen(TE = effect, seTE = se, studlab = study.name, data, sm = "SMD")
    }
  }
}

#Get sign of expected bias side direction:
bias.side.fct2 <- function(outcome, outcome.measure.merged, lrr, var.lrr, smd, var.smd, effect, se){
  alpha = 0.05
  outcome <- unique(outcome)
  
  if(outcome == "DICH"){
    yi <- lrr
    sei <- sqrt(var.lrr)
  }
  
  if(outcome == "CONT"){
    if(all(outcome.measure.merged == "MD")){
      yi = effect
      sei = se
    } else{
      yi <- smd
      sei <- sqrt(var.smd)
    }
  }
  
  if(outcome == "IV"){
    if(outcome.measure.merged == "SMD" | outcome.measure.merged == "MD"){
      yi = effect
      sei = se
    } else {
      if(outcome.measure.merged == "RR" | outcome.measure.merged == "OR"){
        yi = log(effect)
        sei = se
      } else{
        if(outcome.measure.merged == "Hazard Ratio"){
          yi = log(effect)
          sei = se
        } else{
          if(outcome.measure.merged == "Rate Ratio"){
            yi = log(effect)
            sei = se
          } else{
            yi = NA
            sei = NA
          }
        }
      }
    }
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

#Function to get one-sided p-value from publication bias test:
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

#Excess significance test function:
tes.fct2 <- function(data){
  outcome.flag <- unique(data$outcome.flag)
  outcome <- unique(data$outcome.measure.merged)
  
  if(outcome.flag == "DICH"){
    yi <- data$lrr
    sei <- sqrt(data$var.lrr)
  }
  
  if(outcome.flag == "CONT"){
    if(outcome == "MD"){
      yi = data$effect
      sei = data$se
    } else{
      yi <- data$cohensd
      sei <- sqrt(data$var.cohensd)
    }
  }
  
  if(outcome.flag == "IV"){
    if(outcome == "SMD" | outcome == "MD"){
      yi = data$effect
      sei = data$se
    } else {
      if(outcome == "RR" | outcome == "OR"){
        yi = log(data$effect)
        sei = data$se
      } else{
        if(outcome == "Hazard Ratio"){
          yi = log(data$effect)
          sei = data$se
        } else{
          if(outcome == "Rate Ratio"){
            yi = log(data$effect)
            sei = data$se
          } else{
            yi = NA
            sei = NA
          }
        }
      }
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
  
  n <- length(yi) # Number of studies in meta-analysis
  E <- sum(pow) # Expected number of statistically significant result  
  
  A <- (O - E)^2/E + (O - E)^2/(n - E) # Compute chi-square statistic
  pval.chi <- pchisq(A, 1, lower.tail = FALSE) # Compute p-value
  pval.chi <- ifelse(pval.chi < 0.5, pval.chi*2, (1-pval.chi)*2) #Van Aert's test
  
  pval.bin <- pbinom(q = O-1, size = n, prob = E/n, lower.tail = F)
  
  return(c(A = A, Expected = E, pval.chi = pval.chi, pval.bin = pval.bin, O = O, E = E, n = n))
}
#----------------------------------------------------------------------------------------------#

#Dataset processing function to get transformed effect sizes, p-values, event counts with increments, etc. :
pb.process3 <- function(data){
  data <- data %>% mutate(meta.id = group_indices(., id, comparison.id, outcome.id, subgroup.id)) %>%
    group_by(meta.id) %>% mutate(study.id2 = row_number()) %>% ungroup()
  
  #Mark duplicate effects between meta-analyses:
  data <- data %>% group_by(meta.id) %>%
    mutate(n = n())  %>% ungroup() %>% group_by(id) %>% 
    mutate(dupl.id = dupl.finder(effects = effect, names = study.name, metas = meta.id), 
           dupl.remove = dupl.max.finder(duplicate.index = dupl.id, study.number = n, metas = meta.id)) %>% ungroup()
  
  data <- data %>% mutate(lrr = NA,
    var.lrr = NA,
    hedgesg = NA,
    var.hedgesg = NA,
    cohensd = NA,
    var.cohensd = NA,
    smd.ordl = NA,
    var.smd.ordl = NA,
    pval.single = NA, 
    cor.pearson = NA,
    var.cor.pearson = NA,
    z = NA,
    var.z = NA,
    events1c = NA,
    events2c = NA,
    sig.single = NA,
    smd.pool = NA,
    se.smd.pool = NA)

  
  cont.ind <- which(data$outcome.flag == "CONT")
  bin.ind <- which(data$outcome.flag == "DICH")
  IV.ind <- which(data$outcome.flag == "IV")
  #----------------------------------------------------------------------------------------------#
  
  #Outcome.flag == "DICH"
  data[bin.ind,] <- escalc(data = data[bin.ind,], ai = events1, n2i = total2, ci = events2, n1i = total1, 
                           to = "only0",
                           add = 1/2,
                           measure = "RR", append = T, var.names = c("lrr", "var.lrr")) #Calculate Risk Ratios

  data[bin.ind,] <- escalc(data = data[bin.ind,], ai = events1, n2i = total2, ci = events2, n1i = total1, 
                           to = "only0",
                           add = 1/2,
                           measure = "OR2DL", append = T, var.names = c("smd.ordl", "var.smd.ordl")) #Calculate SMD (w. logistic transf.)

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
  data[bin.ind, "pval.single"] <- data[bin.ind, ] %>% 
    mutate(pval.single = 2*(1-pnorm(abs(lrr/sqrt(var.lrr))))) %>% select(pval.single)
  #----------------------------------------------------------------------------------------------#
  
  #Outcome.flag == "CONT"
  data[cont.ind,] <- escalc(data = data[cont.ind,], m1i = mean1, m2i = mean2, 
                            sd1i = sd1, sd2i = sd2, n1i = total1, n2i = total2,
                            measure = "SMD" , append = T, var.names = c("hedgesg", "var.hedgesg")) #Calculate Hedge's g
  
  data[cont.ind, ] <- data[cont.ind, ] %>% 
    mutate(cohensd = (mean1 - mean2)/sqrt((((total1 - 1)*sd1^2) + (total2 - 1)*sd2^2)/(total1 + total2 - 2)),
           var.cohensd = ((total1 + total2)/(total1 * total2)) + (cohensd^2)/(2*(total1 + total2)))

  data[cont.ind, ] <- data[cont.ind, ] %>% mutate(cohensd = case_when(sd1 == 0 ~ NA_real_,
                                                                      sd2 == 0 ~ NA_real_,
                                                                      total2 == 0 ~ NA_real_,
                                                                      total1 == 0 ~ NA_real_,
                                                                      TRUE ~ cohensd))
  
  #p-value for continuous outcomes (Student):
  data[cont.ind, "pval.single"] <- data[cont.ind, ] %>% 
    mutate(t = (mean1 - mean2)/(sqrt((((total1-1)*sd1^2) + (total2-1)*sd2^2)/(total1 + total2 -2))*sqrt((1/total1)+(1/total2))), 
           pval.single = 2*(1-pt(abs(t), df = total1 + total2 - 2))) %>% select(pval.single)
  
  sd1.0 <- which(data$outcome.flag == "CONT" & data$sd1 == 0) 
  sd2.0 <- which(data$outcome.flag == "CONT" & data$sd2 == 0)
  sd.0 <- union(sd1.0, sd2.0)
  data[sd.0, "pval.single"] <- NA #Such that p-value here is not equal to zero
  #----------------------------------------------------------------------------------------------#
  
  #Outcome.flag == "IV"
  data[IV.ind, "pval.single"] <- data[IV.ind, ] %>% mutate(pval.single = 2*(1-pnorm(abs((effect)/se)))) %>% 
    select(pval.single)
  #----------------------------------------------------------------------------------------------#
  
  #Declare smd.ordl, cohensd and SMD of IV flag all as "smd.pool":
  data <- data %>% mutate(smd.pool = case_when(outcome.flag == "DICH" ~ smd.ordl,
                                               outcome.flag == "CONT" ~ cohensd,
                                               outcome.flag == "IV" & outcome.measure.merged == "SMD" ~ effect,
                                               TRUE ~ NA_real_),
                          se.smd.pool = case_when(outcome.flag == "DICH" ~ sqrt(var.smd.ordl),
                                                  outcome.flag == "CONT" ~ sqrt(var.cohensd),
                                                  outcome.flag == "IV" & outcome.measure.merged == "SMD" ~ se,
                                                  TRUE ~ NA_real_))
  #----------------------------------------------------------------------------------------------#
  
  #Transform effect sizes:
  smd.pool.ind <- which(!is.na(data$smd.pool))
  
  #Pearson correlation:
  data[smd.pool.ind, "cor.pearson"]<- data[smd.pool.ind, ] %>% mutate(
    a = ((total1 + total2)^2)/(total1*total2),
    cor.pearson  = smd.pool/sqrt((smd.pool^2) + a),
    var.cor.pearson = ((a^2)*se.smd.pool^2)/(((smd.pool^2)+a)^3)) %>% select(cor.pearson)
  
  data[smd.pool.ind, "var.cor.pearson"]<- data[smd.pool.ind, ] %>% mutate(
    a = ((total1 + total2)^2)/(total1*total2),
    cor.pearson  = smd.pool/sqrt((smd.pool^2) + a),
    var.cor.pearson = ((a^2)*se.smd.pool^2)/(((smd.pool^2)+a)^3)) %>% select(var.cor.pearson)
  
  #Fisher's z-score:
  data[smd.pool.ind, "z"] <- data[smd.pool.ind, ] %>% 
    mutate(z = 1/2 * log((1+cor.pearson)/(1-cor.pearson)),
           var.z = 1/(total1 + total2 - 3)) %>% select(z)
  data[smd.pool.ind, "var.z"] <- data[smd.pool.ind, ] %>% 
    mutate(z = 1/2 * log((1+cor.pearson)/(1-cor.pearson)), 
           var.z = 1/(total1 + total2 - 3)) %>% select(var.z)
  
  #----------------------------------------------------------------------------------------------#
  data$sig.single <- ifelse(data$pval.single < 0.05, 1, 0)
  
  return(data)
}



#Function to use with "data %>% mutate(.. = dupl.finder(..))"
dupl.finder <- function(effects, names, metas){
  
  results <- rep(NA, times = length(effects))
  meta.double.marker <- 1
  
  
  for(u in seq_along(effects)){
    
    if(is.na(results[u])){
      
      double.indices <- which(effects[u] == effects & names[u] == names)
      
      if(length(double.indices) > 1){
        
        if(any(!is.na(results[double.indices]))){ #If any meta-analysis has already duplicates, give all the same double marker
          
          already.marked <- which(!is.na(results[double.indices]))
          meta.double.marker.2 <- unique(results[double.indices[already.marked]])
          meta.ids.2 <- unique(metas[which(results %in% meta.double.marker.2)]) #Collect meta.ids with same double markers
          meta.ids <- metas[double.indices] #Collect metas with the same results
          
          meta.ids.pooled <- union(meta.ids, meta.ids.2)
          meta.indices <- which(metas %in% meta.ids.pooled)
          results[meta.indices] <- meta.double.marker
          
        } else{
          
          meta.ids <- metas[double.indices]
          meta.indices <- which(metas %in% meta.ids)
          results[meta.indices] <- meta.double.marker
          
        }
      } else results[u] <- 0
      
      meta.double.marker <- meta.double.marker + 1
    }}
  return(results)
}


#Find the biggest meta-analysis within a duplicate set, again as "data %>% mutate(.. = dupl.finder(..))":
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
    copas <- copas(meta.obj, gamma0.range = c(2*gamma0 - 2, 2), gamma1.range = c(0, 2*gamma1))
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
  }
  result <- c(est = corr.est, se =  se.corr.est, missing =  N.unpubl)
  return(result)
}

#Functions to quickly scrutinize data with meta.ids:
setting.id <- function(meta.id.tolook) { 
  x <- meta.f %>% filter(meta.id == meta.id.tolook) %>% 
    select(id, comparison.name, outcome.name, subgroup.name)
  return(x)}

getdata <- function(meta.id.tolook) return(data.ext2 %>% filter(meta.id == meta.id.tolook))

getinfo <- function(meta.id.tolook, np) print(getdata(meta.id.tolook) %>% 
                                                select(outcome.flag, effect, se, z, var.z))

meta.id <- function(meta.id.tolook){
  dt <- getdata(meta.id.tolook)
  settings.meta(hakn=FALSE, method.tau="PM", method = "Inverse")
  if(unique(dt$outcome.flag) == "DICH"){
    meta.ex <- metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2, studlab = study.name, sm = "RR", data = dt, method = "Inverse")
  } else{
    if(unique(dt$outcome.flag) == "CONT"){
      meta.ex <-     metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name, data = dt)
} else{
      meta.ex <- metagen(TE = effect, seTE = se, data = dt)
    }
  }
  return(meta.ex)
}

funnel.id <- function(meta.id.tolook){
  funnel(meta.id(meta.id.tolook))
}


limitmeta.id <- function(meta.id.tolook){
  funnel(limitmeta(meta.id(meta.id.tolook)))
}

est.id <- function(meta.id.tolook){
  return(meta.f %>% filter(meta.id == meta.id.tolook) %>% select(est.fixef, est.ranef, est.copas, est.reg))
}
