# by Simon Schwab, 2019

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
  data$outcome.measure[grep("&#223;-agonist", data$outcome.measure)] = "beta-agonist"
  data$outcome.measure[grep("Hedge[s ]*&#180", data$outcome.measure)] = "Hedges' g"
  
  return(data)
}

# Cleans the raw data with some reg expressions
pb.clean = function(data) {
  
  names = unique(data$outcome.measure)
  data$outcome.measure.new = data$outcome.measure
  
  # Define aliases, first element will be used for all aliases
  aliases=list()
  aliases[[1]] = grep("risk ratio|RR|relat[ive]* risk", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[2]] = grep("^mean dif|^MD|change in", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[3]] = grep("Std. mean|SMD", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[4]] = sort(grep("% change|% increase|% rate|risk difference (%)", names, perl = TRUE, ignore.case = TRUE, value = TRUE))
  aliases[[5]] = grep("odds ratio", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[6]] = grep("hazard|^HR$|Survival HR|Change in HR|RR or HR|HR and variance", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[7]] = grep("rate ratio", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[8]] = grep("risk dif", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[9]] = grep("rate dif", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[10]] = grep("prevent[ed ]*fracti", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[11]] = grep("hedges", names, perl = TRUE, ignore.case = TRUE, value = TRUE)

  
  
  for (i in 1:length(aliases)) {
    idx =  data$outcome.measure.new %in% aliases[[i]] # find all aliases
    data$outcome.measure.new[idx] = aliases[[i]][1] # replace consistently with 1st element
  }
  
  return(list(data, aliases))
}

# Creates dataset used for pooling. We need to know what outcomes and subgroups can be combined
# in a pooled analysis. Searches for duplicates in {outcome, subgroup} pairs per m meta-analyses
# and c comparisons.
pb.pool = function(data) {
  
  r.id = unique(data$file.nr)
  mylist = list()
  count = 1
  for ( m in 1:length(r.id) ) {
    review = data[data$file.nr == r.id[m],]
    c.id = unique(review$comparison.nr)
    
    for ( c in 1:length(c.id) ) {
      compar = review[review$comparison.nr == c.id[c],]
      pairs = cbind(compar$outcome.nr, compar$subgroup.nr)
      idx = duplicated(pairs)
      
      if (any(idx)) {
        pairs = unique(array(pairs[idx,], dim=c(sum(idx),2)))
        
        for ( p in 1:nrow(pairs) ) {
          d = compar[compar$outcome.nr==pairs[p,1] & compar$subgroup.nr==pairs[p,2],]
          d$pool.nr = p
          mylist[[count]] = d
          count = count + 1
        }
      }
    }
  }
  return(rbindlist(mylist))
}

# Creates a database for reviews (each line one review).
pb.rev = function(data) {
  
  Nm = length(unique(data$file.nr))
  data.rev = data[!duplicated(data$file.nr), c(1, 3)]
  assert(nrow(data.rev) == Nm)
  rownames(data.rev) = 1:Nm
  
  data.rev$study.count = rep(NA, Nm)
  id = unique(data$file.nr)
  for (i in 1:Nm) {
    idx = data$file.nr == id[i]
    data.rev$study.count[i] = length(unique(data[idx,]$study.name))
  }
  return(data.rev)
}

pb.search = function(keyword) {
  idx1 = grepl(keyword, data$outcome.name)
  idx2 = grepl(keyword, data$comparison.name)
  idx3 = grepl(keyword, data$study.name)
  return(data[idx1|idx2|idx3,c(1,3,7,9,10,12,13,18,22)])
}

