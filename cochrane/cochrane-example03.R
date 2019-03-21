##########################################################
# read the generated (single) CSV file into a data frame:
madata <- read.csv("~/Data/cochrane/cochrane_2018-06-09.csv",
                   colClasses=c("integer", "character", "character",  # file       :  nr, name, doi
                                "integer", "integer",                 # file       :  index, version
                                "integer", "character",               # comparison :  nr, name
                                "integer", "character", "character",  # outcome    :  nr, name, measure
                                "integer", "character",               # subgroup   :  nr, name
                                "character", "integer",               # study      :  name, year
                                "numeric", "numeric",                 # effect, std.err.
                                "integer", "numeric",                 # group 1    :  events, total
                                "numeric", "numeric",                 # group 1    :  mean, se
                                "integer", "numeric",                 # group 2    :  events, total
                                "numeric", "numeric",                 # group 2    :  mean, se
                                "integer"))                           # total N
madata$total1 <- as.integer(madata$total1)
madata$total2 <- as.integer(madata$total2)


######################################
# generate & show some summary stats:
sumstats <- c("reviews"    =length(unique(madata$file.nr)),
              "comparisons"=length(unique(paste(madata$file.nr, madata$comparison.nr, sep="."))),
              "outcomes"   =length(unique(paste(madata$file.nr, madata$comparison.nr, madata$outcome.nr, sep="."))),
              "subgroups"  =length(unique(paste(madata$file.nr, madata$comparison.nr, madata$outcome.nr, madata$subgroup.nr, sep="."))),
              "studies"    =nrow(madata))
print(t(t(sumstats)))


##############################
# show details of data frame:
str(madata)

head(madata)

# check first entry, corresponding to first Cochrane review:
madata[madata$file.nr==1, ]

# example: show Richer-et-al review's data
# (see http://doi.org/10.1002/14651858.CD005220.pub2):
madata[madata$file.name=="MYdwnld_CD005220StatsDataOnly_Version2.csv", c(2,5,7,10,12)]

# show outcome measures and their frequencies:
t(t(sort(table(madata$outcome.measure), decreasing=TRUE)))

# non-zero entries for event counts only for RR/OR/RD/HR outcomes: 
t(t(sort(table(madata$outcome.measure[(madata$events1>0) | (madata$events2>0)]),dec=TRUE)))

# binary outcomes' names:
outcomenames <- unique(madata$outcome.measure[(madata$events1>0) | (madata$events2>0)])
# studies with binary outcomes:
binstudies <- is.element(madata$outcome.measure, outcomenames)
# fraction of binary outcomes:
mean(binstudies)  # 67%
# fraction of zero counts among binary outcomes:
mean((madata$events1[binstudies] == 0)
     | (madata$events1[binstudies] == madata$total1[binstudies])
     | (madata$events2[binstudies] == 0)
     | (madata$events2[binstudies] == madata$total2[binstudies]))
#  26 % (!)


##########################################
# distribution of analysis sizes;
# overall analyses:
tab1 <- table(table(paste(madata$file.nr, madata$comparison.nr, madata$outcome.nr, sep=".")))
tab1 <- c(tab1[1:15], ">15"=sum(tab1[-(1:15)]))
tab1 <- tab1/sum(tab1)
round(cbind("prob"=tab1, "cumulative"=cumsum(tab1)),3)
barplot(tab1, col="red3")

# subgroup analyses:
tab2 <- table(table(paste(madata$file.nr, madata$comparison.nr, madata$outcome.nr, madata$subgroup.nr, sep=".")))
tab2 <- c(tab2[1:15], ">15"=sum(tab2[-(1:15)]))
tab2 <- tab2/sum(tab2)
round(cbind("prob"=tab2, "cumulative"=cumsum(tab2)),3)
barplot(tab2, col="blue3")

# both analysis types:
barplot(rbind(tab1, tab2), beside=TRUE, col=c("red3","blue3"),
        xlab="number of studies included", ylab="frequency")
legend("topright", c("overall analyses","subgroup analyses"),
       pch=15, col=c("red3","blue3"))



