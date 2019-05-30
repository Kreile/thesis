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



comp <- metac %>% select(meta.id, comparison.name) #%>% distinct(comparison.name)

comp$isin <- mapply(FUN = grepl, x = comp$comparison.name, pattern = "sensitivity")
length(which(comp$isin == TRUE))
comp$isin2 <- mapply(FUN = grepl, x = comp$comparison.name, pattern = "Sensitivity")
length(which(comp$isin2 == TRUE))

meta.wout.sensitivity <- metac %>% filter(!grepl(pattern = "sensitivity", comparison.name) | !grepl(pattern = "Sensitivity", comparison.name))


comp <- metac %>% select(meta.id, outcome.name) #%>% distinct(outcome.name)

comp$isin <- mapply(FUN = grepl, x = comp$outcome.name, pattern = "sensitivity")
length(which(comp$isin == TRUE))
comp$isin2 <- mapply(FUN = grepl, x = comp$outcome.name, pattern = "Sensitivity")
length(which(comp$isin2 == TRUE))

meta.wout.sensitivity <- metac %>% filter(!grepl(pattern = "sensitivity", outcome.name) | !grepl(pattern = "Sensitivity", outcome.name))


comp <- metac %>% select(meta.id, sungroup.name) #%>% distinct(sungroup.name)

comp$isin <- mapply(FUN = grepl, x = comp$sungroup.name, pattern = "sensitivity")
length(which(comp$isin == TRUE))
comp$isin2 <- mapply(FUN = grepl, x = comp$sungroup.name, pattern = "Sensitivity")
length(which(comp$isin2 == TRUE))

meta.wout.sensitivity <- metac %>% filter(!grepl(pattern = "sensitivity", sungroup.name) | !grepl(pattern = "Sensitivity", sungroup.name))
