#New Data process function

#Import data
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

if (file.exists(file.path(PATH_RESULTS, "data.RData"))) {
  load(file.path(PATH_RESULTS, "data.RData"))
} else {
  data = pb.readData(path = PATH_DATA, file = FILE)
  tmp = pb.clean(data)
  data = tmp[[1]]
  aliases = tmp[[2]]
  save(data, file =  file.path(PATH_RESULTS, "data.RData"))
}

if (file.exists(file.path(PATH_RESULTS, "data.ext.RData"))) {
  load(file.path(PATH_RESULTS, "data.ext.RData"))
} else {
  data.ext <- pb.process(data)
  save(data.ext, file =  file.path(PATH_RESULTS, "data.ext.RData"))
}

require(metafor)

data.ext2 <- pb.process2(data)

mcor <- function(data){
  if(all(data$outcome.type == "bin")){
    mt <- metacor(cor = cor.phi, n = total1 + total2, studlab = study.name, data = data)
  }else{
    if(all(data$outcome.type == "cont")){
      mt <- metacor(cor = cor.pearson, n = total1 + total2, studlab = study.name, data = data)
    } else{
      if(all(data$outcome.type == "surv")){
        mt <- metagen(TE = effect, seTE =  se, studlab = study.name, data = data)
      } else{
        mt <- NA
      }
    }
  }
  
  return(mt)
}

mcor <- function(data){
  if(all(outcome.type == "bin")){
    mt <- metacor(cor = cor.phi, n = total1 + total2, studlab = study.name)
  }else{
    if(all(outcome.type == "cont")){
      mt <- metacor(cor = cor.pearson, n = total1 + total2, studlab = study.name)
    } else{
      if(all(outcome.type == "surv")){
        mt <- metagen(TE = effect, seTE =  se, studlab = study.name)
      } else{
        mt <- NA
      }
    }
  }
  
  return(mt)
}

ex <- function(){
  if(all(outcome.type == "bin")){
    return("coool")
  } else{
    return("no")
  }
}

ex <- function(data){
  data %>% ifelse(outcome.type == "bin", "jo", "nai")
}

ex <- function(data){
  data %>% mutate(x = ifelse(outcome.type == "bin", "jo", "nai"))
}

ex <- function(outcome.type){
  return(ifelse(outcome.type == "bin", "jo", "nai"))
}

cor.just <- data.ext2 %>% group_by(meta.id) %>% mutate(n = n()) %>% filter(n > 9) %>% 
  filter(outcome.type != "rate")  ,
            est.ran = mcor(.)$TE.random,
            est.trim.fix = trimfill(mcor(.))$TE.fixed,
            est.reg = limitmeta(mcor(.))$TE.adjust,
            est.copas = auto.copas(mcor(.), sig.level = 0.1)[1],
            
            se.fix = mcor(.)$seTE.fixed,
            se.ran = mcor(.)$seTE.random,
            se.trim.fix = trimfill(mcor(.))$seTE.fixed,
            se.reg = limitmeta(mcor(.))$seTE.adjust,
            se.copas = auto.copas(mcor(.), sig.level = 0.1)[2])




#Examplary m.a. for test purpose: 134944 (bin), 158355 (cont)

bin <- data.ext2 %>% filter(meta.id == 134944)
cont <- data.ext2 %>% filter(meta.id == 158355)
bin2 <- data.ext2 %>% filter(meta.id == 16060)


forest(meta.id(134944))
forest(metacor(cor = cor.phi, n = total1 + total2, data = bin))
forest(metacor(cor = cor.fisher, n = total1 + total2, data = bin))


forest(rma.uni(ai = events1, bi = total1 - events1, ci = events2, di = total2 - events2,  
               measure = "RR", method = "FE", data = bin)) #For comparison
forest(rma.uni(yi = smd.pbit, vi = var.smd.pbit, measure = "PBIT", method = "FE", 
               data = bin)) 
forest(rma.uni(yi = cor.phi, vi = var.cor.phi, measure = "PHI", method = "FE", 
               data = bin)) 
forest(rma.uni(yi = cor.fisher, vi = var.cor.fisher, measure = "ZCOR", method = "FE", 
               data = bin)) 

forest(rma.uni(yi = effect, sei = se, measure = "MD", method = "FE", data = cont))
forest(rma.uni(m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, n1i = total1, n2i = total2, 
               measure = "MD", method = "FE", data = cont))
forest(rma.uni(yi = smd, vi = var.smd, measure = "SMD", method = "FE", data = cont))
forest(rma.uni(yi = cor.pearson, vi = var.cor.pearson, measure = "COR", method = "FE", data = cont))
forest(rma.uni(yi = cor.fisher, vi = var.cor.fisher, measure = "ZCOR", method = "FE", data = cont))
