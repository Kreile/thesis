dt <- data.ext2 %>% filter(meta.id == meta.bin$meta.id[33])
forest.meta(metacor(cor = z, n = total1 + total2, studlab = study.name, data = dt, sm = "ZCOR"))
forest.meta(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2,
                    studlab = study.name, sm = "RR",
                    method = "Inverse", data = dt))


dt <- data.ext2 %>% filter(meta.id == meta.cont[which(meta.cont$outcome.measure.merged == "MD")[33],]$meta.id)

forest.meta(metacor(cor = z, n = total1 + total2, studlab = study.name, data = dt, sm = "ZCOR"))
forest.meta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2,
                     mean.c = mean2, sd.c = sd2, sm = "MD", studlab = study.name,
                     data = dt))
forest.meta(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2,
                     mean.c = mean2, sd.c = sd2, sm = "SMD", studlab = study.name,
                     data = dt))

par(mfrow= c(1,2))
forest(rma.uni(yi = z, vi = var.z, measure = "ZCOR", slab = study.name, data = dt))
forest(rma.uni(yi = effect, vi = se^2, measure = "MD", slab = study.name,
               method = "FE", data = dt))
forest(rma.uni(yi = cohensd, vi = var.cohensd, measure = "SMD", slab = study.name,
               method = "FE", data = dt))
dt %>% mutate(se.smd = sqrt(var.smd)) %>% ggplot(aes(x = se.smd, y = smd, col = total1 + total2)) + geom_point()
dt %>% ggplot(aes(x = se, y = effect, col = total1 + total2)) + geom_point()



funnel(metagen(TE = smd, seTE = sqrt(var.smd), data = mtd))
funnel(metagen(TE = effect, seTE = se, data = mtd))
# forest(rma.uni(yi = z, vi = var.z, measure = "ZCOR", slab = study.name,
#                method = "FE", data = dt))


mtd <- data.ext2 %>% filter(meta.id == meta.bin$meta.id[33])
funnel(metagen(TE = smd.pool, seTE = se.smd.pool, data = mtd))
funnel(metacor(cor = z, n = total1 + total2, studlab = study.name, data = dt, sm = "ZCOR"))
funnel(metagen(TE = lrr, seTE = sqrt(var.lrr), data = mtd))

plot(mtd$smd.pool, x = mtd$se.smd.pool, col = factor(mtd$study.id), pch = 20)
plot(mtd$z, x = sqrt(mtd$var.z), col = factor(mtd$study.id), pch = 20)
rank(mtd$z); rank(mtd$smd.pool); rank(mtd$effect)
rank(1/mtd$n); rank(sqrt(mtd$var.z)); rank(sqrt(mtd$var.smd.ordl))


#Variance/study sample size check: 
#Correspondance of variance of different measures with study sample size:
#----------------------------------------------------------------------------------------------#
data.ext2 <- data.ext2 %>% mutate(n = total1 + total2)

#Play through two scenarios; large difference, large error and vice versa:
mean1 <- c(30, 15); mean2 <- c(10, 10)
sd1 <- c(10, 10); sd2 <- c(10, 10)
total1 <- c(20, 100); total2 <- c(20, 100)

escalc(measure = "MD", m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, n1i = total1, n2i = total2)
escalc(measure = "SMD", m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, n1i = total1, n2i = total2)


#----------------------------------------------------------------------------------------------#
# examples <- data.ext2 %>% filter(outcome.measure.merged == "MD")
# examples <- examples[c(4987, 33200, 660),]
examples <- data.ext2 %>% filter(meta.id == meta.cont$meta.id[59])

plot(x = examples$se, y = examples$effect)
abline(lm(effect~ se, data = examples))
plot(x = examples$se.smd.pool, y = examples$smd.pool)
abline(lm(smd.pool ~ se.smd.pool, data = examples))
plot(x = examples$n, y = examples$effect)
abline(lm(effect ~ n, data = examples))
examples$se.z <- sqrt(examples$var.z)
plot(y = examples$z, x = examples$se.z)
abline(lm(z ~ se.z, data = examples))


plot(x = examples$n, y = examples$se)
plot(x = examples$n, y = examples$se.smd.pool)
plot(x = examples$effect, y = examples$smd.pool)

#----------------------------------------------------------------------------------------------#
#pvalue test results:
example.md <- data.ext2 %>% filter(meta.id == meta.cont$meta.id[59])

metabias(metagen(TE = effect, seTE = se, data = example.md))
metabias(metagen(TE = smd.pool, seTE = se.smd.pool, data = example.md))

