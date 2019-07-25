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




