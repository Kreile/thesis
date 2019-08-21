idt <- meta.cont$meta.id[1]

idt <- meta.bin$meta.id[1]


mtd <- data.ext2 %>% filter(meta.id == idt)

# mtd <- escalc(measure = "PBIT", ai = events1, n2i = total2, ci = events2, n1i = total1, 
#               to = "only0",
#               add = 1/2, append = T, var.names = c("smd.pbit", "var.smd.pbit"), data = mtd)

#Issue starts here:
forest(rma.uni(yi = cor.pearson, vi = var.cor.pearson, measure = "COR",
               method = "FE", data = mtd))
forest(rma.uni(yi = z, vi = var.z, measure = "ZCOR",
               method = "FE", data = mtd), layout = "JAMA")

# mtd$smd.ordl[1] <- 1.1547
# mtd$var.smd.ordl[1] <- 0.055
# mtd$total1[1] <- 100; mtd$total2[1] <- 100
# 
# mtd <- mtd %>% mutate(
#   a = ((total1 + total2)^2)/(total1*total2),
#   cor.pearson  = smd.ordl/sqrt((smd.ordl^2) + a),
#   var.cor.pearson = ((a^2)*var.smd.ordl)/(((smd.ordl^2)+a)^3))
# 
# mtd <- mtd %>% 
#   mutate(z = 0.5 * log( (1 + cor.pearson)/(1 - cor.pearson) , base = exp(1)), 
#          var.z= 1/(total1 + total2 - 3))


par(mfrow = c(1,2))
forest(rma.uni(yi = lrr, vi = var.lrr,
               measure = "RR", method = "FE", data = mtd))
forest(rma.uni(yi = smd.ordl, vi = var.smd.ordl, measure = "SMD",
               method = "FE", data = mtd))
forest(rma.uni(yi = smd.ordl, vi = var.smd.ordl, measure = "SMD",
               method = "FE", data = mtd))


par(mfrow = c(1, 1))
funnel(metagen(TE = smd, seTE = sqrt(var.smd), data = mtd))
funnel(metagen(TE = effect, seTE = se, data = mtd))
funnel.id(idt)
funnel(metagen(TE = z, seTE = sqrt(var.z), data = mtd))

# forest.meta(metabin(event.e = events1c, n.e = total1, event.c = events2c, n.c = total2,
#                     studlab = study.name, sm = "RR",
#                     method = "Inverse", data = mtd))
# forest.meta(metacor(cor = z, n = total1 + total2, studlab = study.name, data = mtd, sm = "ZCOR"))


#----------------#
mtd <- data.ext2 %>% filter(meta.id == meta.cont$meta.id[11])

forest(rma.uni(yi = effect, sei = se,, method = "FE", data = mtd))
forest(rma.uni(yi = smd, vi = var.smd,, method = "FE", data = mtd))
forest(rma.uni(yi = cor.pearson, vi = var.cor.pearson,
               method = "FE", data = mtd))
forest(rma.uni(yi = z, vi = var.z, 
               method = "FE", data = mtd))
