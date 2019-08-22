grid.arrange(diff.z.fixef, diff.z.ranef,
             diff.d.fixef, diff.d.ranef,
             diff.log.hazard.ratio.fixef, diff.log.hazard.ratio.ranef, ncol = 2)


#Plots separated by Method:


sig.level <- 0.1
meta.f <- meta.f %>% rowwise() %>% 
  mutate(egger.test = ifelse(pval.egger < sig.level, 1, 0),
         thompson.test = ifelse(pval.thompson < sig.level, 1, 0),
         begg.test = ifelse(pval.begg < sig.level, 1, 0),
         
         schwarzer.test = ifelse(pval.schwarzer < sig.level, 1, 0),
         rucker.test = ifelse(pval.rucker < sig.level, 1, 0),
         harbord.test = ifelse(pval.harbord < sig.level, 1, 0),
         peter.test = ifelse(pval.peter < sig.level, 1, 0),
         I2 = max(c(0, 1 - (k-1)/Q)),
         i2f = factor(ifelse(I2 == 0, 1, 0)),
         se.stat = cut(abs(ifelse(outcome.type == "bin", stat.rucker, stat.thompson)),
                      breaks =  c(0,1,2, 7)),
         se.stat.ct = abs(ifelse(outcome.type == "bin", stat.rucker, stat.thompson)))

meta.bin <- meta.f %>% filter(outcome.type == "bin")
meta.cont <- meta.f %>% filter(outcome.type == "cont")
meta.surv <- meta.f %>% filter(outcome.type == "surv")

#Meta-Analysis and adjusted treatment effect estimate difference:

#Z-score:
diff.z.fixef <- meta.f %>% 
  mutate(fixef = est.z.fixef, ranef = est.z.ranef,
         copas = est.z.copas, regression = est.z.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(fixef) ~  abs(copas),
                           sign(copas) != sign(fixef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(fixef) ~  abs(regression),
                                sign(regression) != sign(fixef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  select(se.stat, fixef, ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.z.score", copas:regression) %>% 
  mutate(difference = fixef - adjusted.z.score) %>% 
  filter(difference > -.75 & difference < .75) %>% 
  ggplot(aes(x = difference, fill = se.stat)) + geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values  = c("seagreen1", "seagreen3", "green4")) +
  xlab(expression(paste("Fixed effects - adjusted ", z, "-score"))) 
#--------------------------------------------------------------------------------------------------------------------#

diff.z.ranef <- meta.f %>% 
  mutate(fixef = est.z.fixef, ranef = est.z.ranef,
         copas = est.z.copas, regression = est.z.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(ranef) ~  abs(copas),
                           sign(copas) != sign(ranef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(ranef) ~  abs(regression),
                                sign(regression) != sign(ranef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  select(se.stat, fixef, ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.z.score", copas:regression) %>% 
  mutate(difference = ranef - adjusted.z.score) %>% 
  filter(difference > -.75 & difference < .75) %>% 
  ggplot(aes(x = difference, fill = se.stat)) + geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values  = c("seagreen1", "seagreen3", "green4")) + ylab("") + 
  xlab(expression(paste("Random effects - adjusted ", z, "-score")))
#--------------------------------------------------------------------------------------------------------------------#

#SMD:
diff.d.fixef <- meta.f %>% 
  mutate(fixef = est.d.fixef, ranef = est.d.ranef,
         copas = est.d.copas, regression = est.d.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(fixef) ~  abs(copas),
                           sign(copas) != sign(fixef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(fixef) ~  abs(regression),
                                sign(regression) != sign(fixef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  select(se.stat, fixef, ranef, copas, regression) %>% 
  gather(key = "method", value = "adjusted.smd", copas:regression) %>% 
  mutate(difference = fixef - adjusted.smd) %>% 
  filter(difference > -1 & difference < 1) %>% 
  ggplot(aes(x = difference, fill = se.stat)) + geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values  = c("seagreen1", "seagreen3", "green4")) + 
  xlab(expression(paste("Fixed effects - adjusted Hedges ", g)))
#--------------------------------------------------------------------------------------------------------------------#

diff.d.ranef <- meta.f  %>% 
  mutate(fixef = est.d.fixef, ranef = est.d.ranef,
         copas = est.d.copas, regression = est.d.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(ranef) ~  abs(copas),
                           sign(copas) != sign(ranef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(ranef) ~  abs(regression),
                                sign(regression) != sign(ranef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  gather(key = "method", value = "adjusted.smd", copas:regression) %>% 
  mutate(difference = ranef - adjusted.smd) %>% 
  filter(difference > -1 & difference < 1) %>% 
  ggplot(aes(x = difference, fill = se.stat)) + geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + theme(legend.position = "none") +  
  scale_fill_manual(values  = c("seagreen1", "seagreen3", "green4")) +  ylab("")
  xlab(expression(paste("Random effects - adjusted Hedges ", g)))
#--------------------------------------------------------------------------------------------------------------------#

#Log hazard ratios:
diff.log.hazard.ratio.fixef <- meta.surv %>% 
  mutate(fixef = est.fixef, ranef = est.ranef,
         copas = est.copas, regression = est.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(fixef) ~  abs(copas),
                           sign(copas) != sign(fixef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(fixef) ~  abs(regression),
                                sign(regression) != sign(fixef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  gather(key = "method", value = "adjusted.log.hazard.ratio", copas:regression)  %>% 
  ggplot(aes(x = (fixef - adjusted.log.hazard.ratio), fill = se.stat)) + 
  geom_histogram(binwidth = 0.1, center = 0) +
  facet_wrap(~method) + theme_bw() + theme(legend.position = "none") +  
  scale_fill_manual(values  = c("seagreen1", "seagreen3", "green4")) + 
                      xlab("Fixed effects - adjusted log hazard ratio")
#--------------------------------------------------------------------------------------------------------------------#

diff.log.hazard.ratio.ranef <- meta.surv %>% 
  mutate(fixef = est.fixef, ranef = est.ranef,
         copas = est.copas, regression = est.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(ranef) ~  abs(copas),
                           sign(copas) != sign(ranef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(ranef) ~  abs(regression),
                                sign(regression) != sign(ranef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  gather(key = "method", value = "adjusted.log.hazard.ratio", copas:regression)  %>% 
  ggplot(aes(x = (ranef - adjusted.log.hazard.ratio), fill = se.stat)) + 
  geom_histogram(binwidth = 0.1, center = 0) + ylab("") +
  facet_wrap(~method) + theme_bw() + guides(fill=guide_legend(title="small study effect test stat.")) +
    scale_fill_manual(values  = c("seagreen1", "seagreen3", "green4")) + 
                        xlab("Random effects - adjusted log hazard ratio")




grid.arrange(diff.z.fixef, diff.z.ranef,
             diff.d.fixef, diff.d.ranef,
             diff.log.hazard.ratio.fixef, diff.log.hazard.ratio.ranef, ncol = 2)
