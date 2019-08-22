mid <- median(meta.f$pval.se)

meta.f <- meta.f %>% mutate(pval.se1 = case_when(outcome.flag == "DICH" ~ pval.harbord,
                                                 outcome.flag == "CONT" ~ pval.egger,
                                                 TRUE ~ NA_real_))


method.names <- c(z.reg = "regression", z.copas = "selection model")

p <- meta.f %>% filter(outcome.flag != "IV") %>% 
  select(pval.se, z.fixef, z.ranef, z.reg, z.copas) %>% 
  gather(key = "method", value = "test.stat", z.reg:z.copas) %>%
  ggplot(aes(x = test.stat, y = z.fixef, colour = pval.se)) + 
  facet_wrap(~ method, labeller = as_labeller(method.names)) + 
  scale_color_gradient2(midpoint=mid, low="firebrick4", mid = "gold",
                        high="chartreuse2", space ="Lab") +
  geom_point(size = .7) + theme_bw() + theme(legend.position = "none") + 
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = .3) +
  xlab("adjusted test statistic") + ylab("fixed effects test statistic")

meta.f %>% filter(outcome.flag != "IV") %>% 
  select(pval.se, z.fixef, z.ranef, z.reg, z.copas) %>% 
  gather(key = "method", value = "test.stat", z.reg:z.copas) %>%
  ggplot(aes(x = test.stat, y = z.ranef#, colour = pval.se1
             )) + 
  facet_wrap(~ method, labeller = as_labeller(method.names)) + 
  #scale_color_gradient2(midpoint=mid, low="firebrick4", mid = "gold",
  #                      high="chartreuse2", space ="Lab", name = expression(paste(italic(p),"-value"))) +
  geom_point(size = 1, color = "gray29") + theme_bw() + 
  ggtitle("Evidence for treatment efficacy") +
  theme(legend.position = "none",
    axis.title=element_text(size=20, margin = margin(r = 20, b = 20)), 
    axis.title.y = element_text(size=40, margin = margin(r = 20)), 
    axis.title.x = element_text(size=40, margin = margin(t = 20)), 
    legend.text = element_text(size = 20), 
        legend.title = element_text(size = 40, face = "bold"), title = element_text(size = 22),
    strip.text.x = element_text(size = 40),
    axis.text= element_text(size=30, color = "black"),
        #rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
    strip.background = element_rect(fill = "transparent", linetype = 0),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 40, margin = margin(b = 20)),
    axis.line = element_line(size = 1, linetype = "solid", color = "black"),
    panel.border = element_blank()) +
  #guides(fill=guide_legend(title="p-value")) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = .3) +
  xlab("adjusted test statistic") + ylab("meta analysis test statistic") 
    

ggsave(p, filename = "tr_tst2.png",  bg = "transparent", width = 16.5, height = 8.5)


p.s <- meta.f %>% filter(outcome.flag != "IV") %>% 
  mutate(fixef = est.d.fixef, ranef = est.d.ranef,
         copas = est.d.copas, regression = est.d.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(ranef) ~  abs(copas),
                           sign(copas) != sign(ranef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(ranef) ~  abs(regression),
                                sign(regression) != sign(ranef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  gather(key = "method", value = "adjusted.smd", copas) %>% 
  mutate(difference = ranef - adjusted.smd) %>% 
  filter(difference > -1 & difference < 1) %>% 
  ggplot(aes(x = difference#, fill = se.stat
             )) + geom_histogram(binwidth = 0.1, center = 0) +
  theme_bw() + #guides(fill=guide_legend(title="test stat.")) + 
  ggtitle("Selection model") +
  xlim(c(-0.7, 1)) +
  theme(legend.position = "none",
        axis.title=element_text(size=20), 
        axis.title.y = element_text(size=40, margin = margin(r = 20)), 
        axis.title.x = element_text(size=40, margin = margin(t = 20)), 
        legend.text = element_text(size = 20), 
        axis.text= element_text(size=30),
        strip.text.x = element_text(size = 30),
        #rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40, margin = margin(b = 20)),
        axis.line = element_line(size = 1, linetype = "solid"),
        panel.border = element_blank()) +
  scale_fill_manual(values  = c("seagreen1", "seagreen3", "green4")) + ylab("count") +
  xlab(expression(paste("meta analyis ",  italic(d)," - adjusted ", italic(d)))) + 
  geom_vline(xintercept = 0, linetype = "dashed" ) +
  geom_text(aes(x = 0.4, y = 600, label = "75% >= 0"), size = 12)


p.r <- meta.f %>% filter(outcome.flag != "IV") %>% 
  mutate(fixef = est.d.fixef, ranef = est.d.ranef,
         copas = est.d.copas, regression = est.d.reg) %>% 
  mutate(copas = case_when(sign(copas) == sign(ranef) ~  abs(copas),
                           sign(copas) != sign(ranef) ~ -abs(copas)),
         regression = case_when(sign(regression) == sign(ranef) ~  abs(regression),
                                sign(regression) != sign(ranef) ~ -abs(regression)),
         fixef = abs(fixef),
         ranef = abs(ranef)) %>% 
  gather(key = "method", value = "adjusted.smd", regression) %>% 
  mutate(difference = ranef - adjusted.smd) %>% 
  filter(difference > -1 & difference < 1) %>% 
  ggplot(aes(x = difference#, fill = se.stat
             )) + geom_histogram(binwidth = 0.1, center = 0) +
  theme_bw() + #guides(fill=guide_legend(title="test stat.")) 
  ggtitle("Regression") +
  xlim(c(-0.7, 1)) +
  theme(legend.position = "none",
        axis.title=element_text(size=20), 
        axis.text= element_text(size=30),
        axis.title.y = element_text(size= 40, margin = margin(r = 20)), 
        axis.title.x = element_text(size= 40, margin = margin(t = 20)), 
        strip.text.x = element_text(size = 30),
        #rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 40, margin = margin(b = 20)),
        axis.line = element_line(size = 1, linetype = "solid"),
        panel.border = element_blank(),
        axis.ticks = element_line(size = 1)) +
  scale_fill_manual(values  = c("seagreen1", "seagreen3", "green4")) +  ylab("") +
  xlab(expression(paste("meta analyis ",  italic(d)," - adjusted ", italic(d)))) + 
   geom_vline(xintercept = 0, linetype = "dashed" ) +
   geom_text(aes(x = 0.5, y = 200, label = "67% > 0"), size = 12)

p <- grid.arrange(p.s, p.r, ncol = 2)

ggsave(p, filename = "tr_tst3.png",  bg = "transparent", width = 16.5, height = 8.5)
