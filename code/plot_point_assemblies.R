reg.dat %>% filter(se.smd > 1.08 & se.smd < 1.15) %>% distinct(meta.id)


reg.dat %>% filter(se.smd > 1.08 & se.smd < 1.15) %>% distinct(id)


x <- reg.dat %>% filter(se.smd > 0.75 & se.smd < 0.85) %>% filter(smd > -0.1 & smd < 0.1)
x %>% distinct(meta.id)
x %>% distinct(id)



plot(x$smd, x = x$se.smd, col = x$meta.id, pch = 20)

x$n <- x$total1 + x$total2

hist(x$n)

x %>% ggplot(aes(x = se.smd, y = smd, color = factor(meta.id))) + geom_point() + theme(legend.position = "none")


x %>% ggplot(aes(x = se.smd, y = smd, color = n)) + geom_point() +
  scale_color_continuous(low = "darkgreen", high = "green")#+ theme(legend.position = "none")

x %>% ggplot(aes(x = se.smd, y = smd, color = factor(outcome.measure.merged))) + geom_point() 
  
x %>% ggplot(aes(x = se.smd, y = smd, color = factor(events1 + events2))) + geom_point() 

x %>% ggplot(aes(x = se.smd, y = smd, color = factor(events1 + events2))) + geom_point() 

  
ct <- x %>% filter(outcome.flag == "DICH")
plot(ct$events1, ct$n)
plot(ct$events1, ct$events2)
hist(ct$events1 + ct$events2, breaks = 55)


x %>% ggplot(aes(x = se.smd, y = smd, color = factor(events1))) + geom_point() 

x %>% ggplot(aes(x = se.smd, y = smd, color = factor(events2))) + geom_point() 
