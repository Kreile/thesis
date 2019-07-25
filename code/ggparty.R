#Load data:
rm(list = ls())
PATH_HOME = path.expand("~") # user home
PATH = file.path(PATH_HOME, 'Data/PubBias')
PATH2 = file.path(PATH_HOME, 'PubBias')
FILE = 'cochrane_2019-07-04.csv'
PATH_DATA = file.path(PATH, 'data')
PATH_CODE = file.path(PATH2, 'code')
PATH_RESULTS = file.path(PATH2, 'results_new')
PATH_FIGURES = file.path(PATH_RESULTS, 'figures')

require(nlme)

source(file.path(PATH_CODE, 'PubBias_functions.R'))

load(file.path(PATH_RESULTS, "data2.RData"))
# data.ext2 <- pb.process3(data)
load(file.path(PATH_RESULTS, "data_used_for_analysis.RData"))

load(file.path(PATH_RESULTS, "meta_analyses_summary_complete.RData"))

load(file.path(PATH_RESULTS, "meta_id_vector.RData"))

load(file.path(PATH_DATA, "PubBias_2019-07-19.RData"))

meta.f$heterogeneity <- ifelse(meta.f$I2 == 0, 1, 0)
meta.f <- meta.f %>% mutate(p.value = case_when(outcome.flag == "DICH" & heterogeneity == 1 ~ pval1.rucker,
                                                outcome.flag == "DICH" & heterogeneity == 0 ~ pval1.rucker.linreg,
                                                outcome.flag == "CONT" & heterogeneity == 1 ~ pval1.thompson,
                                                outcome.flag == "CONT" & heterogeneity == 0 ~ pval1.egger,
                                                outcome.flag == "IV" & heterogeneity == 1 ~ pval1.thompson,
                                                outcome.flag == "IV" & heterogeneity == 0 ~ pval1.egger),
                            
                            classes = case_when(outcome.flag == "DICH" & heterogeneity == 1 ~ 1,
                                                outcome.flag == "DICH" & heterogeneity == 0 ~ 2,
                                                outcome.flag == "CONT" & heterogeneity == 1 ~ 3,
                                                outcome.flag == "CONT" & heterogeneity == 0 ~ 4,
                                                outcome.flag == "IV" & heterogeneity == 1 ~ 5,
                                                outcome.flag == "IV" & heterogeneity == 0 ~ 6))
levels(meta.f$outcome.flag) <- list("CONT" = "CONT", "DICH" = "DICH", "IV" = "IV")

meta.f %>% ggplot(aes(x = p.value)) + facet_wrap(~ classes, nrow = 2) + geom_histogram(boundary = 0, binwidth = 0.1)
#----------------------------------------------------------------------------------------------#
require(ggparty)

tr_tree <- lmtree(p.value ~ outcome.flag + heterogeneity, data = meta.f)
ggparty(tr_tree) +
  geom_edge() +
  geom_edge_label() +
  geom_node_splitvar() +
  geom_node_plot(gglist = list(geom_histogram(aes(x = p.value)),
                               theme_bw(base_size = 10)),
                 shared_axis_labels = TRUE,
                 legend_separator = TRUE,
                 # predict based on variable
                 #predict = "beauty",
                 # graphical parameters for geom_line of predictions
                 #predict_gpar = list(col = "blue",
                 #size = 1.2)
  )


stump <- partynode(id = 1L, 
                   split = partysplit(which(names(meta.f) == "outcome.flag"), index = 1:3), 
                   kids = lapply(2:4, partynode))
meta.f2 <- meta.f %>% filter(!is.na(p.value))
plot(party(stump, meta.f))

table(kidids_node(stump, meta.f), fitted_node(stump, meta.f))
party(stump, meta.f)


ggparty(party = party(stump, meta.f)) +
  geom_edge() +
  geom_edge_label() +
  geom_node_splitvar() +
  geom_node_plot(gglist = list(geom_histogram(aes(x = p.value)),
                               theme_bw(base_size = 10)),
                 shared_axis_labels = TRUE,
                 legend_separator = TRUE,
                 # predict based on variable
                 #predict = "beauty",
                 # graphical parameters for geom_line of predictions
                 #predict_gpar = list(col = "blue",
                 #size = 1.2)
  )

#----------------------------------------------------------------------------------------------#
sp_o <- partysplit(which(names(meta.f) == "outcome.flag"), index = 1:3)
sp_t <- partysplit(which(names(meta.f) == "I2"), breaks = 0, index = 1:2)

n1 <- partynode(id = 1L, split = sp_o, kids = list(
  partynode(2L, split = sp_t, kids = lapply(5:6, partynode), info = "DICH"), 
  partynode(3L, split = sp_t, kids = lapply(7:8, partynode), info = "CONT"), 
  partynode(4L, split = sp_t, kids = lapply(9:10, partynode), info = "IV")))
t2 <- party(n1, data = meta.f)

t2 <- party(n1, data = meta.f,
            fitted = data.frame(
              "(fitted)" = fitted_node(n1, data = meta.f),
              "(response)" = meta.f$p.value,
              check.names = FALSE),
            terms = terms(p.value ~ outcome.flag + I2, data = WeatherPlay))
t2 <- as.constparty(t2)
plot(t2)
# t2 <- as.constparty(t2)

kidids_split(sp_o, meta.f)
table(kidids_node(n1, meta.f), fitted_node(n1, WeatherPlay))

ggparty(t2) +
  geom_edge() +
  geom_edge_label() +
  geom_node_label(aes(label = splitvar), ids = "inner") +
  # identical to  geom_node_splitvar() +
  geom_node_label(aes(label = info), ids = "terminal")






# sp_o <- partysplit(1L, index = 1:3)
# n1 <- partynode(id = 1L, split = sp_o, kids = lapply(2L:4L, partynode))
# t2 <- party(n1,
#             data = WeatherPlay,
#             fitted = data.frame(
#               "(fitted)" = fitted_node(n1, data = WeatherPlay),
#               "(response)" = WeatherPlay$play,
#               check.names = FALSE),
#             terms = terms(play ~ ., data = WeatherPlay)
# )
# t2 <- as.constparty(t2)
# ggparty(t2) +
#   geom_edge() +
#   geom_edge_label() +
#   geom_node_splitvar() +
#   # pass list to gglist containing all ggplot components we want to plot for each
#   # (default: terminal) node
#   geom_node_plot(gglist = list(geom_bar(aes(x = "", fill = play),
#                                         position = position_fill()),
#                                xlab("play")))
# table(kidids_node(n1, WeatherPlay), fitted_node(n1, WeatherPlay))

sp_o <- partysplit(1L, index = 1:3)
n1 <- partynode(id = 1L, split = sp_o, kids = lapply(2L:4L, partynode))
t2 <- party(n1,
            data = meta.f,
            fitted = data.frame(
              "(fitted)" = fitted_node(n1, data = meta.f),
              "(response)" = meta.f$p.value,
              check.names = FALSE),
            terms = terms(p.value ~ outcome.flag, data = meta.f)
)
t2 <- as.constparty(t2)
ggparty(t2) +
  geom_edge() +
  geom_edge_label() +
  geom_node_splitvar() +
  # pass list to gglist containing all ggplot components we want to plot for each
  # (default: terminal) node
  geom_node_plot(gglist = list(geom_histogram(aes(x = p.value),
                                        position = position_fill()),
                               xlab("play")))












data("WeatherPlay", package = "partykit")
sp_o <- partysplit(1L, index = 1:3)
sp_om <- partysplit(3L, index = 1:3)

sp_h <- partysplit(3L, breaks = 75)
sp_w <- partysplit(4L, index = 1:2)

sp_oh <- partysplit(2L, breaks = 0, index = 1:2)
sp_wh <- partysplit(3L, breaks = 0)
sp_sh <- partysplit(4L, breaks = 0)

pnm <- partynode(1L, split = sp_om, kids = list(
  partynode(2L, split = sp_oh, kids = list(
              partynode(5L),
              partynode(6L))),
  partynode(3L, split = sp_oh, kids = list(
    partynode(7L),
    partynode(8L))),
  partynode(4L, split = sp_oh, kids = list(
    partynode(9L),
    partynode(10L)))
))

pnm2 <- partynode(1L, split = sp_om, kids = lapply(2L:4L, partynode))

pym2 <- party(pnm2, meta.f) 
pym <- party(pnm, meta.f) 
# t1m <- as.constparty(pym2)
# ggparty(pym2) +
#   geom_edge() +
#   # geom_edge_label() +
#   # geom_node_splitvar() +
#   geom_node_plot(gglist = list(geom_histogram(aes(x = p.value)),
#                                theme_bw(base_size = 10)),
#                  shared_axis_labels = TRUE,
#                  legend_separator = TRUE
#                  # predict based on variable
#                  #predict = "beauty",
#                  # graphical parameters for geom_line of predictions
#                  #predict_gpar = list(col = "blue",
#                  #size = 1.2)
#   )

ggparty(pym2) +
  geom_edge() +
  geom_edge_label() +
  geom_node_label(aes(label = splitvar), ids = "inner") +
  # identical to  geom_node_splitvar() +
  geom_node_label(aes(label = info), ids = "terminal")



ggparty(pym) +
  geom_edge() +
  geom_edge_label() +
  geom_node_label(aes(label = splitvar), ids = "inner") +
  # identical to  geom_node_splitvar() +
  geom_node_label(aes(label = info), ids = "terminal")


ggparty(party(pnm, meta.f))

pn <- partynode(1L, split = sp_o, kids = list(
  partynode(2L, split = sp_h, kids = list(
    partynode(3L, info = "yes"),
    partynode(4L, info = "no"))),
  partynode(5L, info = "yes"),
  partynode(6L, split = sp_w, kids = list(
    partynode(7L, info = "yes"),
    partynode(8L, info = "no")))))
py <- party(pn, WeatherPlay)




n1 <- partynode(id = 1L, split = sp_o, kids = lapply(2L:4L, partynode))
t2 <- party(n1,
            data = WeatherPlay,
            fitted = data.frame(
              "(fitted)" = fitted_node(n1, data = WeatherPlay),
              "(response)" = WeatherPlay$play,
              check.names = FALSE),
            terms = terms(play ~ ., data = WeatherPlay)
)
t2 <- as.constparty(t2)

ggparty(t2) +
  geom_edge() +
  geom_edge_label() +
  geom_node_splitvar() +
  # pass list to gglist containing all ggplot components we want to plot for each
  # (default: terminal) node
  geom_node_plot(gglist = list(geom_bar(aes(x = "", fill = play),
                                        position = position_fill()),
                               xlab("play")))




data("WeatherPlay", package = "partykit")
require(ggparty)
sp_o <- partysplit(1L, index = 1:3)
sp_h <- partysplit(3L, breaks = 75)
sp_w <- partysplit(4L, index = 1:2)
pn <- partynode(1L, split = sp_o, kids = list(
  partynode(2L, split = sp_h, kids = list(
    partynode(3L, info = "yes"),
    partynode(4L, info = "no"))),
  partynode(5L, info = "yes"),
  partynode(6L, split = sp_w, kids = list(
    partynode(7L, info = "yes"),
    partynode(8L, info = "no")))))
py <- party(pn, WeatherPlay)
plot(py)




data("TeachingRatings", package = "AER")
tr <- subset(TeachingRatings, credits == "more")

tr_tree <- lmtree(eval ~ beauty | minority + age + gender + division + native +
                    tenure, data = tr, weights = students, caseweights = FALSE)
ggparty(tr_tree) +
  geom_edge() +
  geom_edge_label() +
  geom_node_splitvar() +
  geom_node_plot(gglist = list(geom_point(aes(x = beauty,
                                              y = eval,
                                              col = tenure,
                                              shape = minority),
                                          alpha = 0.8),
                               theme_bw(base_size = 10)),
                 shared_axis_labels = TRUE,
                 legend_separator = TRUE,
                 # predict based on variable
                 predict = "beauty",
                 # graphical parameters for geom_line of predictions
                 predict_gpar = list(col = "blue",
                                     size = 1.2)
  )




n1 <- partynode(id = 1L, split = sp_o, kids = lapply(2L:4L, partynode))
t2 <- party(n1,
            data = WeatherPlay,
            fitted = data.frame(
              "(fitted)" = fitted_node(n1, data = WeatherPlay),
              "(response)" = WeatherPlay$play,
              check.names = FALSE),
            terms = terms(play ~ ., data = WeatherPlay)
)
t2 <- as.constparty(t2)

m1 <- partynode(id = 1L, split = partysplit(varid = 1L, index = c(1L,3L,4L)), kids= lapply(1L:3L, partynode))

sp_m <- partysplit(1L, varid = 3)

sp_o <- partysplit(1L, index = 1:3)
sp_h <- partysplit(3L, breaks = 75)
sp_w <- partysplit(4L, index = 1:2)
pn <- partynode(1L, split = sp_o, kids = list(
  partynode(2L, split = sp_h, kids = list(
    partynode(3L, info = "yes"),
    partynode(4L, info = "no"))),
  partynode(5L, info = "yes"),
  partynode(6L, split = sp_w, kids = list(
    partynode(7L, info = "yes"),
    partynode(8L, info = "no")))))
py <- party(pn, WeatherPlay)