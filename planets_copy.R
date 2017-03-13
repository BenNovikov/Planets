Sys.setlocale('LC_ALL', 'ukrainian')
#  install.packages("radiant", repos = "https://radiant-rstats.github.io/minicran/") #rq
#  install.packages("ggplot2")
#  install.packages("ggdendro")
#  install.packages("dplyr")
#  install.packages("tidyr")
library(radiant)
library(ggplot2)
library(ggdendro)
library(dplyr)
library(tidyr)
library(cluster)

# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)
# setting parameter names and cleaning leading spaces
columnSet <- function (columns) {
  columns <- rbind(as.data.frame(t(c("rowid", "Row ID"))), columns)
  colnames(columns) <- c("ColumnID", "Parameter")
  columns$Parameter <- trim.leading(columns$Parameter)
  return(columns)
}
# Multiple plot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## ----
qqplot.data <- function (vec) # argument: vector of numbers
{
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  d <- data.frame(resids = vec)
  ggplot(d, aes(sample = resids)) + 
    stat_qq() + 
    geom_abline(slope = slope, intercept = int, col =2) +
    labs(x = "Теоретичний нормальний розподіл", y = "Розподіл вибірки")
}

# reading data
#setwd("~/Desktop/SofieRFiles")
setwd("F:/SofieRFiles")
planets <- read.csv("raw_confirmed.csv")
# kepler <- read.csv("raw_kepler.csv")
# k2 <- read.csv("raw_k2.csv")
colPlanets <- columnSet(read.csv("raw_conf_col.csv", FALSE, sep=":"))
# colKepler <- columnSet(read.csv("raw_kepler_col.csv", FALSE, sep=":"))
# colK2 <- columnSet(read.csv("raw_k2_col.csv", FALSE, sep=":"))
#write.csv(columns, file = "colPlanets.csv")

# get list of unique discovery methods 
discoveryMethods.x <- sort(as.character(unique(planets$pl_discmethod)), decreasing = FALSE)
discoveryMethods.x
discoveryMethods.ua <- c("астрометрія", 
                 "ETV",
                 "прямі спостереження",
                 "мікролінзування",
                 "модуляція орбітального блиску",
                 "спостереження за пульсарами",
                 "PTV",
                 "доплерівська спектроскопія", 
                 "транзитний метод", 
                 "TTV")

# split cumulative dataframe by discovery methods
discoveryMethods.filename <- gsub(" ", "", paste0("planets", discoveryMethods.x, ".csv"))
discoveryMethods.filename
splittedPlanetsDF <- list()
for (i in 1:length(discoveryMethods.x)) { 
  splittedPlanetsDF[[i]] <- planets[planets$pl_discmethod == discoveryMethods.x[i], ]
  write.csv(splittedPlanetsDF[[i]], file = discoveryMethods.filename[i])
}

# number of discovered planets til 2016 and 2013
planetsData <- arrange(count(na.omit(planets[, c('pl_discmethod', 'pl_disc')])), pl_disc)
pld2016 <- ggplot(planetsData, aes(x = pl_disc, y = freq, fill = pl_discmethod)) + 
  geom_bar(stat = "identity") +
  ggtitle("Кількість відкритих планет. 1989-2016.") +
  labs(x = "рік відкриття", y = "кількість відкритих екзопланет") + 
  scale_x_continuous(breaks = seq(min(planetsData$pl_disc),max(planetsData$pl_disc),9)) + 
  scale_y_continuous(breaks = seq(0,1500,500), limits = c(0, 1550))

planetsData <- planetsData[planetsData$pl_disc < 2014, ]
pld2013 <- ggplot(planetsData, aes(x = pl_disc, y = freq, fill = pl_discmethod)) + 
  geom_bar(stat = "identity") +
  ggtitle("Кількість відкритих планет. 1989-2013.") +
  labs(x = "рік відкриття", y = "кількість відкритих екзопланет") + 
  scale_x_continuous(breaks = seq(min(planetsData$pl_disc),max(planetsData$pl_disc),8)) + 
  scale_y_continuous(breaks = seq(0,150,50), limits = c(0, 155))

multiplot(pld2013 + 
            scale_fill_discrete(name="методи:",
                                breaks=discoveryMethods.x,
                                labels=discoveryMethods.ua) +
            theme(legend.position = c(0.2,0.6), 
                  legend.justification=c(0.1,0.6), 
                  legend.background = element_rect(fill="transparent")),
          pld2016 + 
            theme(legend.position="none"),
          cols=2)

## Number of planets discovered with detection methods
## update planetsData
planetsData <- arrange(count(na.omit(planets[, c('pl_discmethod', 'pl_disc')])), pl_disc)
tempFunc <- function(x) sum(planetsData[planetsData$pl_discmethod == x,]$freq)
discoveryMethods.length <- length(discoveryMethods.x)
discoveryMethods.count <- c(NULL)
for(i in 1:discoveryMethods.length) {
  discoveryMethods.count[i] <- tempFunc(discoveryMethods.x[i])
}
discoveryMethods.data <- as.data.frame(discoveryMethods.ua)
colnames(discoveryMethods.data) <- 'method'
discoveryMethods.data[, 'count'] <- (discoveryMethods.count)

pld01 <- ggplot(data = discoveryMethods.data,aes(x = 'method', y = 'count')) + 
  geom_bar(stat = "identity") 
#+  labs(x = "рік відкриття", y = "кількість відкритих екзопланет") 
pld01

tempSum <- sum(planetsData$freq)
lbl8 <- paste0(discoveryMethods.ua[8],": ", round(tempFunc(discoveryMethods.x[8])/tempSum*100, 1), "%")
lbl9 <- paste0(discoveryMethods.ua[9],": ", round(tempFunc(discoveryMethods.x[9])/tempSum*100, 1), "%")
pld02 <- ggplot(planets, aes(x = "", labels = NULL, fill = planets$pl_discmethod)) +
  geom_bar(width = 1) + 
  coord_polar("y") + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none") + 
  annotate("text", label = lbl8, x = 1, y = 0, size = 6, colour = "black") +
  annotate("text", label = lbl9, x = .5, y = 0, size = 6, colour = "black") +
  scale_fill_grey()

multiplot(pld01 + 
            scale_fill_discrete(name="методи:",
                                breaks=discoveryMethods.x,
                                labels=discoveryMethods.ua) +
            theme(legend.position = c(0.2,0.6), 
                  legend.justification=c(0.1,0.6), 
                  legend.background = element_rect(fill="transparent")),
          pld02 + 
            theme(legend.position="none"),
          cols=2)

#
# Mass & Radius
#
#
## Add density column to planets data
planets[, "pl_density"] <- ((planets$pl_bmassj * 1.898E+30)/((planets$pl_radj * 7.149E+09)^3*4/3*pi))

## check the correlation 
result <- correlation(dataset = planets, c("pl_radj", "pl_bmassj", "pl_density"), method = "pearson", data_filter = "")
summary(result)
plot(result)

## Have a look at the densities
p1 <- ggplot(planets[!is.na(planets$pl_bmassj),], aes(x = pl_bmassj)) +
  geom_histogram(colour="black", binwidth=.5) + 
  ggtitle("Відомі маси екзопланет") +
  labs(x = "Маса планети (в масах Юпітера)", y = "Кількість спостережень")

p2 <- ggplot(planets[!is.na(planets$pl_radj),], aes(x = pl_radj)) +
  geom_histogram(colour="black", binwidth=.5) + 
  ggtitle("Відомі радіуси екзопланет") +
  labs(x = "Радіус планети (в радіусах Юпітера)", y = "Кількість спостережень")

p3 <- ggplot(planets[!is.na(planets$pl_bmasse),], aes(x = pl_bmasse)) +
  geom_histogram(colour="black", binwidth=.5) + 
  ggtitle("Відомі маси екзопланет") +
  labs(x = "Маса планети (в масах Землі, логарифмічна шкала)", y = "Кількість спостережень")

p4 <- ggplot(planets[!is.na(planets$pl_rade),], aes(x = pl_rade)) +
  geom_histogram(colour="black", binwidth=.5) + 
  ggtitle("Відомі радіуси екзопланет") +
  labs(x = "Радіус планети (в радіусах Землі, логарифмічна шкала)", y = "Кількість спостережень")

multiplot(p1, p2, p3 + scale_x_log10(), p4 + scale_x_log10(), cols=2)

p5 <- ggplot(planets[!is.na(planets$pl_density),], aes(x = pl_density)) +
  geom_histogram(colour="black", binwidth=.5) + 
  ggtitle("Обчислені значення густини екзопланет") +
  labs(x = "Густина планети (гр/см3, логарифмічна шкала)", y = "Кількість спостережень") 
# p5 + scale_x_log10()
p6 <- ggplot(planets[!is.na(planets$pl_density),], aes(x = pl_density)) +
  geom_density(colour="black", alpha=.3) +   
  ggtitle("Обчислені значення густини екзопланет") +
  labs(x = "Густина планети (гр/см3)", y = "Густина розподілу") 
# p6
multiplot(p6 + labs(x = "Густина планети (гр/см3)", y = "Густина розподілу"), 
          p5 + labs(x = "Густина планети (гр/см3, логарифмічна шкала)", y = "Кількість спостережень") + 
            scale_x_log10(), cols = 2)

shapiro.test(planets$pl_bmassj)
shapiro.test(planets$pl_radj)
shapiro.test(planets$pl_density)
## Plot using a qqplot
# qqnorm(planets$pl_bmassj);qqline(planets$pl_bmassj, col = 2)
# qqnorm(planets$pl_radj);qqline(planets$pl_radj, col = 2)
# qqnorm(planets$pl_density);qqline(planets$pl_density, col = 2)

q1 <- qqplot.data(planets[!is.na(planets$pl_bmassj), ]$pl_bmassj)
q2 <- qqplot.data(planets[!is.na(planets$pl_radj), ]$pl_radj)
q3 <- qqplot.data(planets[!is.na(planets$pl_density), ]$pl_density)
multiplot(q1 + ggtitle("Відомі маси екзопланет"),  
          q2 + ggtitle("Відомі радіуси екзопланет"),  
          q3 + ggtitle("Обчислена густина екзопланет"), 
          cols = 3)

# ## Homoscedacity test
# var.test(planets$pl_bmassj, planets$pl_density)

# ## Generate Normal
# words1 = rnorm(100)
# ## Have a look at the densities
# plot(density(words1))
# ## Perform the test
# shapiro.test(words1)
# ## Plot using a qqplot
# qqnorm(words1);qqline(words1, col = 2)

planetsData <- na.omit(planets[, c('pl_radj', 'pl_bmassj')])
## Hierarchical Clustering
## methods: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), 
## "median" (= WPGMC) or "centroid" (= UPGMC).
clHierarch <- hclust(dist(planetsData), method = "ward.D")
gd1 <- ggdendrogram(clHierarch, labels = FALSE, leaf_labels = FALSE, rotate = TRUE) + 
  labs(title="Непараметрична кластерізація за методом Варда", 
       subtitle = "10% лінія відсічі встановлює оптимальний рівень - 5 кластерів")
gd1  + geom_hline(aes(yintercept = max(ggplot_build(gd1)$data[[2]]$y) * 0.1))  

cluster <- cutree(clHierarch, k = 5)
## number of planets in each cluster
table(unlist((factor(cluster))))
## plot clustered data
x <- planetsData[, 1]
y <- planetsData[, 2]
Clusters <- unlist(factor(cluster))
pd1 <- ggplot(planetsData, aes(x, y, colour = Clusters)) + 
  geom_point(size = 3) +
  ggtitle("Непараметрична кластерізація на 5 груп за методом Варда. Графік розсіювання") +
  labs(x = "Радіус планети (в радіусах Юпітера)", y = "Маса планети/M*sin(i) (в масах Юпітера)") 
pd1 
#+ scale_x_log10() + scale_y_log10()
#+ facet_wrap(nrow =2, ~ cluster)
#+ facet_grid(. ~ cluster)

planetsData[, "cluster"] <- cluster
planetsData[, "mass_gr"] <- (planetsData$pl_bmassj * 1.898E+30)
planetsData[, "radius_cm"] <- planetsData$pl_radj * 7.149E+09
planetsData[, "volume_cm3"] <- (planetsData$radius_cm)^3*4/3*pi
planetsData[, "density"] <- planetsData$mass_gr/planetsData$volume_cm3
write.csv(planetsData, file = "massRadiusPlanetsData.csv")
#clusterData <- as.data.frame(table(planetsData$cluster))
# distinct(planetsData, cluster)

print(clusterSummary <- summarise(group_by(planetsData, cluster), 
          count=n(), 
          min_mass = min(pl_bmassj), 
          max_mass = max(pl_bmassj),
          median_mass = median(pl_bmassj),
          min_radius = min(pl_radj), 
          max_radius = max(pl_radj),
          median_radius = median(pl_radj),
          min_density = min(density), 
          max_density = max(density),
          median_density = median(density)))
write.csv(clusterSummary, file = "massRadiusClusterSummary.csv")
#### =====



planetsData <- na.omit(planets[, c('pl_radj', 'pl_bmassj')])
planetsData <- planetsData %>% filter(pl_bmassj < 25, pl_radj <= 2)
summarise(planetsData, count = n())
planetsData <- planetsData %>% select(pl_radj, pl_bmassj)
clHierarch <- hclust(dist(planetsData), method = "ward.D")
gd2 <- ggdendrogram(clHierarch, labels = FALSE, leaf_labels = FALSE, rotate = TRUE) + 
  labs(title="Непараметрична кластерізація за методом Варда", 
       subtitle = "10% лінія відсічі встановлює оптимальний рівень - 6 кластерів")
gd2   + geom_hline(aes(yintercept = max(ggplot_build(gd2)$data[[2]]$y) * 0.1))  
multiplot(gd1 + geom_hline(aes(yintercept = max(ggplot_build(gd1)$data[[2]]$y) * 0.1)), 
          gd2 + geom_hline(aes(yintercept = max(ggplot_build(gd2)$data[[2]]$y) * 0.1)), 
          cols = 2)

cluster <- cutree(clHierarch, k = 6)
## number of planets in each cluster
table(unlist((factor(cluster))))
## plot clustered data
x <- planetsData[, 1]
y <- planetsData[, 2]
Clusters <- unlist(factor(cluster))
pd2 <- ggplot(planetsData, aes(x, y, colour = Clusters)) + 
  geom_point(size = 3) +
  ggtitle("Непараметрична кластерізація на 6 груп за методом Варда. Дані очищені від викідів") +
  labs(x = "Радіус планети (в радіусах Юпітера)", y = "Маса планети/M*sin(i) (в масах Юпітера)") 
pd2


planetsData[, "cluster"] <- cluster
planetsData[, "mass_gr"] <- (planetsData$pl_bmassj * 1.898E+30)
planetsData[, "radius_cm"] <- planetsData$pl_radj * 7.149E+09
planetsData[, "volume_cm3"] <- (planetsData$radius_cm)^3*4/3*pi
planetsData[, "density"] <- planetsData$mass_gr/planetsData$volume_cm3
write.csv(planetsData, file = "massRadiusPlanetsData2.csv")

print(clusterSummary <- summarise(group_by(planetsData, cluster), 
                                  count=n(), 
                                  min_mass = min(pl_bmassj), 
                                  max_mass = max(pl_bmassj),
                                  median_mass = median(pl_bmassj),
                                  min_radius = min(pl_radj), 
                                  max_radius = max(pl_radj),
                                  median_radius = median(pl_radj),
                                  min_density = min(density), 
                                  max_density = max(density),
                                  median_density = median(density)))
write.csv(clusterSummary, file = "massRadiusClusterSummary2.csv")



earthLike <- filter(planets, pl_bmasse <= 3, pl_rade <= 1.25)
summarise(earthLike, count = n())
select(earthLike, rowid, pl_name, pl_bmasse, pl_rade)

# ## K-Means Clustering
# clKmeans <- kmeans(planetsData, 5)
# cluster <- clKmeans$cluster
# ## number of planets in each cluster
# table(unlist((factor(cluster))))
# ggplot(planetsData, aes(x, y)) + 
#   geom_point(aes(colour = cluster), size = 3) +
#   ggtitle("K-Means Clustering") +
#   labs(x = "Planet Radius [Jupiter radii]", y = "Planet Mass or M*sin(i) [Jupiter mass]")

planets %>% correlation("st_teff:st_rad") %>% plot
planets %>% correlation("pl_rvamp:pl_insol") %>% plot
# for (i in colnames(planets$st_teff:planets$pl_insol)) {
  i <- planets$pl_eqt
  plot(density(i, na.rm = TRUE))
  shapiro.test(i)
  qqnorm(i);qqline(i, col = 2)
# }

  ## OLS
  result <- regress(dataset = "raw_confirmed", rvar = "st_teff", evar = c("st_mass", "pl_eqt"))
  summary(result, sum_check = c("rmse", "confint"))
  plot(result, plots = "dashboard", lines = "line")
  pred <- predict(result, pred_data = "raw_confirmed")
  print(pred, n = 10)
  store(pred, data = 'raw_confirmed', name = 'predict_reg')
  # write.csv(pred, file = '~/reg_predictions.csv', row.names = FALSE)
  
  
radiant()











