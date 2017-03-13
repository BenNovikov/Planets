Sys.setlocale('LC_ALL', 'ukrainian')
 # install.packages("radiant", repos = "https://radiant-rstats.github.io/minicran/") #rq
 # install.packages("ggplot2")
 # install.packages("ggdendro")
 # install.packages("dplyr")
 # install.packages("tidyr")
 # install.packages("reshape2")
library(radiant)
library(ggplot2)
library(ggdendro)
library(dplyr)
library(tidyr)
library(cluster)
library(reshape2)
library(scales)

##
## constants and functions
##
breaks = 10**(-5:30)            ## log scale breaks
densIron   <- 7.874             ## gr/sm3
densSolid  <- 2
masSuper   <- 10 ## max super Earth mass (Earth mass)

G <- 6.67 * 10^(-11)            ## gravitational constant

cMassEarth    <- 5.9723 * 10^27 ## gr 
cRadiusEarth  <- 6.3780 * 10^8  ## cm
cVolumeEarth  <- 1.0832 * 10^27 ## cm3
cDensityEarth <- 5.5135         ## gr/cm-3

cMassNept     <- 1.0241 * 10^29
cRadiusNept   <- 2.4764 * 10^09
cVolumeNept   <- 6.2540 * 10^28
cDensityNept  <- 1.6376

cMassSat    <- 5.6834 * 10^29
cRadiusSat  <- 6.0268 * 10^09
cVolumeSat  <- 8.2713 * 10^29
cDensitySat <- 0.687

cMassJup      <- 1.89819* 10^30
cRadiusJup    <- 7.1492 * 10^09
cVolumeJup    <- 1.4313 * 10^30
cDensityJup   <- 1.3262

cMassJup2Earth    <- cMassJup/cMassEarth      #317.8323    
cRadiusJup2Earth  <- cRadiusJup/cRadiusEarth  #11.2092
cMassSat2Earth    <- cMassSat/cMassEarth      #95.1627 
cRadiusSat2Earth  <- cRadiusSat/cRadiusEarth  #9.4494 
cMassNept2Earth   <- cMassNept/cMassEarth     #17.148 
cRadiusNept2Earth <- cRadiusNept/cRadiusEarth #3.8827 

cols_bw <- c("#771C19","#AA3929","#E25033","#F27314","#F8A31B",
             "#E2C59F","#B6C5CC","#8E9CA3","#556670","#000000")

myColors <- function() scale_colour_grey()

cUnitsJup <- data.frame(name = 'Земля', 
                        mass    = cMassEarth/cMassJup, 
                        radius  = cRadiusEarth/cRadiusJup,
                        volume  = cVolumeEarth/cVolumeJup,
                        density = cDensityEarth)
cUnitsJup <- rbind(cUnitsJup, c("Нептун", 
                                cMassNept/cMassJup, 
                                cRadiusNept/cRadiusJup,
                                cVolumeNept/cVolumeJup,
                                density = cDensityNept))
cUnitsJup <- rbind(cUnitsJup, c("Юпітер", 1, 1, 1, density = cDensityJup))
cUnitsJup <- rbind(cUnitsJup, c("Сатурн", 
                                cMassSat/cMassJup, 
                                cRadiusSat/cRadiusJup,
                                cVolumeSat/cVolumeJup,
                                density = cDensitySat))

## returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

## setting parameter names and cleaning leading spaces
columnSet <- function (columns) 
  {
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
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) 
  {
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

removeOutliers <- function(x, na.rm = TRUE, ...) 
  {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

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

##
##  end of functions
##

## reading data
#setwd("~/Desktop/SofieRFiles")
setwd("F:/SofieRFiles")
planets <- read.csv("raw_confirmed.csv")
# kepler <- read.csv("raw_kepler.csv")
# k2 <- read.csv("raw_k2.csv")
colPlanets <- columnSet(read.csv("raw_conf_col.csv", FALSE, sep=":"))
# colKepler <- columnSet(read.csv("raw_kepler_col.csv", FALSE, sep=":"))
# colK2 <- columnSet(read.csv("raw_k2_col.csv", FALSE, sep=":"))
#write.csv(columns, file = "colPlanets.csv")

## get list of unique discovery methods 
discoveryMethods.x <- sort(as.character(unique(planets$pl_discmethod)), decreasing = FALSE)
discoveryMethods.x
discoveryMethods.ua <- c("астрометрія", 
                 "зміна термінів затемнення",
                 "прямі спостереження",
                 "мікролінзування",
                 "модуляція орбітального блиску",
                 "спостереження за пульсарами",
                 "зміна термінів пульсації",
                 "доплерівська спектроскопія", 
                 "транзитний метод", 
                 "зміна термінів транзиту")

# split cumulative dataframe by discovery methods
discoveryMethods.filename <- gsub(" ", "", paste0("planets", discoveryMethods.x, ".csv"))
discoveryMethods.filename
splittedPlanetsDF <- list()
for (i in 1:length(discoveryMethods.x)) { 
  splittedPlanetsDF[[i]] <- planets[planets$pl_discmethod == discoveryMethods.x[i], ]
  write.csv(splittedPlanetsDF[[i]], file = discoveryMethods.filename[i])
}

# number of discovered planets til 2016 and 2013
planetsData <- dplyr::arrange(dplyr::select(planets, pl_disc, pl_discmethod), pl_disc)
planetsData <- plyr::count(planetsData, c('pl_disc', 'pl_discmethod'))
pld2016 <- ggplot(planetsData, aes(x = pl_disc, y = freq, fill = pl_discmethod)) +
  geom_bar(stat = "identity") +
  ggtitle("Кількість відкритих планет. 1989-2016.") +
  labs(x = "рік відкриття", y = "кількість відкритих екзопланет") + 
  scale_x_continuous(breaks = seq(min(planetsData$pl_disc),max(planetsData$pl_disc),9)) + 
  scale_y_continuous(breaks = seq(0,1500,500), limits = c(0, 1550)) 

planetsData <- planetsData[planetsData$pl_disc < 2014, ]
pld2013 <- ggplot(planetsData, aes(x = pl_disc, y = freq, fill = pl_discmethod)) + 
  scale_fill_manual(values = cols_bw) + 
  geom_bar(stat = "identity") +
  ggtitle("Кількість відкритих планет. 1989-2013.") +
  labs(x = "рік відкриття", y = "кількість відкритих екзопланет") + 
  scale_x_continuous(breaks = seq(min(planetsData$pl_disc), max(planetsData$pl_disc), 8)) + 
  scale_y_continuous(breaks = seq(0, 150, 50), limits = c(0, 155))

multiplot(pld2013 + 
            scale_fill_manual(values = cols_bw, 
                              name="методи:",
                                breaks=discoveryMethods.x,
                                labels=discoveryMethods.ua) +
            theme(legend.position = c(0.2,0.6), 
                  legend.justification = c(0.1,0.6), 
                  legend.background = element_rect(fill="transparent")),
          pld2016 + 
            theme(legend.position="none"),
          cols=2) 

## Number of planets discovered with detection methods
## update planetsData
planetsData <- dplyr::arrange(dplyr::select(planets, pl_disc, pl_discmethod), pl_disc)
planetsData<- plyr::count(planetsData, c('pl_disc', 'pl_discmethod'))
discoveryMethods.count <- c(NULL)
for(i in 1:length(discoveryMethods.x)) {
  discoveryMethods.count[i] <- sum(planetsData[planetsData$pl_discmethod == discoveryMethods.x[i],]$freq)
}
# discoveryMethods.data <- as.data.frame(discoveryMethods.ua)
# colnames(discoveryMethods.data) <- 'method'
# discoveryMethods.data[, 'count'] <- (discoveryMethods.count)


## scale_y_continuous("кількість відкритих екзопланет", trans = log10_trans(), breaks = breaks)
pld01 <- ggplot(planets, aes(pl_discmethod)) + 
  geom_bar(width = .5) +
  scale_x_discrete("методи пошуку екзопланет", 
                   labels = gsub(" ", "\n", discoveryMethods.ua),
                   position = "bottom") + 
  scale_y_continuous("кількість відкритих екзопланет", 
                     trans = log10_trans(), 
                     breaks = breaks) + 
  annotation_logticks(sides = "l") +
  geom_text(stat='count', aes(label=..count..), vjust = -.3, size = 5) +
  ggtitle("Методи пошуку і кількість відкритих екзопланет") + 
  scale_colour_grey()
pld01

tempSum <- sum(planetsData$freq)
temp8 <- sum(planetsData[planetsData$pl_discmethod == discoveryMethods.x[8],]$freq)
temp9 <- sum(planetsData[planetsData$pl_discmethod == discoveryMethods.x[9],]$freq)
lbl8 <- paste0(discoveryMethods.ua[8],": ", round(temp8/tempSum*100, 2), "%")
lbl9 <- paste0(discoveryMethods.ua[9],": ", round(temp9/tempSum*100, 2), "%")
lbl0 <- paste0("інші методи: ", round((tempSum - temp8 - temp9)/tempSum*100, 2), "%")

pld02 <- ggplot(planets, aes(x = "", labels = NULL, fill = planets$pl_discmethod)) +
  geom_bar(width = 1) + 
  coord_polar("y") + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none") +   
  annotate("text", label = lbl0, x = 1.3, y = 0, size = 6, colour = "black") +
  annotate("text", label = lbl8, x = .9, y = 0, size = 6, colour = "black") +
  annotate("text", label = lbl9, x = .5, y = 0, size = 6, colour = "black") +
  scale_fill_grey() +
  ggtitle("Кількість відкритих планет. 1989-2016.") +
  theme(plot.title = element_text(hjust = 1))
pld02

lbls <- seq(1, max(planets$pl_pnum, na.rm = TRUE), by = 1)
pld03 <- ggplot(planets, aes(pl_pnum)) + 
  geom_bar(width = .5) +
  scale_x_discrete("кількість планет в системі", 
                   breaks = lbls, 
                   limits = lbls,
                   labels = lbls, 
                   position = "bottom") + 
  scale_y_continuous("кількість відкритих екзопланет", 
                     trans = log10_trans(), 
                     breaks = breaks) + 
  annotation_logticks(sides = "l") +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -.3, size = 5) +
  ggtitle("Розподіл за кількістю планет в зоряних системах")
pld03

multiplot(pld01 + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)),
          pld03 + theme(plot.margin = unit(c(.2, 1, 2.6, .5), "cm")),
          cols=2)

##
## Mass & Radius
##

## Add density column to planets data
planets[, "pl_density"] <- ((planets$pl_bmassj * cMassJup)/((planets$pl_radj * cRadiusJup)^3*4/3*pi))

## check the correlation 
result <- correlation(dataset = planets, c("pl_radj", "pl_bmassj", "pl_density"), method = "pearson", data_filter = "")
summary(result)
plot(result)

## Have a look at the densities
# ggplot(d, aes(x = x, y = y)) + geom_point() + 
#   coord_trans(y = "log10", x = "log10") + 
#   scale_y_continuous(tran = log10_trans(), 
#                      breaks = trans_breaks("log10", function(x) 10^x), 
#                      labels = trans_format("log10", math_format(10^.x))) + 
#   scale_x_continuous(trans = log10_trans(),
#                      breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x)))
##====
## all known planets masses and radii
print(c("Planets with mass", as.numeric(summarise(planets[!is.na(planets$pl_bmassj),], count = n()))))
p1 <- ggplot(planets[!is.na(planets$pl_bmassj),], aes(x = pl_bmassj)) +
  geom_histogram(colour="black", binwidth=.2) + 
  ggtitle("Відомі маси екзопланет") +
  labs(x = "Маса планети (в масах Юпітера)", y = "Кількість спостережень")  

print(c("Planets with radius", as.numeric(summarise(planets[!is.na(planets$pl_radj),], count = n()))))
p2 <- ggplot(planets[!is.na(planets$pl_radj),], aes(x = pl_radj)) +
  geom_histogram(colour="black", binwidth=.05) + 
  ggtitle("Відомі радіуси екзопланет") +
  labs(x = "Радіус планети (в радіусах Юпітера)", y = "Кількість спостережень")

p3 <- ggplot(planets[!is.na(planets$pl_bmasse),], aes(x = pl_bmasse)) +
  geom_histogram(colour="black", binwidth=.075) + 
  ggtitle("Відомі маси екзопланет") +
  labs(x = "Маса планети (в масах Землі, логарифмічна шкала)", y = "Кількість спостережень") + 
  scale_x_log10(breaks = breaks, labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(sides = "b")

p4 <- ggplot(planets[!is.na(planets$pl_rade),], aes(x = pl_rade)) +
  geom_histogram(colour="black", binwidth=.025) + 
  ggtitle("Відомі радіуси екзопланет") +
  labs(x = "Радіус планети (в радіусах Землі, логарифмічна шкала)", y = "Кількість спостережень") + 
  scale_x_log10(breaks = breaks, labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(sides = "b")

multiplot(p1, p2, p3, p4, cols=2)

##====
print(paste("Planets with density", as.numeric(summarise(planets[!is.na(planets$pl_density),], count = n()))))
p5 <- ggplot(planets[!is.na(planets$pl_density),], aes(x = pl_density)) +
  geom_histogram(colour="black", binwidth=.1) + 
  ggtitle("Обчислені значення густини екзопланет") +
  labs(x = "Густина планети (г/см3, логарифмічна шкала)", y = "Кількість спостережень") 
# p5 + scale_x_log10()
p6 <- ggplot(planets[!is.na(planets$pl_density),], aes(x = pl_density)) +
  geom_density(colour="black", alpha=.3) +   
  ggtitle("Обчислені значення густини екзопланет") +
  labs(x = "Густина планети (г/см3)", y = "Щільність розподілу") 
# p6
p61 <- ggplot(planets[!is.na(planets$pl_density),], aes(x = pl_density)) +
  geom_histogram(colour="black", binwidth=.1) +   
  ggtitle("Обчислені значення густини екзопланет") +
  labs(x = "Густина планети (г/см3)", y = "Кількість спостережень") 
p61

multiplot(p61 , 
          p5 + scale_x_log10(), cols = 2)

# g	=	G * M / r2
# =	(6.67 * 10-11) * (5.98 * 1024) / (6.378 * 106)2
# =	9.81 meters/second2
earthGravity <- G * cMassEarth / cRadiusEarth ^ 2
earthGravity  * 10.5 / 2 ^ 2
earthDensity <- cMassEarth / (cRadiusEarth ^ 3 * 4 / 3 * pi)

earth <- data.frame(mass = 1, density = earthDensity)
p7 <- ggplot(planets[!is.na(planets$pl_density),], aes(x = pl_density, y = pl_bmasse)) + 
  geom_point(size = 3, color = "#b6c5cc") +
  labs(title="Співвідношення густини та маси планет", 
       subtitle = "вертикальна лінія - густина заліза, горизонтальна - верхня границя маси \"суперземель\" ") + # ggtitle("") +
  labs(x = "Густина планети (г/см3)", y = "Маса планети/M*sin(i) (в масах Землі)") + 
  scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x)), 
                limits = c(1E-3, max(planets$pl_density, na.rm = TRUE))) + 
  scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x)), 
                limits = c(1E-2, max(planets$pl_bmasse, na.rm = TRUE))) + 
  annotation_logticks(sides = "bltr") +
  geom_hline(yintercept = masSuper, colour="black") +
  geom_vline(xintercept = densIron, colour="black")+ 
  geom_point(data = cUnitsJup[1:4,], 
             mapping = aes(x = as.numeric(density),y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "black", size = 3, shape = 19) +
  geom_text(data = cUnitsJup[1:4,], 
            aes(x = as.numeric(density), 
                y = as.numeric(mass) * cMassJup / cMassEarth, 
                label = name), 
            colour = "black", vjust = -.5, hjust = -.2, size = 5)

# geom_point(data = earth, aes(density, mass), colour="green", size=4) + 
p7 

head(sort(planets$pl_density, decreasing = T, na.last = NA), n = 50)

## Shapiro-Wilk test for normality
shapiro.test(planets$pl_bmassj)
shapiro.test(planets$pl_radj)
shapiro.test(planets$pl_density)
## Plot using a qqplot
# qqnorm(planets$pl_bmassj);qqline(planets$pl_bmassj, col = 2)
# qqnorm(planets$pl_radj);qqline(planets$pl_radj, col = 2)
# qqnorm(planets$pl_density);qqline(planets$pl_density, col = 2)
## QQPlot compare to normal distribution
q1 <- qqplot.data(planets[!is.na(planets$pl_bmassj), ]$pl_bmassj)
q2 <- qqplot.data(planets[!is.na(planets$pl_radj), ]$pl_radj)
q3 <- qqplot.data(planets[!is.na(planets$pl_density), ]$pl_density)
multiplot(q1 + ggtitle("Відомі маси екзопланет"),  
          q2 + ggtitle("Відомі радіуси екзопланет"),  
          q3 + ggtitle("Обчислена густина екзопланет"), 
          cols = 3)

# ## clear density with boxplot
# planetsData <- na.omit(planets[, c('pl_density', 'pl_bmassj')])
# x <- removeOutliers(planetsData$pl_density)
# par(mfrow = c(1, 2))
# boxplot(planetsData$pl_density)
# boxplot(x)
# shapiro.test(x)
# q4 <- qqplot.data(x)

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

## 
## initial data preparing
##
columnsPlanetsData <- c('rowid', 'pl_radj', 'pl_bmassj')
colAddPlanetsData  <- 
  c('rowid', 'pl_rade', 'pl_bmasse', 'pl_dens', 'pl_orbsmax', 'pl_eqt', 'st_age', 'pl_discmethod')
prepPlanetsData <- function() {
  result <- na.omit(planets[, columnsPlanetsData])
  result <- merge(x = result, y = planets[, colAddPlanetsData], by.x = 'rowid', by.y = 'rowid' )
}
##

## ====
## clear data with boxplot
planetsData <- prepPlanetsData()
x <- removeOutliers(planetsData$pl_radj)
y <- removeOutliers(planetsData$pl_bmassj)
# png()
par(mfrow = c(1, 4))
boxplot(planetsData$pl_radj)
boxplot(x)
boxplot(planetsData$pl_bmassj)
boxplot(y)
mtext(
  "Вхідні дані до і після усунення екстремальних значень (радіуси планет - 2 графіки зліва, маси - справа)",
  side = 3, 
  line = -3, 
  outer = TRUE, 
  cex = 1
)
# dev.off()

c(max(x, na.rm = TRUE),max(y, na.rm = TRUE))
planetsData <- planetsData %>% filter(pl_radj <= 2.2, pl_bmassj <= 3.47)
c(max(planetsData$pl_radj, na.rm = TRUE),max(planetsData$pl_bmassj, na.rm = TRUE))
summarise(planetsData, count = n())
boxplotData <- planetsData
##
##


####
##
## Clustering
##
####
##
##
planetsData <- prepPlanetsData()
## Hierarchical Clustering
## methods: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), 
## "median" (= WPGMC) or "centroid" (= UPGMC).
clHierarch <- hclust(dist(planetsData[, 2:3]), method = "ward.D")
gd1 <- ggdendrogram(clHierarch, labels = FALSE, leaf_labels = FALSE, rotate = TRUE) + 
  labs(title="Попередня ієрархічна кластерізація за Уордом", 
       subtitle = "Лінія виключення найближчих (зліва) зв'язків між групами")
gd1  + geom_hline(aes(yintercept = max(ggplot_build(gd1)$data[[2]]$y) * 0.07))  

cluster <- cutree(clHierarch, k = 5)
## number of planets in each cluster
table(unlist((factor(cluster))))
## plot clustered data
x <- planetsData$pl_radj
y <- planetsData$pl_bmassj
Clusters <- unlist(factor(cluster))
pd1 <- ggplot(planetsData, aes(x, y, colour = Clusters)) + 
  geom_point(size = 3) +
  ggtitle("Попередня кластерізація за методом Уорда. Графік розсіювання") +
  labs(x = "Радіус планети (в радіусах Юпітера)", y = "Маса планети/M*sin(i) (в масах Юпітера)")  + 
  scale_colour_grey()
pd1 
#+ scale_x_log10() + scale_y_log10()
#+ facet_wrap(nrow =2, ~ cluster)
#+ facet_grid(. ~ cluster)

addPlanetsData <- function() {
  mas <- planetsData$pl_bmassj * cMassJup
  rad <- planetsData$pl_radj * cRadiusJup
  vol <- rad^3 * 4/3 * pi
  den <- mas/vol
  
  planetsData[, "cluster"]    <- cluster
  planetsData[, "density"]    <- den
  planetsData[, "mass_gr"]    <- mas
  planetsData[, "radius_cm"]  <- rad
  planetsData[, "volume_cm3"] <- vol

  planetsData
}

planetsData <- addPlanetsData()
write.csv(planetsData, file = "massRadiusPlanetsData.csv")
#clusterData <- as.data.frame(table(planetsData$cluster))
# distinct(planetsData, cluster)

sumPlanetsData <- function() {
  result <- summarise(group_by(planetsData, cluster), 
            count=n(), 
            min_mass = min(pl_bmassj), 
            max_mass = max(pl_bmassj),
            median_mass = median(pl_bmassj),
            min_radius = min(pl_radj), 
            max_radius = max(pl_radj),
            median_radius = median(pl_radj),
            min_density = min(density), 
            max_density = max(density),
            median_density = median(density))
  result
}

print(clusterSummary <- sumPlanetsData())
write.csv(clusterSummary, file = "massRadiusClusterSummary.csv")

## ====
## clear data with clustering 
planetsData <- planetsData %>% filter(cluster != 1)
c(max(planetsData$pl_radj, na.rm = TRUE),max(planetsData$pl_bmassj, na.rm = TRUE))
summarise(planetsData, count = n())
clusterData <- planetsData


##====
## use one of these lines:
planetsData <- clusterData ## uncomment to use clustering cleared data
# planetsData <- boxplotData ## uncomment to use boxplot cleared data

## max radius(Jupiter Radii), max mass(Jupiter Masses)
c(max(planetsData$pl_radj, na.rm = TRUE),max(planetsData$pl_bmassj, na.rm = TRUE))
as.numeric(summarise(planetsData, count = n()))

planetsData <- planetsData[, unique(c(columnsPlanetsData, colAddPlanetsData))]
method = "ward.D"
# method = "ward.D2"
clHierarch <- hclust(dist(planetsData[, 2:3]), method = method)
gd2 <- ggdendrogram(clHierarch, labels = FALSE, leaf_labels = FALSE, rotate = TRUE) + 
  labs(title = paste0("Остаточна кластерізація за методом Уорда"), 
       subtitle = "Лінія відокремлення слабких (справа) зв'язків між групами")
gd2   + geom_hline(aes(yintercept = max(ggplot_build(gd2)$data[[2]]$y) * 0.13))  

multiplot(gd1 + geom_hline(aes(yintercept = max(ggplot_build(gd1)$data[[2]]$y) * 0.07)), 
          gd2 + geom_hline(aes(yintercept = max(ggplot_build(gd2)$data[[2]]$y) * 0.13)), 
          cols = 2) 

cluster <- cutree(clHierarch, k = 4)

## number of planets in each cluster
table(unlist((factor(cluster))))
## plot clustered data
x <- planetsData$pl_radj
y <- planetsData$pl_bmassj

Clusters <- unlist(factor(cluster))
pd2 <- ggplot(planetsData, aes(x, y, colour = Clusters)) + 
  geom_point(size = 3) +
  ggtitle("Кластерізація за методом Уорда. Дані, очищені від викідів") +
  labs(x = "Радіус планети (в радіусах Юпітера)", y = "Маса планети/M*sin(i) (в масах Юпітера)")   

pd2 + 
  geom_text(data = cUnitsJup, 
            aes(x = as.numeric(radius), y = as.numeric(mass), label = name, colour = "black"), 
            colour = "black",
            vjust = -1,
            hjust = .8,
            size = 5) +
  geom_point(data = cUnitsJup, 
             mapping = aes(x = as.numeric(radius), y = as.numeric(mass)),
             colour = cols_bw[7],
             size = 5) + 
  geom_point(data = cUnitsJup, 
             mapping = aes(x = as.numeric(radius), y = as.numeric(mass)),
             colour = "white",
             size = 2) + 
  scale_colour_grey()

pd2 + scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bltr") + 
  geom_text(data = cUnitsJup, 
            aes(x = as.numeric(radius), y = as.numeric(mass), label = name, colour = "black"), 
            colour = "black",
            vjust = 1.1,
            hjust = -.2,
            size = 5) +
  geom_point(data = cUnitsJup, 
             mapping = aes(x = as.numeric(radius), y = as.numeric(mass)),
             colour = cols_bw[7],
             size = 5) + 
  geom_point(data = cUnitsJup, 
             mapping = aes(x = as.numeric(radius), y = as.numeric(mass)),
             colour = "white",
             size = 2) +  
  scale_colour_grey()

planetsData <- addPlanetsData()
write.csv(planetsData, file = "massRadiusPlanetsData2.csv")

print(clusterSummary <- sumPlanetsData())
write.csv(clusterSummary, file = "massRadiusClusterSummary2.csv")

## log relation
##
pd3 <- ggplot(planetsData, aes(x, y, colour = Clusters)) + 
  geom_point(size = 3) +
  ggtitle("Кластеризовані дані без викідів. Логарифмічні шкали") +
  labs(x = "Радіус планети (в радіусах Юпітера)", 
       y = "Маса планети/M*sin(i) (в масах Юпітера)",
       subtitle = "Кластер №1 об'єднує Земле- та Нептуноподібні планети"
  ) + 
  scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides = "bltr") 

pd3 + 
  geom_text(data = cUnitsJup, 
            aes(x = as.numeric(radius), y = as.numeric(mass), label = name, colour = "black"), 
            colour = "black",
            vjust = -1,
            hjust = 1.2,
            size = 5) +
  geom_point(data = cUnitsJup, 
             mapping = aes(x = as.numeric(radius), y = as.numeric(mass)),
             colour = cols_bw[7],
             size = 5) + 
  geom_point(data = cUnitsJup, 
             mapping = aes(x = as.numeric(radius), y = as.numeric(mass)),
             colour = "white",
             size = 2) + 
  scale_colour_grey() 


#####
##
## cluster 1 selected
##
#####
## 
clusterDataAddType <- function(clusterData) {
  clusterData.len <- as.numeric(summarise(clusterData, count = n()))
  clusterData[, "pl_type"] <- vector(length = clusterData.len)
  
  for(i in 1:clusterData.len) {
    # if(is.na(clusterData$pl_dens[i])) {
    #   clusterData$pl_dens[i] = round(clusterData$density[i], 3)
    # }  
    
    if(clusterData$density[i] > densSolid && 
       clusterData$pl_rade[i] <= 2 &&
       clusterData$pl_bmasse[i] <= 3) {
      clusterData$pl_type[i] = "solid"
    } else if(clusterData$pl_rade[i] < 4 &&
              clusterData$pl_bmasse[i] < 20) {
      clusterData$pl_type[i] = "subneptune"
    } else if(clusterData$pl_rade[i] < 6 &&
              clusterData$pl_bmasse[i] < 30) {
      clusterData$pl_type[i] = "neptunian"
    } else {
      clusterData$pl_type[i] = "gas"
    }
  }
  clusterData
}

clusterData <- planetsData[planetsData$cluster == 1,]
clusterData <- clusterDataAddType(clusterData)
summarise(clusterData, count = n())

pd3 <- ggplot(clusterData, aes(x = pl_rade, y = pl_bmasse)) + 
  geom_point(size = 3) +
  ggtitle("Планети першої групи. Графік розсіювання.") +
  labs(x = "Радіус планети (в радіусах Землі)", 
       y = "Маса планети/M*sin(i) (в масах Землі)") + 
  scale_x_continuous() + 
  scale_y_continuous() + 
  geom_point(data = clusterData[clusterData$pl_type == "solid",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = "green",
             size = 2) + 
  geom_point(data = clusterData[clusterData$pl_type == "subneptune",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = "blue",
             size = 2) + 
  geom_point(data = clusterData[clusterData$pl_type == "neptunian",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = "grey",
             size = 2) + 
  geom_point(data = clusterData[clusterData$pl_type == "gas",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = "white",
             size = 2) +
  geom_text(data = cUnitsJup[1:2,], 
            aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                y = as.numeric(mass) * cMassJup / cMassEarth, 
                label = name, 
                colour = "black"), 
            colour = "black",
            vjust = 1,
            hjust = -.2,
            size = 5) +
  geom_point(data = cUnitsJup[1:2,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "black",
             size = 5) + 
  geom_point(data = cUnitsJup[1:2,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "white",
             size = 2)# 

multiplot(pd3 + labs(subtitle = "Умовні кольори планет: зелений - 'Землі', блакитний - субнептунові, 
сірий - транснептунові, білий - транснептунові/газові гіганти" ),
  pd3 + labs(subtitle = "Умовні кольори планет: зелений - 'Землі', блакитний - субнептунові, 
сірий - транснептунові, білий - транснептунові/газові гіганти" ) + 
  scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides = "bltr"),
  cols = 2) 

pd31 <- ggplot(clusterData, aes(x = pl_rade, y = pl_bmasse)) + 
  geom_point(size = 3) +
  ggtitle("Планети першої групи. Графік розсіювання.") +
  labs(x = "Радіус планети (в радіусах Землі)", 
       y = "Маса планети/M*sin(i) (в масах Землі)") + 
  scale_x_continuous() + 
  scale_y_continuous() + 
  geom_point(data = clusterData[clusterData$pl_type == "solid",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = "white",
             size = 2) + 
  geom_point(data = clusterData[clusterData$pl_type == "subneptune",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = cols_bw[2],
             size = 2) + 
  geom_point(data = clusterData[clusterData$pl_type == "neptunian",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = cols_bw[4],
             size = 2) + 
  geom_point(data = clusterData[clusterData$pl_type == "gas",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = cols_bw[6],
             size = 2) +
  geom_text(data = cUnitsJup[1:2,], 
            aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                y = as.numeric(mass) * cMassJup / cMassEarth, 
                label = name, 
                colour = "black"), 
            colour = "black",
            vjust = 1,
            hjust = -.2,
            size = 5) +
  geom_point(data = cUnitsJup[1:2,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "black",
             size = 5) +
  geom_point(data = cUnitsJup[1:2,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "white",
             size = 2) # 


multiplot(pd31 + labs(subtitle = "Умовні 'Землі', суб-, транс-Нептуни, транснептунові/газові гіганти" ), 
          pd31 + labs(subtitle = "Логарифмічна шкала") +
            scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            annotation_logticks(sides = "bltr") + 
            scale_colour_grey(),
          cols = 2) 


####
##
##
##
## Subneptunian planets M-R relation
##
##
##
####

smallPlanets <- clusterData[clusterData$pl_type != "gas", ]
smallPlanets <- smallPlanets[smallPlanets$pl_type != "neptunian", ]
as.numeric(summarise(smallPlanets, count = n()))

# ## clear density with boxplot
# planetsData <- na.omit(planets[, c('pl_density', 'pl_bmassj')])
x <- removeOutliers(smallPlanets$density)
par(mfrow = c(1, 2))
boxplot(smallPlanets$density)
boxplot(x)
shapiro.test(x)
q4 <- qqplot.data(smallPlanets$density)

p8 <- ggplot(smallPlanets[!is.na(smallPlanets$density),], aes(x = density)) +
  geom_histogram(colour="black", binwidth = .5) + 
  ggtitle("Вхідна вибірка густини") +
  labs(x = "Густина планети (г/см3)", y = "Кількість спостережень") 

p9 <- ggplot(smallPlanets[!is.na(x),], aes(x = density)) +
  geom_histogram(colour="black", alpha=.2) +   
  ggtitle("Дані без викідів") +
  labs(x = "Густина планети (г/см3)", y = "Кількість спостережень") 

smallPlanets <- smallPlanets[!is.na(x),]

shapiro.test(smallPlanets$density)
q41 <- qqplot.data(smallPlanets$density)

multiplot(p8, q4, p9, q41, cols = 4)

pd32 <- ggplot(smallPlanets, aes(x = pl_rade, y = pl_bmasse)) + 
  geom_point(size = 3) +
  ggtitle("Планети першої групи. Графік розсіювання.") +
  labs(x = "Радіус планети (в радіусах Землі)", 
       y = "Маса планети/M*sin(i) (в масах Землі)") + 
  scale_x_continuous() + 
  scale_y_continuous() + 
  geom_point(data = smallPlanets[smallPlanets$pl_type == "solid",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = "white",
             size = 2) + 
  geom_point(data = smallPlanets[smallPlanets$pl_type == "subneptune",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = cols_bw[2],
             size = 2) + 
  geom_point(data = smallPlanets[smallPlanets$pl_type == "neptunian",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = cols_bw[4],
             size = 2) + 
  geom_point(data = smallPlanets[smallPlanets$pl_type == "gas",], 
             mapping = aes(x = as.numeric(pl_rade), 
                           y = as.numeric(pl_bmasse)),
             colour = cols_bw[6],
             size = 2) +
  geom_text(data = cUnitsJup[1:2,], 
            aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                y = as.numeric(mass) * cMassJup / cMassEarth, 
                label = name, 
                colour = "black"), 
            colour = "black",
            vjust = -.5,
            hjust = .8,
            size = 5) +
  geom_point(data = cUnitsJup[1:2,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "black",
             size = 5) +
  geom_point(data = cUnitsJup[1:2,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "white",
             size = 2)# 


multiplot(pd32 + labs(subtitle = "Субнептунові планети (R < 4Rз, M < 20Mз)" ), 
          pd32 + labs(subtitle = "Логарифмічна шкала") +
            scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            annotation_logticks(sides = "bltr") + 
            scale_colour_grey(),
          cols = 2) 

as.numeric(summarise(smallPlanets, count = n()))

# List to compare to Wolfgang
smallList <- select(planets, rowid, pl_name, pl_orbper, pl_bmasse, pl_rade, pl_disc, rowupdate, pl_discmethod)
smallList <- smallList[smallList$rowid %in% smallPlanets$rowid, ]
smallList <- smallList[ order(smallList$rowupdate, decreasing=TRUE), ]
write.csv(smallList, file = "smallList.csv")

####
##
## sampling for transneptunes models test
##
####
data <- smallPlanets[smallPlanets$pl_type != "neptunian", ]
# data <- data[data$pl_type != "gas", ]
sampleSize <- as.numeric(summarise(data, count = n()))
sampleSize

####
##
## Models
##
####
x <- round(data$pl_rade, 4)  
y <- round(data$pl_bmasse, 4)
par(mfrow = c(1, 2))
shapiro.test(x)
shapiro.test(y)
q51 <- qqplot.data(x)
q52 <- qqplot.data(y)
pd51 <- ggplot(data, aes(x = pl_rade)) +
  geom_histogram(colour="black", alpha=.2) +   
  ggtitle("Радіуси субнептунів") +
  labs(x = "Радіус (радіусів Землі)", y = "Кількість спостережень") 

pd52 <- ggplot(data, aes(x = pl_bmasse)) +
  geom_histogram(colour="black", alpha=.2) +   
  ggtitle("Маси субнептунів") +
  labs(x = "Маса (мас Землі)", y = "Кількість спостережень") 

multiplot(pd51, q51, pd52, q52, cols = 4)

## R2 = 1 - SSE / SST
## R2adj = 1 - MSE / MST
## MST = SST/(n - 1) 
## SSE = sum((y - f) ^ 2)         # y - observed data, f - predicted value
## SSR = sum((f - mean(y)) ^ 2)   # y - observed data, f - predicted value
## SST = sum((y - mean(y)) ^ 2)   # y - observed data
## MSE = SSE/(n - 2)              # n - sample size, 2 is parameters number
##  R2adjust = 1 - [[(T - 1) / (T - k)](1 - R2)]
R2 <- function(x, y, b1, b2, bi) {
  n = length(x)
  f = (b2 + b1 * x ^ bi)
  SSE = sum((y - f) ^ 2)
# 
# are the variables distributed normally?
#
#   shapiro.test((y - f))
#   par(mfrow = c(1, 3))
#   gr1 <- plot(f, (y - f)) 
#   gr2 <- qqnorm((y - f)) 
#   gr3 <- hist((y - f))
#   
  MSE = SSE/(n - 2)
  SST = sum((y - mean(y)) ^ 2)
  MST = SST/(n - 1)
  SSR = sum((f - mean(y)) ^ 2)
  R2    = 1 - SSE / SST
  if(R2 > 0) R2adj = 1 - MSE / MST else R2adj = 1 - ((n - 1)/(n - 2)*(1 - R2))
  Pval  = SSR/MSE
  c(R2, R2adj, Pval, gr1)
}

####
## Wolfgang et al., 2015. M/Mз = 2.7*(R/Rз)^1.3
x <- pRad <- round(data$pl_rade, 4)   # radius
y <- pMas <- round(data$pl_bmasse, 4) # temp data to fit model
b1 <- 2.7
b2 <- 0
bi  <- 1.3
R2(pRad, pMas, b1, b2, bi)

par(mfrow = c(1, 1))
plot(pRad, pMas, xlab = "", ylab = "")

y <- round(b1 * pRad ^ bi, 4)
SSRWolf <- sum((pMas - y) ^ 2)
SSRWolf
title(paste0("mass ~ ", b2, " + ", b1, " * radius ^ ", round(bi,2),  ". SSE:", round(SSRWolf, 1)),
      xlab = "Радіус планети (в радіусах Землі)", ylab = "Маса планети/M*sin(i) (в масах Землі)")
lines(nepRad <- seq(min(pRad), max(pRad), len = nrow(data)), b1 * nepRad ^ bi, type = "l")

####
## Lissauer et al. (2011) M = R ^ 2.06
x <- pRad <- round(data$pl_rade, 4)   # radius
y <- pMas <- round(data$pl_bmasse, 4) # temp data to fit model
b1 <- 1
b2 <- 0
bi  <- 2.06
R2(pRad, pMas, b1, b2, bi)

par(mfrow = c(1, 1))
plot(pRad, pMas, xlab = "", ylab = "")

y <- round(b1 * pRad ^ bi, 4)
SSRLiss <- sum((pMas - y) ^ 2)
SSRLiss
title(paste0("mass ~ ", b2, " + ", b1, " * radius ^ ", round(bi,2),  ". SSE:", round(SSRLiss, 1)),
      xlab = "Радіус планети (в радіусах Землі)", ylab = "Маса планети/M*sin(i) (в масах Землі)")
lines(nepRad <- seq(min(pRad), max(pRad), len = nrow(data)), b1 * nepRad ^ bi, type = "l")

####
## Wu & Lithwick (2013) M = 3 * R
pRad <- round(data$pl_rade, 4)   # radius
pMas <- round(data$pl_bmasse, 4) # temp data to fit model
b1 <- 3
b2 <- 0
bi  <- 1
R2(pRad, pMas, b1, b2, bi)

par(mfrow = c(1, 1))
plot(pRad, pMas, xlab = "", ylab = "")

y <- round(b1 * pRad ^ bi, 4)
SSRWuLi <- sum((pMas - y) ^ 2)
SSRWuLi
title(paste0("mass ~ ", b2, " + ", b1, " * radius ^ ", round(bi,2),  ". SSE:", round(SSRWuLi, 1)),
      xlab = "Радіус планети (в радіусах Землі)", ylab = "Маса планети/M*sin(i) (в масах Землі)")
lines(nepRad <- seq(min(pRad), max(pRad), len = nrow(data)), b1 * nepRad ^ bi, type = "l")

####
## (Marcy et al. 2014) M = 2,69 * R ^ 0.93
pRad <- round(data$pl_rade, 4)   # radius
pMas <- round(data$pl_bmasse, 4) # temp data to fit model
b1 <- 2.69
b2 <- 0
bi  <- 0.93
R2(pRad, pMas, b1, b2, bi)

par(mfrow = c(1, 1))
plot(pRad, pMas, xlab = "", ylab = "")

y <- round(b1 * pRad ^ bi, 4)
SSRMarcy <- sum((pMas - y) ^ 2)
SSRMarcy
title(paste0("mass ~ ", b2, " + ", b1, " * radius ^ ", round(bi,2),  ". SSE:", round(SSRMarcy, 1)),
      xlab = "Радіус планети (в радіусах Землі)", ylab = "Маса планети/M*sin(i) (в масах Землі)")
lines(nepRad <- seq(min(pRad), max(pRad), len = nrow(data)), b1 * nepRad ^ bi, type = "l")


# + geom_smooth(method = "lm", span = 1.5)

##
## simple least squares regression
##
## test
findModel <- function(lowLim = 0.01, hiLim = 10, stepSize = 0.01) {
  b1 = 1
  b2 = 1
  # par(mfrow = c(2, 4))
  res1 <- 0; res2 <- 0; resi <- 0; rese <- 99999
  for(i in seq(lowLim, hiLim, by = stepSize)) {
    ## model
    model <- nls(massVal ~ b2 + b1 * (radiVal)^i, start = list(b1 = b1, b2 = b2))
    # Plot the chart with new data by fitting it to a prediction from 100 data points.
    new.data <- data.frame(radiVal = seq(min(radiVal), max(radiVal), len = sampleSize))
    # lines(new.data$radiVal, predict(model, newdata = new.data))
    res <- (predict(model, radiVal) - massVal) ^ 2
    # print(sum(resid(model)^2))
    SSRMine <- round(sum(res), 3)
    
    # Plot these values.
    # plot(radiVal, massVal, xlab = "", ylab = "")    
    # title(paste0("radius ~ b2 + b1 * mass ^ ", round(i,2),  ". SSE:", 
    #              print(SSRMine <- round(sum(resid(model)^2), 3))))
    # Get the confidence intervals on the chosen values of the coefficients.
    # print(confint(model))
    # summary(model)
    b1 <- round(sum(confint(model)[1,])/2, 4)
    b2 <- round(sum(confint(model)[2,])/2, 4)
    # title(paste0("R ~ ", b2, " + ", b1, " * M ^ ", round(i,3),  ". SSE:", round(SSRMine, 4)),
    #       xlab = "Радіус планети (в радіусах Землі)", ylab = "Маса планети/M*sin(i) (в масах Землі)")
    if(SSRMine < rese) {
      rese <- SSRMine
      resi <- i
      res1 <- b1
      res2 <- b2
      print(c("Power = ", i, " New SSE: ", SSRMine))
      summary(model)
    }
  }
  result <- c(res1, res2, resi, rese)
  result
}

checkModel <- function(model) {
  b1 <- model[1] 
  b2 <- model[2]   ##intercept
  i  <- model[3]   ##power
  rese <- model[4] ##SSE
  y  <- round(b2 + b1 * radiVal ^ i, 4)
  SSRMine <- sum((massVal - y) ^ 2)
  print(c("Model SSE: ", rese, " SSE: ", SSRMine))
  
  #plot it
  par(mfrow = c(1, 1))
  plot(radiVal, massVal, xlab = "", ylab = "")  
  title(paste0("radius ~ ", b2, " + ", b1, " * mass ^ ", round(i,2),  ". SSE:", round(SSRMine, 1)),
        xlab = "Радіус планети (в радіусах Землі)", ylab = "Маса планети/M*sin(i) (в масах Землі)")
  lines(x <- seq(min(radiVal), max(radiVal), len = nrow(data)), b2 + b1 * x ^ i, type = "l")
  result <- c(b1, b2, i, SSRMine)
}

# model <- findModel(.1, 10, 0.1)
# checkModel(model)
massVal <- data$pl_bmasse
radiVal <- data$pl_rade

model <- findModel(-2, -.01, 0.01)
result <- checkModel(model)
result
b1 <- result[1]
b2 <- result[2]
bi <- result[3]
SSRMine <- result[4]




#checkModel(c(1.3, .7, 0, 0))
#   model <- nls(radiVal ~ b1 * massVal ^ i + b2, start = list(b1 = b1, b2 = b2))
#   # Plot the chart with new data by fitting it to a prediction from 100 data points.
#   new.data <- data.frame(massVal = seq(min(massVal), max(massVal), len = 100))
#   lines(new.data$massVal, predict(model, newdata = new.data))
#   title(paste0("radius ~ b2 + b1 * mass ^ ", round(i,2),  ". SSE:", 
#                print(SSRMine <- round(sum(resid(model)^2), 3))))


pd5 <- ggplot(data, aes(x = pl_rade, y = pl_bmasse)) + 
  geom_point(size = 3) +
  ggtitle(paste0("Плутоноподібні екзопланети. Вибірка: ", sampleSize, " планет")) +
  labs(x = "Радіус планети (в радіусах Землі)", 
       y = "Маса планети/M*sin(i) (в масах Землі)",
       subtitle = paste0("Модель: M/Mз = 3 * (R/Rз) ^ 1. SSE = ", SSRWuLi)) + 
  scale_x_continuous() + 
  scale_y_continuous() + 
  geom_text(data = cUnitsJup[1:2,], 
            aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                y = as.numeric(mass) * cMassJup / cMassEarth, 
                label = name), 
            colour = "black", vjust = -.5, hjust = 1.2, size = 5) +
  geom_point(data = cUnitsJup[1:2,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "blue", size = 3) +
  geom_point(data = cUnitsJup[1:2,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "white", size = 2) 

multiplot(pd5 +
            geom_point(data = data, 
                       mapping = aes(data$pl_rade, 
                                     y = 3 * na.omit(data$pl_rade) ^ 1), 
                       colour = "blue", size = 1), 
          pd5 +
            geom_point(data = data, 
                       mapping = aes(data$pl_rade, 
                                     y = 2.7 * na.omit(data$pl_rade) ^ 1.3), 
                       colour = "blue", size = 1) + 
            scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            annotation_logticks(sides = "bltr"),
          cols = 2)

####
##
## Fit model to clusters 1 & 2
##
####
clusterData <- planetsData[planetsData$cluster < 3,]
clusterData <- clusterDataAddType(clusterData)
summarise(clusterData, count = n())
# data <- smallPlanets[smallPlanets$pl_type != "neptunian", ]
data <- na.omit(clusterData[, c('rowid','pl_rade', 'pl_bmasse', 'cluster')])
as.numeric(summarise(data, count = n()))
sampleSize <- as.numeric(summarise(data, count = n()))

# massVal <- data$pl_bmasse
# radiVal <- data$pl_rade
# 
# model <- findModel(.1, 10, 0.1)
# result <- checkModel(model)
# result
b1 <- 3
b2 <- 0
bi <- 1
SSRMine <- round(SSRWuLi, 2)


# x <- clusterData$pl_rade
# y <- clusterData$pl_bmasse
# model <- lm(y ~ x, data = clusterData)
# model
# summary(model)
# anova(model)

pd5 <- ggplot(data, aes(x = pl_rade, y = pl_bmasse)) + 
  geom_point(size = 3, color = cols_bw[6 + data$cluster]) +
  ggtitle(paste0("Кластери 1-2. Вибірка з ", sampleSize, " планет")) +
  labs(x = "Радіус планети (в радіусах Землі)", 
       y = "Маса планети/M*sin(i) (в масах Землі)",
       subtitle = paste0("Модель: M =",b2, " + ", b1, " * (R/Rз) ^ ", bi, ".")) + 
  scale_x_continuous() + 
  scale_y_continuous() + 
  geom_text(data = cUnitsJup[1:4,], 
            aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                y = as.numeric(mass) * cMassJup / cMassEarth, 
                label = name), 
            colour = "black", vjust = -1, hjust = 1.2, size = 5) +
  geom_point(data = cUnitsJup[1:4,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "blue", size = 3) +
  geom_point(data = cUnitsJup[1:4,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "white", size = 2)
pd5 +
  geom_point(data = data, 
             mapping = aes(x = data$pl_rade, y = b2 + b1 * data$pl_rade ^ bi), 
             colour = "red", size = 1)

multiplot(pd5 +
            geom_point(data = data, 
                       mapping = aes(data$pl_rade, 
                                     y = b2 + b1 * na.omit(data$pl_rade) ^ bi), 
                       colour = "blue", size = 1), 
          pd5 +
            geom_point(data = data, 
                       mapping = aes(data$pl_rade, 
                                     y = b2 + b1 * na.omit(data$pl_rade) ^ bi), 
                       colour = "blue", size = 1)  + 
            scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            annotation_logticks(sides = "bltr"),
          cols = 2)  


#####
##
## Earth like planets
##
#####

data <- smallPlanets[smallPlanets$pl_type == "solid", ]
sampleSize <- as.numeric(summarise(data, count = n()))

plTrans <- planets[planets$pl_discmethod == "Transit", ]
plOBM   <- planets[planets$pl_discmethod == "Orbital Brightness Modulation", ]
mean(plTrans$pl_bmasse, na.omit = T)

####
##
## Models
##
## M = 4/3 * pi * R^3 * [1 + (1 - 3/5 * n)(2/3 * pi * R^2)^n]
## n = 0.513 corresponds to the H2O modified polytrope, 
## n = 0.528 corresponds to Fe, and 
## n = 0.544 corresponds to MgSiO3).
## Seager et al., 2007


pd4 <- ggplot(data, aes(x = pl_rade, y = pl_bmasse)) + 
  geom_point(size = 3) +
  ggtitle(paste0("Планети земної групи. Вибірка: ", sampleSize, " спостережень")) +
  labs(x = "Радіус планети (в радіусах Землі)", 
       y = "Маса планети/M*sin(i) (в масах Землі)",
       subtitle = paste0("Лінійна регресія. 95% дов.інтервал. R2: 0.973, Скоригов. R2: 0.968")) + 
  scale_x_continuous() + 
  scale_y_continuous() + 
  geom_text(data = cUnitsJup[1,], 
            aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                y = as.numeric(mass) * cMassJup / cMassEarth, 
                label = name), 
            colour = "black", vjust = -1, hjust = 1.2, size = 5) +
  geom_point(data = cUnitsJup[1,], 
             mapping = aes(x = as.numeric(radius) * cRadiusJup / cRadiusEarth, 
                           y = as.numeric(mass) * cMassJup / cMassEarth), 
             colour = "green", size = 3) + 
  geom_smooth(method = "lm", span = 1.5)
  # geom_abline(mapping = NULL, data = data, slope = 0.355, intercept = 0.563, na.rm = FALSE, show.legend = NA)
  #  

multiplot(pd4, pd4 + 
            scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
            annotation_logticks(sides = "bltr"),
          cols = 2)

pd4 + 
  scale_x_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = breaks, labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides = "bltr")
### doesn't work cause nonlinear
  ## OLS
par(mfrow = c(2, 3))
result <- regress(dataset = data, rvar = "pl_rade", evar = "pl_bmasse")
summary(result, sum_check = c("rmse", "confint"))
plot(result, plots = "dashboard", lines = "line")
pred <- predict(result, pred_data = data)
print(pred, n = 10)


# radiant()

x <- 0:12
y <- sin(pi/5 * x)
op <- par(mfrow = c(3,3), mar = .1+ c(2,2,3,1))
for (tp in c("p","l","b",  "c","o","h",  "s","S","n")) {
  plot(y ~ x, type = tp, main = paste0("plot(*, type = \"", tp, "\")"))
  if(tp == "S") {
    lines(x, y, type = "s", col = "red", lty = 2)
    mtext("lines(*, type = \"s\", ...)", col = "red", cex = 0.8)
  }
}
par(op)
