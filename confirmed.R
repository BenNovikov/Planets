# install.packages("ggplot2")
library(ggplot2)
library(plyr)

# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)
# setting parameter names and cleaning leading spaces
columnSet <- function (columns) {
  columns <- rbind(as.data.frame(t(c("rowid", "Row ID"))), columns)
  colnames(columns) <- c("ColumnID", "Parameter")
  columns$Parameter <- trim.leading(columns$Parameter)
  return(columns)
}

# reading data
#setwd("~/Desktop/SofieRFiles")
setwd("F:/SofieRFiles")
planets <- read.csv("raw_confirmed.csv")
pKepler <- read.csv("raw_kepler.csv")
pK2 <- read.csv("raw_k2.csv")
colPlanets <- columnSet(read.csv("raw_conf_col.csv", FALSE, sep=":"))
colKepler <- columnSet(read.csv("raw_kepler_col.csv", FALSE, sep=":"))
colK2 <- columnSet(read.csv("raw_k2_col.csv", FALSE, sep=":"))

# leglist.x <- levels(planets$pl_discmethod)
# lablist.x <- gsub(" ", " \n", leglist.x)
# barplot(table(unlist(planets$pl_discmethod)), 
#         beside = TRUE,
#         col = rainbow(10),
#         legend.text = leglist.x,
#         args.legend = list(x = "top"), 
#         cex.names = .7)  


Methods <- planets$pl_discmethod
ggplot(planets, aes(x = "", labels = NULL, fill = Methods)) +
  geom_bar(width = 1) + 
  coord_polar("y") + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_fill_grey()

Inclination <- planets$pl_orbincl
ggplot(planets, aes(x = Inclination)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") + # Overlay with transparent density plot
  geom_vline(aes(xintercept=mean(Inclination, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)

ggplot(planets, aes(x = Methods, y = Inclination, color = Inclination)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line 
#  (by default includes 95% confidence region)

#tabNumberByMethods <- table(planets$pl_pnum, planets$pl_discmethod)
#tabNumberByMethods <- table(planets$pl_discmethod, planets$pl_pnum)

## http://flowingdata.com/2010/11/23/how-to-make-bubble-charts
radius <- sqrt(planets$pl_radj/pi)
s <- symbols(planets$pl_discmethod, planets$pl_pnum, circles = radius, inches = 0.25, fg = "white", 
        bg = "red", xaxt = 'n', 
        main = "Number of Planets in System Discovered to their Radius", 
        xlab = "",
        ylab = "Number of Planets in System")
#axis(1, at = 1:10, labels = FALSE)
text(x = 1:10, par("usr")[3] - 0.2, labels = lablist.x, srt = 0, pos = 1, xpd = TRUE)
#s + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

hist(planets$pl_bmassj,
     main='Confirmed Planets by their Mass',
     xlab = 'Planet Mass or M*sin(i)[Jupiter mass]',
     ylab = 'Frequency',
     col = rainbow(6))

methodDistanceToInclination <- function(method) {
  planetsByMethod <- planets[planets$pl_discmethod == method, ]
  Distance <- planetsByMethod$st_dist #Distance [pc]
  Inclination <- planetsByMethod$pl_orbincl
  plotTitle <- as.character(paste(method, "Method. Orbit Inclination To Distance"))
  print(c(plotTitle, length(as.vector(table(unlist(planetsByMethod$rowid))))))
  result <- ggplot(planetsByMethod, aes(x = Distance, y = Inclination, color = Inclination)) + 
    scale_fill_discrete(name="Distance") +
    geom_point(colour = "red", size = 5, shape = 1, stroke = 2, show.legend = F) +
    scale_x_log10() +
    ggtitle(plotTitle) +
    labs(x = "Distance to Planet [pc]", y = "Inclination [deg]") +
    theme(plot.title = element_text(color="#666666", face = "bold", size=12, hjust = 0))
  result
}
allDiscMethods <- as.vector(sort(as.character(unique(planets$pl_discmethod, ordered = F)), decreasing = F))
methodDistanceToInclination(allDiscMethods[1])
methodDistanceToInclination(allDiscMethods[2])
methodDistanceToInclination(allDiscMethods[3])
methodDistanceToInclination(allDiscMethods[4])
methodDistanceToInclination(allDiscMethods[5])
methodDistanceToInclination(allDiscMethods[6])
methodDistanceToInclination(allDiscMethods[7])
methodDistanceToInclination(allDiscMethods[8])
methodDistanceToInclination(allDiscMethods[9])
methodDistanceToInclination(allDiscMethods[10])

methodRadiusToDistance <- function(method) {
  planetsByMethod <- planets[planets$pl_discmethod == method, ]
  Distance <- planetsByMethod$st_dist #Distance [pc]
  planetRadius <- planetsByMethod$pl_radj
  plotTitle <- as.character(paste(method, "Method. Planet Radius To Distance"))
  print(c(plotTitle, length(as.vector(table(unlist(planetsByMethod$rowid))))))
  result <- ggplot(planetsByMethod, 
                   aes(x = planetRadius, 
                       y = Distance, 
                       color = Distance)) + 
    scale_fill_discrete(name="Distance") +
    geom_point(colour = "blue", 
               size = planetRadius*2, 
               shape = 1, 
               stroke = 1,
               show.legend = F) + 
    ggtitle(plotTitle) +
    labs(x = "Planet Radius [Jupiter radii]", y = "Distance [pc]") +
    theme(plot.title = element_text(color="#666666", face = "bold", size=12, hjust = 0))
  result
}
methodRadiusToDistance(allDiscMethods[1])
methodRadiusToDistance(allDiscMethods[2])
methodRadiusToDistance(allDiscMethods[3])
methodRadiusToDistance(allDiscMethods[4])
methodRadiusToDistance(allDiscMethods[5])
methodRadiusToDistance(allDiscMethods[6])
methodRadiusToDistance(allDiscMethods[7])
methodRadiusToDistance(allDiscMethods[8])
methodRadiusToDistance(allDiscMethods[9])
methodRadiusToDistance(allDiscMethods[10])
