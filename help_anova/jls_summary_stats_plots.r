library("tidyverse")
library("readxl")
library("GGally")
library("gridExtra")

environment <- read_excel("HL_RRC_Environmental.xlsx", na="NA")


## Summary stats:
apply(environment %>% select(Temp, SpCond, pH, Sal), 2, summary)


## Correlation matrices for both sites combined:
cor(environment %>% select(Temp, SpCond, pH, Sal), use="complete.obs")
## Correlation matrices for the sites, separately:
## Henley Lake:
cor(environment %>% filter(Site=="HenleyLake") %>% select(Temp, SpCond, pH, Sal), use="complete.obs")
## Rice Rivers:
cor(environment %>% filter(Site=="RiceRivers") %>% select(Temp, SpCond, pH, Sal), use="complete.obs")


## Scatterplots, correlations, and density plots for the 4 response
## variables, by site.  (Uses GGally package.)
ggpairs(environment %>% drop_na(), columns=6:9, aes(color=Site))


## Boxplots of each response variable by site (uses gridExtra package).
plotTemp <- ggplot(environment %>% drop_na(), aes(x=Site, y=Temp, fill=Site)) + geom_boxplot()
plotSpCond <- ggplot(environment %>% drop_na(), aes(x=Site, y=SpCond, fill=Site)) + geom_boxplot()
plotpH <- ggplot(environment %>% drop_na(), aes(x=Site, y=pH, fill=Site)) + geom_boxplot()
plotSal <- ggplot(environment %>% drop_na(), aes(x=Site, y=Sal, fill=Site)) + geom_boxplot()
grid.arrange(plotTemp, plotSpCond, plotpH, plotSal, ncol=2)


## Individual scatterplots of the various pairs of variables, with different
## colors for the 2 sites:
ggplot(environment %>% drop_na(), aes(x=Temp, y=SpCond, color=Site)) + geom_point()
ggplot(environment %>% drop_na(), aes(x=Temp, y=pH, color=Site)) + geom_point()
ggplot(environment %>% drop_na(), aes(x=Temp, y=Sal, color=Site)) + geom_point()
## 3 more pairings possible...
