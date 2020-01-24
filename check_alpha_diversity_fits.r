## #################################
## Check lines on rib panel.

## Set up x values to evaluate.
x <- seq(525, 5000, by=10)

## Check line for Rice Rivers points on rib panel.
yrr <- 31 + (116*x) - (27.5*x*x)

## Check line for Henley Lake points on rib panel.
yhl <- 14.2 + (28.4*x) - (9.16*x*x) + (0.253*x*x*x)


par(mfrow=c(2,2))  ## 2 rows, 2 cols of plots
plot(x, yrr, type="l", xlab="ADD", ylab="Diversity", main="Rice Rivers - Rib")
plot(x, yhl, type="l", xlab="ADD", ylab="Diversity", main="Henley Lake - Rib")
## #################################



## #################################
## Check lines on scapula panel.

## Set up x values to evaluate.
x <- seq(525, 5000, by=10)

## Check line for Rice Rivers points on rib panel.
yrr <- 53.7 + (105*x) - (42.3*x*x) + (49.3*x*x*x) - (34.7*x*x*x*x)

## Check line for Henley Lake points on rib panel.
yhl <- 7.88 + (0.00665*x)


plot(x, yrr, xlab="ADD", ylab="Diversity", main="Rice Rivers - scapula")
plot(x, yhl, xlab="ADD", ylab="Diversity", main="Henley Lake - scapula")
## #################################
