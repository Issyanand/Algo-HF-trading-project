# 821 project
# Authors: Issy Anand, Julien Bessette

#PROBLEM 1
data <- read.csv('_GSPC.csv')
data <- data$Adj.Close
rets <- c()
for (i in 1:(length(data)-1)){
  return <- log(data[i+1]/data[i])
  rets <- c(rets, c(return))
}
library(MASS)
fit <- fitdistr(rets, 'Normal')
para <- fit$estimate
hist(rets, breaks=21, main='Histogram of log returns', xlab='log returns', prob = TRUE)
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
#Empirical distribution clearly not normal, frequencies still high around the tails of the fitted distribution

#PROBLEM 2
#a) Using the given data as portfolio formation period, will use data downloaded after final date as testing period
#b)
require(tseries)
rets_ts <- ts(rets, frequency = 365)
ar1 <- arma(rets_ts, order=c(1,0))
#lag1 autocorrelation of an AR(1) model = coefficient, therefore there is not significant negative lag-1 autocorrelation

#c) Very confused about this part, couldnt find the slide he was refering to
#Is the continuous time mean reverting model the Ornstein-Uhlenbeck model?
#If so, to extract this model from an AR(1) model with intercept a, coefficient b:
#OU model: dxt = theta*(mu-xt)*dt + sigma*dWt
#where mu = -a/b, theta = -b
mu = -ar1$coef[2]/ar1$coef[1]
theta = -ar1$coef[1]
sigma = sd(ar1$residuals[2:22])
