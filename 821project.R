# 821 project
# Authors: Issy Anand, Julien Bessette

library(qrmtools)
library(MASS)
library(tseries)
library(forecast)

#PROBLEM 1
data <- read.csv('_GSPC.csv')
data <- data$Adj.Close
rets <- returns(data)
rets <- rets[2:length(rets)]

#rets <- c()
#for (i in 1:(length(data)-1)){
#  return <- log(data[i+1]/data[i])
  #rets <- c(rets2, c(return))
#}

fit <- fitdistr(rets, 'Normal')
para <- fit$estimate
hist(rets, breaks = 30, main='Histogram of log returns', xlab='log returns', prob = TRUE)
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
#Empirical distribution clearly not normal, frequencies still high around the tails of the fitted distribution
qqnorm(rets, main = "Normal Q-Q plot of the S&P500 log returns")
qqline(rets)
# data have more extreme values that would be expected if they care from normal distribution 
# it has heavier tails

#PROBLEM 2
#a) Using the given data as portfolio formation period, will use data downloaded after final date as testing period
#b)
rets_ts <- ts(rets, frequency = 365)
ar1 <- arma(rets_ts, order=c(1,0)) # we use the arma function to estimate the paramtere of the AR. the arima function would result in slightly different results
par(mfrow=c(1,2))
Acf(rets_ts, main = "ACF of the log-returns") # for an AR(1) model, PACF should drop close to 0 for lags 2, 3, ...
# doesn't look really significant in this case
Pacf(rets_ts, main = "PACF of the log-returns")
#lag1 autocorrelation of an AR(1) model = coefficient, therefore there is not significant negative lag-1 autocorrelation
# from the two previous graphs, w can't be sure that an order of 1 fits well the data

#c) Very confused about this part, couldnt find the slide he was refering to
#Is the continuous time mean reverting model the Ornstein-Uhlenbeck model?
#If so, to extract this model from an AR(1) model with intercept a, coefficient b:
#OU model: dxt = theta*(mu-xt)*dt + sigma*dWt
#where mu = -a/b, theta = -b
mu = -ar1$coef[2]/ar1$coef[1]
theta = -ar1$coef[1]
sigma = sd(ar1$residuals[2:22])


#PROBLEM 3
plot(rets_ts)
abline(h = mu, col = "red")
abline(h = mu + sigma, col = "blue")
abline(h = mu - sigma, col = "blue")
abline(h = 11 * mu / 10, lty = 2)
abline(h = 9 * mu / 10, lty = 2)




#sim <- arima.sim(n = 5000, list(ar = c(0.3)), sd = sqrt(vari), mean = 0)
#Acf(sim)
#Pacf(sim)

#sim <- arima.sim(n = 5000, list(ma = c(0.7)), sd = sqrt(vari), mean = 0)
#Acf(sim)
#Pacf(sim)

create_df_BB <- function(returns, window, no_sd){
  mov_avg <- rollapply(returns, window, mean)
  mov_sd <- rollapply(returns, window, sd)
  up_bound <- mov_avg + no_sd * mov_sd
  low_bound <- mov_avg - no_sd * mov_sd
  df <- data.frame("Rolling Average" = mov_avg, "Upper bound" = up_bound,
                   "Lower bound" = low_bound)
  return(df)
}




