library(qrmtools)
library(MASS)
library(tseries)
library(forecast)

#PROBLEM 1
#data <- read.csv('_GSPC.csv')
data <- read.csv('MF821_sp500data.csv')
data <- data$Adj.Close
rets <- diff(log(data),lag=1)
plot(data,type='l',ylab = 'S&P500 price',xlab = 'days')
title('S&P500 Price from 02/01/2020 to 02/04/2020')
abline(v=45, col="red", lty=2, lwd=3)

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
time <- c(1:45)/252
time2 <- c(1:65)/252
S_ave <- lm(data[1:45] ~ time)
line1 <- S_ave$coefficients[1]+S_ave$coefficients[2]*time2
lines(line1,col='red')
line1

S_ave$coefficients
#c) 
Y <- data[1:45] - S_ave$fitted.values
delta_Y <- diff(Y,lag=1)
ar1 <- arima(delta_Y, order=c(1,0,0))
#par(mfrow=c(1,2))
Acf(ts(delta_Y), main = "ACF of delta_y")
Pacf(ts(delta_Y), main = "PACF of delta_y")
#The lag one autocorrelation is given bu the ar1 coefficient, in this case there is siginificant negative autocorrelation
#for an AR(1) model the PACF should drop close to zero for lags 2,3,...
#From the Pacf the later lag partial autocorrelations arent overly significant but
#we cant be 100% sure an ar(1) is a good fit for the


#d) Very confused about this part, couldnt find the slide he was refering to
#Is the continuous time mean reverting model the Ornstein-Uhlenbeck model?
#If so, to extract this model from an AR(1) model with intercept a, coefficient b:
#OU model: dxt = kappa*(theta-xt)*dt + sigma*dWt
#where theta = -a/b, kappa = -b
kappa = as.numeric(1-ar1$coef[1])
#Is theta right?
theta = as.numeric((1/kappa)*(mean(delta_Y+kappa*Y[1:(length(Y)-1)])))
sigma = sd(Y)

par(mfrow=c(1,1))
#PROBLEM 3
library(zoo)
preds = predict.lm(S_ave,data.frame('time' = c(1:length(data[1:45])) / 252))
rolling_preds = rollapply(preds,3, mean)
Y_u = rolling_preds + 3*sigma
Y_l = rolling_preds - 3*sigma
plot(data[1:length(data[1:45])], type = 'l')
lines(Y_u, col = "red")
lines(Y_l, col = "blue")

create_df_BB <- function(rets, window, no_sd, inner_sd){
  moving_avg <- rollapply(rets, window, mean)
  moving_sd <- rollapply(rets, window, sd)
  up_bound <- moving_avg + no_sd * moving_sd
  low_bound <- moving_avg - no_sd * moving_sd
  #adding the inner bands where we close our position
  inner_up_bound <- moving_avg + inner_sd * moving_sd
  inner_low_bound <- moving_avg - inner_sd * moving_sd
  df <- data.frame("Rolling_Ave" = moving_avg, "Upper_bound" = up_bound,
                   "Lower_bound" = low_bound, "Inner_upper_bound" = inner_up_bound, 
                   "Inner_lower_bound" = inner_low_bound)
  return(df)
}
bb = create_df_BB(data[45:67], 3, 0.75, 0.1)
plot(data[45:65], type = 'l')
lines(bb$Lower_bound, col = 'red')
lines(bb$Upper_bound, col = 'red')

lines(bb$Rolling_Ave, col  = 'blue')
#inner bands
lines(bb$Inner_upper_bound, col  = 'green')
lines(bb$Inner_lower_bound, col  = 'green')
#PROBLEM 4

#################
M <- 1000 #size of price mesh
s_max <- 1000
s_min <- -1000
hs <- (2*s_max)/M
T <- (length(data)-44)/252 #time period
N <- 500 #size of time mesh
ht <- T/N
#create price series
s <- seq(s_min,s_max,hs)
#create time series
t <- seq(0,T,ht)
a <- (1-ht*sigma**2/(hs**2))
l <- (-kappa*(theta-s)*(ht/(2*hs))+(ht*sigma**2/(2*hs**2)))
u <- (kappa*(theta-s)*(ht/(2*hs))+(ht*sigma**2/(2*hs**2)))
A <- diag(a, nrow=M-1,M-1)
for (i in 1:M-2){
  A[i,i+1] <- u[i+1];
  A[i+1,i] <- l[i+2]
}

fun_y_preds <- function(time_input){
  return(as.numeric(S_ave$coefficients[1]) + as.numeric(S_ave$coefficients[2]) * (time_input+44/252))
}


find_H <- function(N, s, c){
  
  mat <- matrix(data = 0, nrow = M - 1, ncol = N)
  #H_T <- s + fun_y_preds(T) - c
  H_T <- s-c
  H_T <- H_T[2:(length(H_T) - 1)]
  HN <- H_T
  for (j in 1:(N)){
   # H_T <- s + fun_y_preds((T - (j) * ht)) - c
    H_T <- s - c
    b_end <- u[M] * (H_T[M + 1])
    b_start <- l[1] * H_T[1]
    H_T <- H_T[2:(length(H_T) - 1)]
    b <- matrix(data = c(b_start, rep(0, M - 3), b_end), nrow = M - 1, ncol = 1)
    
    HN <- A %*% HN + b
    mat_new <- matrix(data = 0, nrow = length(H_T), ncol = 1)
    for (w in 1:(length(H_T))){
      mat[w,N-j+1] <- max(HN[w, 1], H_T[w])
      mat_new[w, 1] <- max(HN[w, 1], H_T[w])
    }
    HN <- mat_new
  }
  final_mat <- matrix(data = 0, nrow = length(500), ncol = 1)
  for (a in 1:500){
    final_mat[a]<-(which(H_T>=mat[,a])[1])
  }
  return(s[final_mat]+S_ave$coefficients[1]+S_ave$coefficients[2]*t[1:500])
}

value <- find_H(N, s, c = 100)
value


find_G <- function(N, s, c){
  
  mat <- matrix(data = 0, nrow = M - 1, ncol = N)
  #H_T <- s + fun_y_preds(T) - c
  H_T <- s-c
  H_T <- H_T[2:(length(H_T) - 1)]
  HN <- H_T
  for (j in 1:(N)){
    # H_T <- s + fun_y_preds((T - (j) * ht)) - c
    H_T <- s - c
    b_end <- u[M] * (H_T[M + 1])
    b_start <- l[1] * H_T[1]
    H_T <- H_T[2:(length(H_T) - 1)]
    b <- matrix(data = c(b_start, rep(0, M - 3), b_end), nrow = M - 1, ncol = 1)
    
    HN <- A %*% HN + b
    mat_new <- matrix(data = 0, nrow = length(H_T), ncol = 1)
    for (w in 1:(length(H_T))){
      mat[w,N-j+1] <- max(HN[w, 1], H_T[w])
      mat_new[w, 1] <- max(HN[w, 1], H_T[w])
    }
    HN <- mat_new
  }
  H <- mat[,500]
  mat_G <- matrix(data = 0, nrow = M - 1, ncol = N)
  #H_T_G <- H - s - fun_y_preds(T) - c
  H_T_G <- H - s - c
  H_T_G <- H_T_G[2:(length(H_T_G) - 1)]
  HN_G <- H_T_G
  for (j in 1:(N)){
    # H_T_G <- H - s - fun_y_preds((T - (j) * ht)) - c
    H_T_G <- H - s - c
    b_end_G <- u[M] * (H_T_G[M + 1])
    b_start_G <- l[1] * H_T_G[1]
    H_T_G <- H_T_G[2:(length(H_T_G) - 1)]
    b_G <- matrix(data = c(b_start_G, rep(0, M - 3), b_end_G), nrow = M - 1, ncol = 1)
    
    HN_G <- A %*% HN_G + b_G
    mat_new_G <- matrix(data = 0, nrow = length(H_T_G), ncol = 1)
    for (w in 1:(length(H_T_G))){
      mat_G[w,N-j+1] <- max(HN_G[w, 1], H_T_G[w])
      mat_new_G[w, 1] <- max(HN_G[w, 1], H_T_G[w])
    }
    HN_G <- mat_new_G
  }
  final_mat_G <- matrix(data = 0, nrow = length(500), ncol = 1)
  for (a in 1:500){
    final_mat_G[a]<-(which(H_T_G<=mat_G[,a])[1])
  }
  return(s[final_mat_G]+S_ave$coefficients[1]+S_ave$coefficients[2]*t[1:500])
}

value2 <- find_G(N, s, c=100)
value2

moving_avg <- rollapply(data[45:67],3,mean)
n <- as.integer(500/21)
up_bound <- matrix(data=0,nrow=21,ncol=1)
low_bound <- matrix(data=0,nrow=21,ncol=1)
for (i in 1:21){
  up_bound[i] <- value[n*i]
  low_bound[i] <- value2[n*i]
}
df_OB <- data.frame("Rolling_Ave"=moving_avg,"Upper_bound"=up_bound,"Lower_bound"=low_bound)
SP500_price <- data[45:65]
time_days <- 1:21
ggdat <- data.frame(SP500_price, time_days)
plot(time_days,SP500_price,type='l',ylim=c(1000,4000))
lines(df_OB$Upper_bound,col='red')
lines(df_OB$Lower_bound, col='red')
