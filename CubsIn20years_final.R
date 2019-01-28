
###############################################
# Import Libraries
###############################################
library(ggplot2)
library(corrplot)
library(tseries)
library(fBasics)
library(zoo)
library(TTR)
library(forecast)
library(quantmod)
library(lmtest)
library(TSA)
library(plyr)
library(fUnitRoots)
library(seasonal)
library(astsa)
library(fGarch)

path <- "/Users/Queena/Documents/CDMacademic/CSC425/Cubs"
setwd(path)

###############################################
# Data Import and Format Correction
###############################################
Cubs = read.csv("Cubs Game Logs.csv",header=T, sep=',',stringsAsFactors=FALSE)

# Date Format Correction
pattern_right = "^[1-9]"
t_right = grepl(pattern_right, Cubs$Date)
Cubs$Date[which(t_right == FALSE)]

pattern_org = "^[A-Z]" #pattern not in common format
t_org = grepl(pattern_org, Cubs$Date)
Cubs$Date[which(t_right==FALSE&t_org==TRUE)]
Cubs$Date[which(t_right==FALSE&t_org==FALSE)]

for (i in 1:nrow(Cubs)){
  date = as.character(Cubs$Date[i])
  data = sub('susp',"",date)
  if (grepl("^[A-Z]", date)==TRUE){
    month = strsplit(date,' ')[[1]][1]
    day = strsplit(date,' ')[[1]][2]
    day = gsub("[^0-9]", "", day)
    Cubs$Date[i] = paste(c(day,month),collapse='-')
  }
}

Cubs$DATE = as.Date(with(Cubs, paste(Year, Date, sep='-')),
                    "%Y-%d-%b")

any(is.na(Cubs$DATE)) 

###############################################
# Preprocess
###############################################

# Issue: unevenly spaced time series -- to differenct degree
# an unevenly (or unequally or irregularly) spaced time series is 
# a sequence of observation time and value pairs (tn, Xn) 
# with strictly increasing observation times. 
# A common approach to analyzing unevenly spaced time series is 
# to transform the data into equally spaced observations using some form of interpolation
# good: not highly irregular, th new created biases not that significant

# sub 1: In one season, it's not evenly distributed.
# In the last 10 years, there are 43 days with two games one day
# For example, there are 4 days in 2018 with two games 
dailygame = as.data.frame(count(Cubs,'DATE'))
subset_morethan1 = subset(dailygame,dailygame$freq!=1)
subset_morethan1

# sub 2: different total games count in some years
yearlygame = as.data.frame(count(Cubs, 'Year'))
yearlygame 
years = seq(1999,2018)

# solution to sub 1: 
correction_DATE_Year = rep(yearlygame$Year, yearlygame$freq)
correction_DATE_Day = as.Date(NA)

for (year in years){
  day = 1
  month = 4
  start_day = paste(c(year,month,day),collapse='-')
  daysseq = seq(as.Date(start_day), by = 1, 
               length=yearlygame$freq[yearlygame$Year==year])
  correction_DATE_Day = c(correction_DATE_Day, daysseq)
}

correction_DATE = correction_DATE_Day[-1]

Cubs = cbind(Cubs,correction_DATE)
head(Cubs)

# sub 3: How about time between seasons -- almost 6 month, same as the in-season length.
# solution: try interventional analysis

# create time series -- disconnected, and connected
cubs_dc = zoo(x=Cubs$SO, correction_DATE)
cubs_c = ts(Cubs$SO,
            frequency = 162,
            start = c(1999))

plot(cubs_dc)
plot(cubs_c)
print(cubs_c, calendar =T)

###############################################
# Basic Exploration 
###############################################
# basic statistics
summary(cubs_dc)
basicStats(cubs_dc)

# Normality exploration
hist(cubs_dc, xlab="Strike Out Count", prob=TRUE, main="Histogram of Strike Out Count (1999-2018)")
curve(dnorm(x, mean=mean(cubs_dc), sd=sd(cubs_dc)), add=TRUE, col = "red")
## Normal quantile plot
qqnorm(cubs_dc, main="Strike Out Count Normal Q-Q Plot (1999-2018)")
qqline(cubs_dc) 
## the Jarque Bera test (5% significance level)
normalTest(cubs_dc, method="jb")

# plot
plot(cubs_dc, main="Strike Out Count (1999-2018)")

ggplot(Cubs, aes(x=Gtm, y=SO, colour=as.factor(Year))) + geom_line() + geom_smooth() + geom_smooth(col='blue') + 
  ggtitle("Strike Out Count (1999-2018) with Loess Smoothing") +xlab('Time') +
  geom_smooth(col='blue')

# White Noise test
Box.test(cubs_dc, lag = 8)

###############################################
# Basic Exploration -- Stationarity test
###############################################
# Stationarity test
# ACF
acf(cubs_c,
    type = c("correlation"), sub="1999-2018 (lagmax = 100)",
    lag.max = 100)
# PACF
pacf(cubs_c, sub="1999-2018")
# EACF
eacf(cubs_c)
t = eacf(cubs_c, ar.max = 150, ma.max = 150)
t_char = as.data.frame(t[['symbol']])

# Create each year ACF
# It's more like white noise if just look at each year's series
par(mfrow=c(4,4))
for (year in seq(1999,2006)){
  year = as.integer(year)
  acfs = acf(Cubs$SO[Cubs$Year==year],
             type = c("correlation"),lag.max = 50,
             plot = FALSE)
  plot(acfs, main = paste(c(year)))
  pacfs = pacf(Cubs$SO[Cubs$Year==year],
               lag.max = 50, plot = FALSE)
  plot(pacfs, main = paste(c(year)))
}
for (year in seq(2007,2014)){
  year = as.integer(year)
  acfs = acf(Cubs$SO[Cubs$Year==year],
             type = c("correlation"),lag.max = 50,
             plot = FALSE)
  plot(acfs, main = paste(c(year)))
  pacfs = pacf(Cubs$SO[Cubs$Year==year],
               lag.max = 50, plot = FALSE)
  plot(pacfs, main = paste(c(year)))
}
for (year in seq(2015,2018)){
  year = as.integer(year)
  acfs = acf(Cubs$SO[Cubs$Year==year],
             type = c("correlation"),lag.max = 50,
             plot = FALSE)
  plot(acfs, main = paste(c(year)))
  pacfs = pacf(Cubs$SO[Cubs$Year==year],
               lag.max = 50, plot = FALSE)
  plot(pacfs, main = paste(c(year)))
}
par(mfrow=c(1,1))

# Decompose
# classic
cubs_c %>% decompose() %>%
  autoplot() + xlab("Year") +
  ggtitle("Cubs Strikeout Count 1999-2018")

cubs_c_decompose = decompose(cubs_c)
plot(cubs_c_decompose[["seasonal"]])

#########################################
## Try different differencing
#########################################
# get the difference needed to make it stationary
nsdiffs(cubs_c)

par(mfrow=c(2,2))

d1 = diff(cubs_c)
plot(d1, main='Difference by lag 1')
acf(d1, lag.max = 200)
pacf(d1)
Box.test(d1) # not white noise
adfTest(d1) # small pvalue 0.01

d9 = diff(cubs_c, 9)
plot(d9, main='Difference by lag 9')
acf(d9, lag.max = 200)
pacf(d9)
Box.test(d9) # white noise
adfTest(d9) # small pvalue 0.01

d10 = diff(cubs_c, 10)
plot(d10,main='Difference by lag 10')
acf(d10, lag.max = 200)
pacf(d10)
Box.test(d10) # white noise
adfTest(d10) # small pvalue 0.01


par(mfrow=c(2,2))
d37 = diff(cubs_c, 37)
plot(d37,main='Difference by lag 37')
acf(d37, lag.max = 200)
pacf(d37)
Box.test(d37) # not white noise 0.02065
adfTest(d37) # small pvalue 0.01

d162 = diff(cubs_c, 162)
plot(d162,main='Difference by lag 162')
acf(d162, lag.max = 200)
pacf(d162)
Box.test(d162) # not white noise
adfTest(d162) # small pvalue 0.01

## Based on Dr. McDonald's instruction
dd9 = diff(d9)
acf(dd9, lag.max = 100)
pacf(dd9)
Box.test(dd9, type='Ljung-Box') # not white noise
adfTest(dd9) # small pvalue 0.01

dd10 = diff(d10)
acf(dd10, lag.max = 100)
pacf(dd10)
Box.test(dd10, type='Ljung-Box') # not white noise
adfTest(dd10) # small pvalue 0.01

par(mfrow=c(3,1))
plot(cubs_c,main='Original')
plot(d9,main='9-Lag Difference')
plot(d10, main = '10-Lag Difference')

acf(cubs_c, lag.max = 200)
acf(d9, lag.max = 200)
acf(d10, lag.max = 200)

par(mfrow=c(1,1))

#########################################
## Naive Modeling -- ARIMA
#########################################
fit1_c = auto.arima(cubs_c, seasonal=TRUE, trace = TRUE,
                 allowdrift = TRUE)
cubs_c10 = ts(Cubs$SO,
              start = c(1999),
              frequency = 10)
cubs_c9 = ts(Cubs$SO,
              start = c(1999),
              frequency = 9)
fit1_c10 = auto.arima(cubs_c10, seasonal=TRUE, trace = TRUE,
                    allowdrift = TRUE)
fit1_c9 = auto.arima(cubs_c9, seasonal=TRUE, trace = TRUE,
                      allowdrift = TRUE)
summary(fit1_c10)
summary(fit1_c)
AIC(fit1_c9)

## Get inspiration based on automatic selection
# To select the best season leading to highest AIC in auto.arima
# To avoid double run this process, all the code below are marked as comment
cubs_c1 = ts(Cubs$SO,
             start = c(1999),
             frequency = 1)
fit_original = auto.arima(cubs_c1, seasonal=TRUE, trace = TRUE,
                          allowdrift = TRUE)
cub_cs_original = cubs_c1

#for (i in seq(2,150)){
#  cub_cs = ts(Cubs$SO,
#              start = c(1999),
#              frequency = i)
#  fit_new = auto.arima(cub_cs, seasonal=TRUE, trace = TRUE,
#                       allowdrift = TRUE)
#  if (AIC(fit_new)<AIC(fit_original)){
#    fit_original = fit_new
#    cub_cs_original = cub_cs
#  }
#}

# The result of the above code is fit_original
summary(fit_original)
summary(cub_cs_original)

## Based on the result, the best season selection is 37.
# There are two models
# ARIMA(1,1,2) -- based on original seasonality (162)
# ARIMA(1,1,2)(0,0,1)[37] -- based on the selection

#########################################
## ARIMA(1,1,2): fit1_c
#########################################
summary(fit1_c)

# Significant Test
coeftest(fit1_c)

# Residual Test: it's not normally distributed white noise
### Normality -- skewed and with kurtosis issue
par(mfrow=c(2,2))
plot(residuals(fit1_c), main = 'Residual of ARIMA(1,1,2)')
hist(residuals(fit1_c), main = 'Residual of ARIMA(1,1,2)')
qqnorm(residuals(fit1_c), main = 'Residual of ARIMA(1,1,2)') 
qqline(residuals(fit1_c), col = 2)
jarque.bera.test(fit1_c$residuals)
skewness(fit1_c$residuals)
kurtosis(fit1_c$residuals)
### White noise test -- it's white noise
Box.test(fit1_c$residuals,lag=8,type='Ljung')
acf(fit1_c$residuals, main = 'Residual of ARIMA(1,1,2)')

## backtest
source("backtest.R")
backtest(fit1_c, cubs_c, trunc(length(cubs_c)*.8), 1)

#########################################
## ARIMA(1,1,2)(0,0,1)[37]: fit_original
#########################################
# summary(fit_original)

fit2 = arima(cubs_c, order = c(1,1,2), 
             seasonal=list(order= c(0,0,1), period=37))
summary(fit2)

# Significant Test
coeftest(fit2)

# Residual Test: it's not normally distributed white noise
### Normality -- skewed and with kurtosis issue
plot(residuals(fit2), main = 'Residual of ARIMA(1,1,2)(0,0,1)[37]')
hist(residuals(fit2), main = 'Residual of ARIMA(1,1,2)(0,0,1)[37]')
qqnorm(residuals(fit2), main = 'Residual of ARIMA(1,1,2)(0,0,1)[37]') 
qqline(residuals(fit2), col = 2)
jarque.bera.test(fit2$residuals)
skewness(fit2$residuals)
kurtosis(fit2$residuals)
### White noise test -- it's white noise
Box.test(fit2$residuals,lag=8,type='Ljung')
acf(fit2$residuals, main = 'Residual of ARIMA(1,1,2)(0,0,1)[37]')

## backtest
# backtest(fit_original, cubs_c, trunc(length(cubs_c)*.8), 1)


#########################################
## Other manual seasonality modeling
#########################################
# Based on the plot, there is kinda of 5-year seasonality
# Best model is still ARIMA(1,1,2)

cubs_c5y = ts(Cubs$SO,
             start = c(1999),
             frequency = 162*5)
fit_c5y = auto.arima(cubs_c5y, seasonal=TRUE, trace = TRUE,
                     allowdrift = TRUE)

# In original traditional modeling approach, 
# there are seasonality limit -- for example, maximum lag support is 350
fit_c5y_2 = arima(cubs_c, order = c(1,1,2), 
                  seasonal=list(order= c(0,0,1), period=162*5))
summary(fit_c5y_2)
 
############################################
## GARCH based on residuals -- doesn't work
############################################
# Inspiration
par(mfrow = c(2,1))
plot(diff(cubs_c),main='Plot of ordinary differencing')
plot(cubs_c,main='Plot of original strickout counts')
par(mfrow = c(1,1))

residuals = fit2$residuals

gfit.ts1 = garch(residuals) 
gfit.ts1
coeftest(gfit.ts1)

red1 = gfit.ts1$residuals[-1]
plot(red1)

skewness(red1)
kurtosis(red1)   # Much less
jarque.bera.test(red1)   # Cannot reject normality


gfit.ts2 = garchFit(~ garch(1, 1), data=residuals, trace=F)
gfit.ts2
red2 = gfit.ts2@residuals[-1]
plot(red2)

skewness(red2)
kurtosis(red2)   # Only a few less
jarque.bera.test(red2)   # reject normality

## Normality didn't change a lot. So this is not the right approach to 
## process the residuals

#####################################################
## GARCH -- in more original dataset, diff(cubs_c)
#####################################################
# How about apply the GARCH on data directly
# Inspiration
par(mfrow = c(2,1))
plot(fit2$residuals,main='Plot of ordinary differencing')
plot(cubs_c,main='Plot of original strickout counts')
par(mfrow = c(1,1))

cubs_diff = diff(cubs_c)

# absolute value
cubs_abdiff = abs(cubs_diff)
plot(cubs_abdiff)
acf(cubs_abdiff)
pacf(cubs_abdiff)

# squared value
cubs_squared = cubs_diff*cubs_diff
plot(cubs_squared)
acf(cubs_squared)
pacf(cubs_squared)

# Start with exploration GARCH(1,1)
gfit.ts3 = garchFit(~ garch(1, 1), data=cubs_diff, trace=F)
gfit.ts3

red3 = gfit.ts3@residuals[-1]
plot(red3)

skewness(red3)   # Much less 0.02445281 (compared with 0.6137563)
kurtosis(red3)   # Much less 0.765704 (compared with 1.283936)
jarque.bera.test(red3)   #  Reject normality

#plot(gfit.ts3)

# Looks the right to thing
# Start with exploration GARCH(2,2) to find the best parameter for GARCH
gfit.ts4 = garchFit(~ garch(2, 2), data=cubs_diff, trace=F)
summary(gfit.ts4)

# alpha 1 only 
gfit.ts5 = garchFit(~ garch(1, 0), data=cubs_diff, trace=F)
gfit.ts5
gfit.ts5@fit$ics[1]
backtest(gfit.ts5, cubs_diff, trunc(length(cubs_diff)*.8), 1)

red5 = gfit.ts5@residuals[-1]
plot(red5)

skewness(red5)   # Much less 0.02445281 (compared with 0.6137563)
kurtosis(red5)   # Much less 0.765704 (compared with 1.283936)
jarque.bera.test(red5)   #  Reject normality

# plot(gfit.ts5)


# How about have GARCH and ARIMA together
gfit.ts6 = garchFit(~ arma(1,1) + garch(1, 0), data=cubs_diff, trace=F)
gfit.ts6
gfit.ts7 = garchFit(~ arma(1,1) + garch(2, 2), data=cubs_diff, trace=F)
gfit.ts7
gfit.ts8 = garchFit(~ arma(1,1) + garch(1, 2), data=cubs_diff, trace=F)
gfit.ts8
gfit.ts9 = garchFit(~ arma(0,1) + garch(1, 0), data=cubs_diff, trace=F)
gfit.ts9

red6 = gfit.ts6@residuals[-1]
plot(red5)

skewness(red6)   # Much less 0.02445281 (compared with 0.6137563)
kurtosis(red6)   # Much less 0.765704 (compared with 1.283936)
jarque.bera.test(red6)   #  Reject normality

#plot(gfit.ts5)

#########################################
## GARCH -- in cubs_c
#########################################
# How about apply the GARCH on data directly
# Inspiration
par(mfrow = c(2,1))
plot(fit2$residuals,main='Plot of ordinary differencing')
plot(cubs_c,main='Plot of original strickout counts')
par(mfrow = c(1,1))

gfit.ts10 = garchFit(~ garch(2, 2), data=cubs_c, trace=F)
summary(gfit.ts10)

gfit.ts11 = garchFit(~ garch(1, 1), data=cubs_c, trace=F)
summary(gfit.ts11)

gfit.ts12 = garchFit(~ arma(1,2) + garch(1, 1), data=cubs_c, trace=F)
summary(gfit.ts12)

gfit.ts13 = garchFit(~ arma(1,1) + garch(1, 1), data=cubs_c, trace=F)
summary(gfit.ts13)

red13 = gfit.ts13@residuals[-1]
plot(red13)

skewness(red13)   # higher 0.6207731 (compared with 0.6137563)
kurtosis(red13)   # higher 1.30887 (compared with 1.283936)
jarque.bera.test(red13)



fit_nodif = arima(cubs_c, order = c(1,0,1))
summary(fit_nodif)
coeftest(fit_nodif)
backtest(fit_nodif, cubs_c, trunc(length(cubs_c)*.8), 1)

fit_try = auto.arima(cubs_c, d=0)

#########################################
## Naive Modeling -- ets
#########################################
# Exponential smoothing model
t= ets(cubs_c1, model="ZZZ", damped=NULL, alpha=NULL, beta=NULL,
       gamma=NULL, phi=NULL, lambda=NULL, biasadj=FALSE,
       additive.only=FALSE, restrict=TRUE,
       allow.multiplicative.trend=FALSE)
print(t)
summary(t)
plot(residuals(t))

# 

cubs_c1 %>% mstl() %>%
  autoplot() + xlab("Week")
cubs_c %>%  stlf() %>%
  autoplot() + xlab("Week")













fit <- stl(cubs_c, t.window=163, s.window="periodic",
           robust=TRUE)
fit %>% seasadj() %>% naive() %>%
  autoplot() + ylab("New orders index") +
  ggtitle("Naive forecasts of seasonally adjusted data")
fit %>% forecast(method="naive") %>%
  autoplot() + ylab("New orders index")







#########################################
## Manual Programming
#########################################
cubs_c37 = ts(Cubs$SO,
             frequency = 37,
             start = c(1999))
fit_nodiff4 = arima(cubs_c, order = c(1,0,1),
                    seasonal = list(order=c(0,0,1), period=37))
summary(fit_nodiff4)
