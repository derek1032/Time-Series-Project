options(stringsAsFactors = FALSE)
# Load data of spot and 1-month forward exchange rates
data = read.csv("proj15_spot_forward_exchange_rate.csv")
# Data setup
colnames(data) = c('date','AUD','AUD_f','JPY','JPY_f','GBP','GBP_f')
data$date = as.Date(data$date,"%m/%d/%Y")
# Convert the quoted form to US dollar/currency
data[,2:5] = 1/data[,2:5]
# Take logarithm
data[,2:7] = log(data[,2:7])

# Extract spot and forward exchange rate
spot = data.frame(date=data$date,AUD=data$AUD,JPY=data$JPY,GBP=data$GBP)
forward = data.frame(date=data$date,AUD_f=data$AUD_f,JPY_f=data$JPY_f,GBP_f=data$GBP_f)

# Plot log spot exchange rate
plot(spot[,1:2],typ='l',ylab='AUD',main='Exchange Rate')
plot(spot[,c(1,3)],typ='l',ylab='JPY',main='Exchange Rate')
plot(spot[,c(1,4)],typ='l',ylab='GBP',main='Exchange Rate')

par(mfrow=c(1,2))
# ACF and PACF
acf(spot$AUD);pacf(spot$AUD)
acf(spot$JPY);pacf(spot$JPY)
acf(spot$GBP);pacf(spot$GBP)
par(mfrow=c(1,1))
# ADF test
library(fUnitRoots)
adfTest(spot$AUD,type=c("nc"));adfTest(spot$JPY,type=c("nc"));adfTest(spot$GBP,type=c("nc"))

# Extract data before 2013-07-01
AUD_diff = diff(spot[spot$date<="2013-07-01",2])
JPY_diff = diff(spot[spot$date<="2013-07-01",3])
GBP_diff = diff(spot[spot$date<="2013-07-01",4])

date = spot[spot$date<="2013-07-01",1]

# Plot log return of exchange rate
plot(date[-1],AUD_diff,typ='l',ylab='AUD',main='Difference of Log Exchange Rate (AUD)')
plot(JPY_diff,typ='l',ylab='JPY',main='Difference of Log Exchange Rate (JPY)')
plot(GBP_diff,typ='l',ylab='GBP',main='Difference of Log Exchange Rate (GBP)')

# ACF and PACF of log return of exchange rate
par(mfrow=c(1,2))
acf(AUD_diff,main='Log Return of Exchange Rate(AUD)');pacf(AUD_diff,main='Log Return of Exchange Rate(AUD)')
acf(JPY_diff,main='Log Return of Exchange Rate(JPY)');pacf(JPY_diff,main='Log Return of Exchange Rate(JPY)')
acf(GBP_diff,main='Log Return of Exchange Rate(GBP)');pacf(GBP_diff,main='Log Return of Exchange Rate(GBP)')
par(mfrow=c(1,1))
# ADF test after differencing ##(No unit root)
adfTest(AUD_diff);adfTest(JPY_diff);adfTest(GBP_diff) 

hist(AUD_diff,breaks=20,ylab='Frequency',xlab='',main='Difference of Log Exchange Rate') 


install.packages('TSA')
library(TSA)
# Identify a good ARMA order using EACF
eacf(AUD_diff, ar.max = 8, ma.max = 8)
eacf(JPY_diff, ar.max = 8, ma.max = 8)
eacf(GBP_diff, ar.max = 8, ma.max = 8)

## ARCH effect analysis
t.test(AUD_diff)
t.test(JPY_diff)
t.test(GBP_diff)

# ACF of diff and squared diff of log exchange rate
par(mfrow=c(1,2))
acf(AUD_diff,main='Log Return of Exchange Rate(AUD)',ylim=c(-0.2,0.4)) 
acf(AUD_diff^2,main='Squared Log Return of Exchange Rate(AUD)',ylim=c(-0.2,0.4)) 
acf(JPY_diff,main='Log Return of Exchange Rate(JPY)',ylim=c(-0.2,0.4)) 
acf(JPY_diff^2,main='Squared Log Return of Exchange Rate(JPY)',ylim=c(-0.2,0.4)) ## not too much dependence
acf(GBP_diff,main='Difference of Exchange Rate(GBP)',ylim=c(-0.2,0.4)) 
acf(GBP_diff^2,main='Squared Log Return of Exchange Rate(GBP)',ylim=c(-0.2,0.4))  
par(mfrow=c(1,1))

# Ljung-Box test
Box.test(AUD_diff,lag=12,type=("Ljung-Box"))
Box.test(AUD_diff^2,lag=12,type=("Ljung-Box"))
Box.test(JPY_diff,lag=12,type=("Ljung-Box"))
Box.test(JPY_diff^2,lag=12,type=("Ljung-Box"))
Box.test(GBP_diff,lag=12,type=("Ljung-Box"))
Box.test(GBP_diff^2,lag=12,type=("Ljung-Box"))

library(fGarch)
# Fit data into Garch(1,1), mean equation is a constant 
m1 = garchFit(~garch(1,1),data=AUD_diff,trace=F)
m2 = garchFit(~garch(1,1),data=JPY_diff,trace=F)
m3 = garchFit(~garch(1,1),data=GBP_diff,trace=F)
summary(m1);summary(m2);summary(m3)

library(rugarch)
# IGARCH for AUD
spec = ugarchspec(variance.model=list(model="iGARCH",garchOrder=c(1,1)),
                         mean.model=list(armaOrder=c(0,0))) # mean equation=constant
m1_I = ugarchfit(spec=spec,data=AUD_diff)
m1_I
# APARCH
m1_ap = garchFit(~1+aparch(1,1), data=AUD_diff, trace=F)

# EGARCH for AUD
egarch11.spec = ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),
                           mean.model=list(armaOrder=c(0,0)))# mean equation=constant
m1_E = ugarchfit(egarch11.spec,data= AUD_diff)
#m2_E = ugarchfit(egarch11.spec,data= JPY_diff)
#m3_E = ugarchfit(egarch11.spec,data= GBP_diff)
summary(m1_E)

m1_ged = garchFit(~garch(1,1),data=AUD_diff,trace=F,cond.dist=c("ged"))
m2_ged = garchFit(~garch(1,1),data=JPY_diff,trace=F,cond.dist=c("ged"))
m3_ged = garchFit(~garch(1,1),data=GBP_diff,trace=F,cond.dist=c("ged"))
summary(m1_ged);summary(m2_ged);summary(m3_ged)
# In the case of AUD/US, alpha1 and beta1 are significant at the level 0.01
# We think GARCH(1,1) model can be used to predcit log return of exchange rate


m1_res = m1@residuals
m1_res_std = m1@residuals/volatility(m1)
m2_res = m2@residuals
m2_res_std = m2@residuals/volatility(m2)
m3_res = m3@residuals
m3_res_std = m3@residuals/volatility(m3)

# Observe volatility and log return of exchange rate 
plot(AUD_diff,typ='l',ylab='USD/AUD')
lines(volatility(m1),col='red')
legend(160,0.1,c('log return of USD/AUD','volatility'),col=c(1,2),lwd=c(2,2))

# ACFs 
par(mfrow=c(2,2)) 
acf(m1_res,main='Residual (AUD)') 
acf(m1_res^2,main='Residual Squared (AUD)') 
acf(m1_res_std,main='GARCH(1,1) Std Residual (AUD)') 
acf(m1_res_std^2,main='GARCH(1,1) Std Residual Squared (AUD)') 

acf(m2_res,main='Residual (JPY)') 
acf(m2_res^2,main='Residual Squared (JPY)') 
acf(m2_res_std,main='GARCH(1,1) Std Residual (JPY)') 
acf(m2_res_std^2,main='GARCH(1,1) Std Residual Squared (JPY)') 

acf(m3_res,main='Residual (GBP)') 
acf(m3_res^2,main='Residual Squared (GBP)') 
acf(m3_res_std,main='GARCH(1,1) Std Residual (GBP)') 
acf(m3_res_std^2,main='GARCH(1,1) Std Residual Squared (GBP)') 

# PACFs 
pacf(m1_res,main='Residual (AUD)',ylim=c(-0.1,0.5)) 
pacf(m1_res^2,main='Residual Squared (AUD)',ylim=c(-0.1,0.5)) 
pacf(m1_res_std,main='GARCH(1,1) Std Residual (AUD)',ylim=c(-0.1,0.5)) 
pacf(m1_res_std^2,main='GARCH(1,1) Std Residual Squared (AUD)',ylim=c(-0.1,0.5)) 

pacf(m2_res,main='Residual (JPY)',ylim=c(-0.15,0.5)) 
pacf(m2_res^2,main='Residual Squared (JPY)',ylim=c(-0.15,0.5)) 
pacf(m2_res_std,main='GARCH(1,1) Std Residual (JPY)',ylim=c(-0.15,0.5)) 
pacf(m2_res_std^2,main='GARCH(1,1) Std Residual Squared (JPY)',ylim=c(-0.15,0.5)) 

pacf(m3_res,main='Residual (GBP)',ylim=c(-0.1,0.5)) 
pacf(m3_res^2,main='Residual Squared (GBP)',ylim=c(-0.1,0.5)) 
pacf(m3_res_std,main='GARCH(1,1) Std Residual (GBP)',ylim=c(-0.1,0.5)) 
pacf(m3_res_std^2,main='GARCH(1,1) Std Residual Squared (GBP)',ylim=c(-0.1,0.5)) 

par(mfrow=c(1,1)) 
plot(m1_res_std,typ='l',ylab='',main='Standardized Residuals (AUD)')
plot(m2_res_std,typ='l',ylab='',main='Standardized Residuals (JPY)')
plot(m3_res_std,typ='l',ylab='',main='Standardized Residuals (GBP)')

n = nrow(spot) # number of full data
m = nrow(spot[spot$date<="2013-07-01",]) 
realized_return = as.data.frame(matrix(rep(0,(n-m)*3),nrow=n-m))
colnames(realized_return)=c('AUD','JPY','GBP')
data_temp=c()
for(j in 1:3)
{
  for (i in 1:(n-m))
  {
    data_temp = diff(spot[1:m-1+i,j+1])
    mdl = garchFit(~garch(1,1),data=data_temp,trace=F)
    pred = predict(mdl,n.ahead = 1)[,3] # predict the volatility of the next day
    mdl_res_std = mdl@residuals/volatility(mdl)
    epsilon = mdl_res_std[length(mdl_res_std)]
    if (mdl@fit$coef[1]+pred * epsilon > (forward[m-1+i,j+1] - spot[m-1+i,j+1]))
      realized_return[i,j] = spot[m+i,j+1] - forward[m-1+i,j+1]
    else
      realized_return[i,j] = forward[m-1+i,j+1] - spot[m+i,j+1]
  }
}
rm(i,j,mdl,pred,epsilon,data_temp)
colMeans(realized_return, na.rm = FALSE, dims = 1)


# Using Generalized Error as underlying distributions of epsilon_t
realized_return_ged = as.data.frame(matrix(rep(0,(n-m)*3),nrow=n-m))
colnames(realized_return_ged)=c('AUD','JPY','GBP')
data_temp=c()

for(j in 1:3)
{
  for (i in 1:(n-m))
  {
    data_temp = diff(spot[1:m-1+i,j+1])
    mdl = garchFit(~garch(1,1),data=data_temp,trace=F,cond.dist=c("ged"))
    pred = predict(mdl,n.ahead = 1)[,3] # predict the volatility of the next day
    mdl_res_std = mdl@residuals/volatility(mdl)
    epsilon = sample(mdl_res_std, 1, replace = TRUE)
    if (mdl@fit$coef[1]+pred * epsilon > (forward[m-1+i,j+1] - spot[m-1+i,j+1]))
      realized_return_ged[i,j] = spot[m+i,j+1] - forward[m-1+i,j+1]
    else
      realized_return_ged[i,j] = forward[m-1+i,j+1] - spot[m+i,j+1]
  }
} 
rm(i,j,mdl,pred,epsilon,data_temp)

# OLS methods 
realized_return_ols = as.data.frame(matrix(rep(0,(n-m)*3),nrow=n-m))
colnames(realized_return_ols)=c('AUD','JPY','GBP')
alpha_temp = c()
beta_temp = c()
for(j in 1:3)
{
  for (i in 1:(n-m))
  {
    y = diff(spot[1:(m-1+i),j+1])
    X = forward[1:(m-2+i),j+1]-spot[1:(m-2+i),j+1] 
    mdl = lm(y~X)
    if (mdl$coefficients[1]+ mdl$coefficients[2]*((forward[m-1+i,j+1] - spot[m-1+i,j+1]))> (forward[m-1+i,j+1] - spot[m-1+i,j+1]))
      realized_return_ols[i,j] = spot[m+i,j+1] - forward[m-1+i,j+1]
    else
      realized_return_ols[i,j] = forward[m-1+i,j+1] - spot[m+i,j+1]
  }
}
rm(i,j,y,X,b)
colMeans(realized_return_ols, na.rm = FALSE, dims = 1)

## IGARCH
install.packages("rugarch")
library(rugarch)
realized_return_AUD = c()
for (i in 1:(n-m))
{
  data_temp = diff(spot[1:m-1+i,2])
  spec = ugarchspec(variance.model=list(model="iGARCH",garchOrder=c(1,1)),
                    mean.model=list(armaOrder=c(0,0))) # mean equation=constant
  mdl = ugarchfit(spec=spec,data=data_temp)
  # standardized residuals
  residual_std = mdl@fit$residuals/mdl@fit$sigma
  # predict the volatility of the next day
  pred1 = as.numeric(sqrt(mdl@fit$coef[2]+mdl@fit$coef[3]*mdl@fit$residuals[length(mdl@fit$residuals)]^2
                    +mdl@fit$coef[4]*mdl@fit$sigma[length(mdl@fit$sigma)]^2))
  # select the last standardized residual as the next one 
  epsilon = residual_std[length(residual_std)]
  # sample(residual_std,size=1,replace=TRUE)
  if (mdl@fit$coef[1]+pred1 * epsilon > (forward[m-1+i,2] - spot[m-1+i,2]))
    realized_return_AUD[i] = spot[m+i,2] - forward[m-1+i,2]
  else
    realized_return_AUD[i] = forward[m-1+i,2] - spot[m+i,2]
  
}
rm(spec,mdl,residual_std,pred1,epsilon)
realized_return_AUD = data.frame(rtn = realized_return_AUD)
mean(realized_return_AUD$rtn)

# Relevant measures
TB = read.csv("TB30.csv")
SP500 = read.csv("sp500.csv")
TB30 = cbind(TB[,2],TB[,2],TB[,2])
excess_return_3 =  realized_return-TB30
excess_return_mkt = SP500[,2] - TB[,2]
excess_return_aud = realized_return_AUD$rtn-TB$t30ret
colMeans(excess_return_3, na.rm = FALSE, dims = 1)

#t.test(excess_return_3[,1],mu=0)
#t.test(excess_return_3[,2],mu=0)
#t.test(excess_return_3[,3],mu=0)

sharpe_ratio_aud = mean(excess_return_3$AUD) /sd(excess_return_3$AUD)
sharpe_ratio_jpy = mean(excess_return_3$JPY) /sd(excess_return_3$JPY)
sharpe_ratio_gbp = mean(excess_return_3$GBP) /sd(excess_return_3$GBP)
sharpe_ratio_mkt = mean(excess_return_mkt) /sd(excess_return_mkt)
sharpe_ratio_aud_igarch = mean(excess_return_aud) /sd(excess_return_aud)
cat('Sharpe Ratios:',sharpe_ratio_aud,sharpe_ratio_jpy,sharpe_ratio_gbp)
# Count the winning months
win = c(sum(realized_return$AUD>0),sum(realized_return$JPY>0),sum(realized_return$GBP>0))
sum(realized_return_AUD$rtn>0)
# Count the lsoing months
lose = c(sum(realized_return$AUD<0),sum(realized_return$JPY<0),sum(realized_return$GBP<0))
sum(realized_return_AUD$rtn<0)
cat('USD/AUD',sharpe_ratio_aud,'USD/JPY',sharpe_ratio_jpy,'USD/GBP',sharpe_ratio_gbp,'MARKET',sharpe_ratio_mkt)
# Plot realized returns of GARCH(1,1) and OLS
par(mfrow=c(3,1))
plot(realized_return[,1],type='l',xlab='month',
     ylab='Return of AUD',col='red',main='Return comparison of GARCH(1,1) and OLS')
lines(realized_return_ols[,1],type='l',col='blue')
legend(10,0.045,c('GARCH','OLS'),col=c(2,4),lwd=c(2,2),bty='n',cex=1)
plot(realized_return[,2],type='l',xlab='month',
     ylab='Return of JPY',col='red',main='Return comparison of GARCH(1,1) and OLS')
lines(realized_return_ols[,2],type='l',col='blue')
legend(10,0.045,c('GARCH','OLS'),col=c(2,4),lwd=c(2,2),bty='n',cex=1)
plot(realized_return[,3],type='l',xlab='month',
     ylab='Return of GBP',col='red',main='Return comparison of GARCH(1,1) and OLS')
lines(realized_return_ols[,3],type='l',col='blue')
legend(10,0.045,c('GARCH','OLS'),col=c(2,4),lwd=c(2,2),bty='n',cex=1)
par(mfrow=c(1,1))

## daily data
data_daily = read.csv("proj15_daily_exchange_rate.csv")

# Data setup
data_daily = data_daily[!is.na(data_daily[,2]),1:4]
colnames(data_daily) = c('date','AUD','JPY','GBP')
data_daily$date = as.Date(data_daily$date,"%d-%B-%y")

data_daily[,2:4] = 1/data_daily[,2:4]

data_daily[,2:4]=log(data_daily[,2:4])

# Extract data before 2013-07-01
AUD_diff = diff(data_daily[data_daily$date<="2013-07-01",2])
JPY_diff = diff(data_daily[data_daily$date<="2013-07-01",3])
GBP_diff = diff(data_daily[data_daily$date<="2013-07-01",4])

plot(AUD_diff,typ='l',ylab='',main='AUD/US')
plot(JPY_diff,typ='l',ylab='',main='JPY/US')
plot(GBP_diff,typ='l',ylab='',main='GBP/US')

eacf(AUD_diff)

library(rugarch)
egarch11.spec = ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),
                           mean.model=list(armaOrder=c(0,0)))# mean equation=constant
m1_E = ugarchfit(egarch11.spec,data= AUD_diff)
m2_E = ugarchfit(egarch11.spec,data= JPY_diff)
m3_E = ugarchfit(egarch11.spec,data= GBP_diff)
library(fGarch)
m1_ap = garchFit(~1+aparch(1,1), data=AUD_diff, trace=F)
m2_ap = garchFit(~1+aparch(1,1), data=JPY_diff, trace=F)
m3_ap = garchFit(~1+aparch(1,1), data=GBP_diff, trace=F)
residual_m1_E = m1_E@fit$residuals
std_residual_m1_E=m1_E@fit$residuals/m1_E@fit$sigma
plot(std_residual_m1_E,typ='l',ylab='Standardized Residuals')
par(mfrow=c(2,2))
acf(residual_m1_E,main='Residual');acf(residual_m1_E^2,main='Squared Residual')
acf(std_residual_m1_E,main='EGARCH(1,1) Standardized Residual');acf(std_residual_m1_E^2,main='EGARCH(1,1) Squared Standardized Residual')
Box.test(m1_E@fit$residuals,lag=12,type=("Ljung-Box"))
Box.test(std_residual_m1_E^2,lag=12,type=("Ljung-Box"))


