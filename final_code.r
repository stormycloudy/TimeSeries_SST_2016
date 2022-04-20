# Codes used in these project

# Re-run with caution: all figures can be re-produced, but some lines need to be commented out or un-commented. 
# Properties of models may need manually typing to show.
# ------------------------------------------------------------------------------------------------------------------------

D <- read.csv("/N/u/xuecliu/Karst/Downloads/rpro/sst_nino34.csv", header=T) 
SST=D$SST
SOI=D$SOI
summary(D)

plot(D$SST,xlab="Time (month)",ylab="Sea Surface Temperature (deg C)",type="l")
plot(D$SOI,xlab="Time (month)",ylab="Southern Oscillation Indices",type="l")

#Spectral analysis

install.packages("FitAR")
library(FitAR)

ts1=spectrum(SST,log="no")
1/ts1$freq[which.max(ts1$spec)] #the 36th freq=0.833, period=12

ts2=spectrum(SOI,log="no")
1/ts2$freq[which.max(ts2$spec)] #period = 144

par(mfrow=c(1,3))
acf(SST,lag.max = 30)
acf(diff(SST),lag.max = 30)
#acf((diff(SST)))
pacf(diff(SST),lag.max = 30)
par(mfrow=c(1,3))
acf(SOI,lag.max = 30)
acf(diff(SOI),lag.max = 30)
pacf(diff(SOI),lag.max = 30)

#AR & MA
m1=Arima(diff(SST), order=c(1,0,0))
m1=Arima(diff(SST), order=c(0,0,1))
m1=Arima(diff(SST), order=c(0,0,18))
ar18=Arima(diff(SST), order=c(18,0,0))
m1 = Arima(diff(SST), order=c(26,0,0) )
m2 = Arima(diff(SST), order=c(0,0,26) )
ar36 = Arima(diff(SST), order=c(36,0,0) )
m1 = Arima(diff(SST), order=c(0,0,36) )
m1 = Arima(diff(SST), order=c(26,0,2) )
dev.off()
par(mfrow=c(1,2))
  acf(ar18$residuals,lag.max = 40)
  acf(ar36$residuals,lag.max = 40)

dev.off()
par(mfrow=c(2,3))
  #plot(m1$residuals)
  acf(m1$residuals,lag.max = 40)
  acf(m1$residuals,lag.max = 400)
  LBQPlot(m1$residuals,lag.max=400,StartLag=1)
  acf(m2$residuals,lag.max = 40)
  acf(m2$residuals,lag.max = 400)
  LBQPlot(m2$residuals,lag.max=400,StartLag=1)
  P1=spectrum(m1$residuals,log="no")
  1/P1$freq[which.max(P1$spec)]

dev.off()
lmfit = lm(SST~c(1:406))
plot(lmfit)
plot(lmfit$residuals)
AIC(lmfit)
BIC(lmfit)

#ARIMA
m3 = Arima(SST, order=c(26,1,0) )
m3 = Arima(SST, order=c(0,1,26) )
m3 = Arima(SST, order=c(0,2,26) )
m3 = Arima(SST, order=c(26,2,0) )


par(mfrow=c(2,2))
  acf(m3$residuals,lag.max = 400)
  LBQPlot(m3$residuals,lag.max=40,StartLag=1)
  LBQPlot(m3$residuals,lag.max=400,StartLag=1)
  P3=spectrum(m3$residuals,log="no")
  1/P3$freq[which.max(P1$spec)]
  plot(m3$residuals)

  par(mfrow=c(1,1))
  pacf(m3$residuals)  
  plot(m3$residuals)
  
Box.test(m3$residuals,lag=20)  

par(mfrow=c(1,1))
fit <- Arima(SST, order=c(26,0,26) )
plot(fit$x,col="black",type="l") #original
lines(fitted(fit),col="red") #fitted

dev.off()
#SARIMA
library(astsa)
ms1 = sarima(SST, 0,1,26, 0,1,1, 6, details=FALSE)
ms1 = sarima(SST, 26,1,0, 0,1,1, 12, details=FALSE)
fit1 <-Arima(SST, order=c(26,1,0), seasonal=list(order=c(0,1,1),period=12))
par(mfrow=c(2,2))
  plot(ms1$fit$residuals,type="l")
  acf(ms1$fit$residuals,lag.max = 400)
  P2=spectrum(ms1$fit$residuals,log="no")
    1/P2$freq[which.max(P1$spec)] # =6
  LBQPlot(ms1$fit$residuals,lag.max=400,StartLag=1)

ms2 = sarima(SST, 0,1,26, 0,1,1, 12, details=FALSE)
par(mfrow=c(2,2))
  plot(ms2$fit$residuals,type="l")
  acf(ms2$fit$residuals,lag.max = 400)
  P2=spectrum(ms2$fit$residuals,log="no")
    1/P2$freq[which.max(P2$spec)]
  LBQPlot(ms2$fit$residuals,lag.max=400,StartLag=1)

fit2 <-Arima(SST, order=c(0,1,26), seasonal=list(order=c(0,1,1),period=12))
par(mfrow=c(3,2))
  acf(fit2$residuals,lag.max = 40)
  acf(fit2$residuals,lag.max = 400)
  LBQPlot(fit2$residuals,lag.max=40,StartLag=1)
  LBQPlot(fit2$residuals,lag.max=400,StartLag=1)
  plot(fit2$residuals)
  P2=spectrum(fit2$residuals,log="no")
    1/P2$freq[which.max(P2$spec)]
  plot(fit2$x,col="black",type="l") #original
  lines(fitted(fit2),col="red") #fitted

dev.off()
sarima.for(SST,5,0,1,26, 0,1,1,12)#sarima forecast 5 points

dev.off()
#Two Time Series analysis
ccf(SST,SOI,30, main="SST vs SOI", ylab="CCF")
B = c(0,1)
soi_1=filter(SOI,B,sides=1) #shift 1 lag back based on CCF
soi_1=soi_1[2:406]
sst_1=SST[2:406] 

#Granger test
library(vars)
x = cbind(sst_1,soi_1)
var <- VAR(x, p = 30, type = "const")
causality(var, cause = "soi_1")$Granger

# AUTOMATICALLY SEARCH FOR THE MOST SIGNIFICANT RESULT
for (i in 1:30)
{
  print(causality(VAR(x, p = i, type = "const"), cause = "soi_1")$Granger)
  print(AIC(VAR(x, p = i, type = "const")))
}

# the other direction:
for (i in 1:30)
{
  print(causality(VAR(x, p = i, type = "const"), cause = "sst_1")$Granger)  
}

#ARIMAX
mx = Arima(sst_1, order=c(0,1,26), xreg=soi_1 )
mx
par(mfrow=c(2,2))
  acf(mx$residuals,lag.max = 400)
  LBQPlot(mx$residuals,lag.max=400,StartLag=1)
  plot(mx$residuals)
  px = spectrum(mx$residuals,log="no")
  #plot(mx$x,col="black",type="l") #original
  #lines(fitted(mx),col="red") #fitted
  1/px$freq[which.max(px$spec)]


dev.off()    
#SARIMAX  
msx = sarima(sst_1, 0,1,26, 0,1,1, 12, xreg=soi_1, details=FALSE )
msx = Arima(sst_1, order=c(0,1,26), seasonal=list(order=c(0,1,1),period=12), xreg=soi_1)
msx
par(mfrow=c(2,2))
  acf(msx$residuals,lag.max = 400)
  plot(msx$residuals)
  psx = spectrum(msx$residuals,log="no")
  LBQPlot(msx$residuals,lag.max=400,StartLag=1)
  #hsx = spectrum(msx$fit$residuals,spans=c(20,20))
  
par(mfrow=c(2,1))
  plot(fit2$x,col="black",type="l") #original
  lines(fitted(fit2),col="red",lty=2) #fitted  
  plot(msx$x,col="black",type="l") #original
  lines(fitted(msx),col="blue",lty=2) #fitted  
  
  
install.packages("forecast")
library(forecast)


#VAR model
install.packages("vars")
library(vars)
x = cbind(SST,SOI)
summary(VAR(x, p=1))
# p = 1 indicates lag = 1
summary(fit2<-VAR(x, p=2, type="both"))
#  type = c("const", "trend", "both", "none")
acf(resid(fi2t), 52)
prediction = predict(fit, n.ahead = 24, ci = 0.95)
fanchart(prediction)
\end{document}


