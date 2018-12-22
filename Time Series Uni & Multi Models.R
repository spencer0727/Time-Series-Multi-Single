

This analysis is used to describe, a univariate and multivairate model using the Box-Jenkins methods for time series analysis. This is written as an R Markdown file. 

A few comments first, I will describe only a full process  for the univariate analysis.
ADF Differenced – No Unit Root
ADF Normal – Unit Root
ARIMA – 1,1,1
Residual Independence – Yes  


To begin the analysis, I have chosen to perform a univariate analysis on the following currency pair. The Brazilian Real, Mexican Peso and the United States Dollar. Due to the plethora of available data but also due to the influence the Dollar has in global economics. All the data I gathered is from FRED (Federal Reserve Economic Data).

```{r}
rm(list=ls())
library(quantmod)
library(tseries)
library(forecast)
library(AER)
library(tinytex)
library(mFilter)
library(readxl)
library(dynlm)
library(urca)


BRLUSD <- read_excel("TIme Series/TERM PAPER/S/BRLUSD.xls")
BRLUSD <- BRLUSD[,2]
BRLUSD <-ts(BRLUSD, start = 2000, frequency = 12)



```


 The data imported is monthly data so a frequency of 12 is needed. 


```{r}
plot(BRLUSD, ylab = 'BRLUSDFOREX')
plot(diff(BRLUSD), ylab = 'BRLUSDFOREX')
```


Next, I plotted both a graph of both the data but also of the differenced series just to see if there were any obvious trends that need to be accounted for. Maybe there is short term trends although it is not so obvious in order to say yes or no. 


Now, we need to check if a random walk is present we do this first by examining an auto correlation of the differenced series. The max lag number is chose is a little arbitary although it does not really seem to make a difference 

```{r}
acf(diff(BRLUSD))
acf(BRLUSD, lag.max = 200)
```


The results do not appear to be very conclusive. Even if we use the lagged data. Also, there appears to be a trend that 
begins downward, reverses, and then continues to go down towards the end of the series. Again, the trend does not appear to be completely random and does look deterministic. 

Now I will use an ADF test, if not stationary, we will take the difference series until stationary thus the series will be integrated d times. The process I am using is the ARIMA(p,d,q)(Box-Jenkins model). 

In order to begin analyzing this time series data, it important to first see if can identify any trends and if so how can we deal with them. This will show if a unit root is present, if so we can estimate the first differences




```{r}
(ur.df(BRLUSD,type = "trend", lags = 12,selectlags = "AIC")@cval)
(ur.df(BRLUSD,type = "trend", lags = 12,selectlags = "AIC")@teststat)
(summary(ur.df(BRLUSD,type = "trend",lags = 12,selectlags = "AIC")))
```

From the results, the test statistic was -1.437024 which is larger than the 95% -3.43 which shows that at this level we can not reject the null of a being present unit root. 

However, the critical values for phi3 is 6.49, and from the test we yield 1.2243 which is smaller than thecritical value, hence, we can safely conclude that it is not a drift and hence a unit root is present. The unit root and joint hyothesis indicates that maybe a trend is present. 
Thus, we can now check for a unit root with respect to drift. 


```{r}
ur.df(BRLUSD ,type = "drift",lags = 12,selectlags = "AIC")@cval
ur.df(BRLUSD,type = "drift",lags = 12,selectlags = "AIC")@teststat
summary(ur.df(BRLUSD,type = "drift",lags = 12,selectlags = "AIC"))
```

Clearly idicated from the results we cannot reject the null hypothesis. The critical values for the test were found to be large in order to reject the null hypothesis at a 95% confidence interval. 

However, the critical values for tau2 and phi1 at a 95% confidence level is -2.88 and 4.63 respecitvely, from the results we yielded -1.3137 and 1.1413. Thus there is still a unit root.Next we estimate the OLS excluding a trend or drift,

```{r}
ur.df(BRLUSD,type = "none",lags = 12,selectlags = "AIC")@cval
ur.df(BRLUSD,type = "none",lags = 12,selectlags = "AIC")@teststat
summary(ur.df(BRLUSD,type = "none",lags = 12,selectlags = "AIC"))
```

Our results finally show with a test statistic of 0.4005 which is larger than the required -1.95 that a unit root is present therefore, we can now take the first difference and run the ADF tests once more. Because we have found a unit root we know that our data is not stationary. 

Starting with a differenced ADF with respect to trend
```{r}
ur.df(diff(BRLUSD,type = "trend", lags = 12,selectlags = "AIC"))@cval
ur.df(diff(BRLUSD,type = "trend", lags = 12,selectlags = "AIC"))@teststat
summary(ur.df(diff(BRLUSD,type = "trend",lags = 12,selectlags = "AIC")))

```


From the results we get a test statistic of -8.3616 which is well below the required level of -1.95 which is the thereshold for the 5 percent level.

Next we can run the differenced test with respect to the drift

```{r}
ur.df(diff(BRLUSD,type = "drift", lags = 12,selectlags = "AIC"))@cval
ur.df(diff(BRLUSD,type = "drift", lags = 12,selectlags = "AIC"))@teststat
summary(ur.df(diff(BRLUSD,type = "drift",lags = 12,selectlags = "AIC")))
```

We are able to reject the null hypothesis of a unit root present with a test stat of -8.3616 which is below the required -1.95 for the 95% confidence interval testing.

Finally we can run the differenced ADF test, with neither a trend nor a drift.

```{r}
ur.df(diff(BRLUSD,type = "none", lags = 12,selectlags = "AIC"))@cval
ur.df(diff(BRLUSD,type = "none", lags = 12,selectlags = "AIC"))@teststat
summary(ur.df(diff(BRLUSD,type = "none",lags = 12,selectlags = "AIC")))

```

we can see that the test stat is -8.3616 which is below the required -1.95 at a 95% confidence interval. We see that there is no unit root in the first diffrence BRLUSD is I(1) and requires differencing once to make it stationary. 




Next we for the final pieces of the ARIMA model we need to find the values of p and q. For an AR(p) process we expect ACF to go down with the PACF going up with increases in p. For an MA(q) process we expect PACF to go down with the PACF going up with increases in q. 

```{r}
par(mfrow=c(1,2))
taperedacf(diff(BRLUSD),lag.max = 12, bty = 'l', nsim = 100, plot = TRUE, level = 95)
taperedpacf(diff(BRLUSD),lag.max = 12, bty = 'l', nsim = 100, plot = TRUE, level = 95)
par(mfrow=c(1,1))
```

From the results it is clear that there is a significant spike in PACF and ACF around the first lag so using the acf and pacf functions, both have significant peaks at the first lag.Thus, AR(1) and MA(1) would be implied. Therefore, and ARIMA model of order (1,1,1) would make the most sense or an ARIMA of (1,0,1) of the differenced series. 



Now that we have decided a that an ARIMA model of (1,1,1) makes the most sense we should now estimate the parameters of the model. After defining the function we can now run an ACF funciton of the residuals to check for correlation.

```{r}
resBRLUSD <- arima((BRLUSD), order = c(1,1,1), method = "ML")
Acf(residuals(resBRLUSD))
plot(residuals(resBRLUSD))
```

From the results we can see that the residuals are not autocorrelated and are indepdent. Finally I need to perform the Box Test to see if residuals are truely independent. With the degrees of freedom set to 2 becuase both p=1 and q=1. 

```{r}
Box.test(residuals(resBRLUSD), fitdf=2, lag=12,  type="Ljung")
```
Because we have a p value that is larger than 0.05, the results yield a p-value of 0.2678 thus the residuals are independent. It looks to like white noise with respect to the residuals. 


Do a Jarque Bara test to see if normal. The results are a p-value = 7.475e-10, therefore further support for non-normality which is more conducive for white noise. 
```{r}
jarque.bera.test(resBRLUSD$residuals)
```

Therefore we can conclude that the model is an Arima is (1,1,1)



MULTI-VARIATE

For my multi variate analysis I chose to perform examine the relationship between the Swiss Franc/United States Dollar, Mexican Peso/United States Dollar, and the United States Federal Funds rate for the first model. For the second model I chose to examine the relationship between the Swiss Franc/United States Dollar and the Mexican Peso/United States Dollar. 

I thought I would be interesting to examine the influence of the United States Federal Funds Rate to trading partners. I,e, is there any change in relationship between Switzerland and Mexico from the United States Federal Reserve. 

```{r}
rm(list=ls())
library(quantmod)
library(tseries)
library(urca)
library(dynlm)
library(forecast)
library(AER)
library(tinytex)
library(mFilter)
library(readxl)
library(vars)

FEDFUN <- read_excel("TIme Series/TERM PAPER/S/FEDFUNDS.xls")
MXNUSD <- read_excel("TIme Series/TERM PAPER/S/MXNUSD.xls")
CHFUSD <- read_excel("TIme Series/TERM PAPER/S/CHFUSD.xls")

FEDFUN <- FEDFUN[,2]
FEDFUN <-ts(FEDFUN, start = 2000, frequency = 12)


MXNUSD <- MXNUSD[,2]
MXNUSD <-ts(MXNUSD, start = 2000, frequency = 12)

CHFUSD <- CHFUSD[,2]
CHFUSD <-ts(CHFUSD, start = 2000, frequency = 12)
```


First, I need to import my libraries and data. After I need to transform my data into time series data. 

Next, I will check for cointegration between the series. If cointegration is found than it would be very important to run an error correction model in order to sort out any problems in the data. 

```{r}
po.test(cbind(CHFUSD,MXNUSD,FEDFUN))
```

From the results of the Phillips-Ouliaris Cointegration Test, it is clear that the p-value is 0.15 which is above the required 0.05 for a 95% confidence interval to reject the null hypothesis of cointegration between the variables. This is an important step in the process to make sure that the results are not spurious.

Now that we have ruled out cointegration, now we can perform a VAR test in order find the relationship between regression models. Fist we do a VAR select test in order to select the number of lags. Then a coeftest gives you the values of the coefficients. And then finally a serial test in order to check for the autocorrelation of the residuals

```{r}
var.data <- ts.intersect(diff(CHFUSD), diff(MXNUSD), diff(FEDFUN))
VARselect(var.data)$selection
var.sw <- VAR(var.data, p=2)
coeftest(var.sw)

serial.test(var.sw)
```


From the results we can see that at a 95% confidence interval that the Swiss Franc/United States exchange rate affects the United States Federal Funds rate at of lag of 1. I found this fairly odd considering the fact that I believed the relationship should be the other way around. 

Additionally, the P-Value is very high of 0.1819 so therefore we cannot reject the null hypothesis of no autocorrelation in the residuals.


For the second model I chose to examine the relationship between the Swiss Franc/United States Dollar and the Mexican Peso/United States Dollar.


First lets import and transform the data.

```{r}
rm(list=ls())
library(quantmod)
library(tseries)
library(urca)
library(dynlm)
library(forecast)
library(AER)
library(tinytex)
library(mFilter)
library(readxl)
library(vars)

CHFUSD <- read_excel("TIme Series/TERM PAPER/S/CHFUSD.xls")
MXNUSD <- read_excel("TIme Series/TERM PAPER/S/MXNUSD.xls")

MXNUSD <- MXNUSD[,2]
MXNUSD <-ts(MXNUSD, start = 2000, frequency = 12)

CHFUSD <- CHFUSD[,2]
CHFUSD <-ts(CHFUSD, start = 2000, frequency = 12)
```

Again we can now check for cointegration using the  Phillips-Ouliaris Cointegration Test. 



```{r}
po.test(cbind(CHFUSD,MXNUSD))
```

Again we find that the P-Value is above 0.05 therefore no cointegration exists.

Now that we have ruled out cointegration, lets perform a VAR test in order find the relationship between regression models. Fist we do a VAR select test in order to select the number of lags. Then a coeftest gives you the values of the coefficients. And then finally a serial test in order to check for the autocorrelation of the residuals.

```{r}
var.data <- ts.intersect(diff(CHFUSD), diff(MXNUSD))
VARselect(var.data)$selection
var.sw <- VAR(var.data, p=2)
coeftest(var.sw)

serial.test(var.sw)
```



The results of the test show that there is a significant relationship between the CHFUSD/MXNUSD of the second lag at a 5% confidence interval. Additionally the P-Value is 0.3092 which is very high therefore we cannot the reject null hypothesis of no autocorrelation. 
