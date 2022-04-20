# TimeSeries_SST_2016
Final class project for STAT 650 - Time series analysis. 

## A Time Series Analysis Study on Sea Surface Temperature (SSTs) Prediction

### Abstract

A 34-year monthly sea surface temperature (SST) time series dataset is used in order to find out a best model that fits it and can make reasonable prediction. Multiple AR, MA, ARMA, ARIMA, ARIMAX, seasonal ARIMA and seaonal ARIMAX models are formulated and tested with di↵erent terms and moments.
Another time series data southern oscillation index (SOI) is utilized as the external regressor in the ARIMAX and SARIMAX models. Based on the AIC and BIC
values calculated for each model and the diagnosis analyses, the best fitting model is determined as SARIMA(0,1,26)(0,1,1)12 model with a AIC value of 19.18. The model is then fitted and compared with real data and the fittng curve shows a good agreement. Though the model does not precisely predict the abnormal increase in temperature caused by the extreme El Ni˜no event occurred in last year, the real data curve is still within the confidence interval.


Code can be found in the `final_code.r` file, and figures and research report is `STAT650FinalProject.pdf`.
