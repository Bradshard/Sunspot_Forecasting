data_sunspot = 
  read.table("E:/2021-2022 Fall - Spring Books/project/Monthly_mean_total_sunspot_number.txt",
             header= FALSE,
             stringsAsFactors = FALSE)

names(data_sunspot) <- c("Year", "Month",
                        "PercentYearCompleted",
                        "MMSN",
                        "MonthlyMeanStandardDeviation",
                        "NumberOfObservation")

#sum(data_sunspot$MMSN[which(data_sunspot$Year == 1995)])

#data_sunspot$MMSN[which(data_sunspot$Year == 1995 & data_sunspot$Month == 5)]

x = unique(data_sunspot$Year[which(data_sunspot$Year == 1995)])
yearly_data = c()
years = c()
while(x <= 2021){
  for(i in 1:max(data_sunspot$Month)){
    yearly_data = c(yearly_data,data_sunspot$MMSN[which(data_sunspot$Year == x & data_sunspot$Month == i)])
    years = c(years,x)
  }
  x = x + 1
}

new_Data = as.data.frame(cbind(yearly_data,years))
names(new_Data) = c("MMSN","Year")

head(new_Data$MMSN)


library(ggplot2)


p <- ggplot(new_Data, aes(x = Year, y = MMSN)) +
  geom_hex(binwidth = c(0.05, 0.3)) +
  geom_point(size = 2, shape = 4)+
  geom_smooth(color = "grey10", method = "loess", span = 0.5) +
  scale_fill_distiller(palette = "BrBG", direction = 1) +
  ggtitle(label = "Mean Sunspot Number vs Year",
          subtitle = "Year by Year") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p8 = p + theme_classic()

p8
## year sunspot monthly data.


# Detrended with linear filter
lin.mod <- lm(new_Data$MMSN ~ time(new_Data$Year))
lin.trend <- lin.mod$fitted.values
linear <- ts(lin.trend,start = c(1995,1), frequency = 12)
lin.cycle <- new_Data$MMSN - linear

linear

par(mfrow = c(1,2), mar = c(2.2,2.2,1,1), cex = 0.8)
# plot two graphs side by side and squash margins
plot.ts(new_Data$MMSN, ylab = "MMSN")
lines(linear,col = "red")
legend("topleft", legend = c("data","trend"), lty = 1,
       col = c("black", "red"), bty = "n")
plot.ts(lin.cycle, ylab = "")
legend("topright", legend = c("cycle  "), lty = 1, col = c("black"),
       bty = "n")


# Detrend data with the Hodrick-Prescott filter

library(mFilter)

hp.decom <- hpfilter(new_Data$MMSN, freq = 129600, type = "lambda") # 129600 from the lambda ^4 * 6.25 
# 6.25 is from lambda annual = 1/4^4 * 1600
hp.decom.cycle <- ts(hp.decom$cycle,start = c(1995,1), frequency = 12)
new_Data.MMSN <- ts(new_Data$MMSN,start = c(1995,1), frequency = 12)
hp.decom.trend <- ts(hp.decom$trend, start = c(1995,1), frequency = 12)


par(mfrow = c(2, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(new_Data.MMSN, ylab = "") # plot time series
lines(hp.decom.trend, col = "red") # include HP trend
legend("topleft", legend = c("data", "HPtrend"), lty = 1,
       col = c("black", "red"), bty = "n")
plot.ts(hp.decom.cycle, ylab = "") # plot cycle

legend("topleft", legend = c("HPcycle"), lty= 1, col = c("black"),
       bty = "n")

length(hp.decom.cycle) # It is to know how many parts to separate the system into.

# Rescaled Range Analysis and Hurst Exponent

library("pracma")

hurst_r_s_measure <- c()
log_n <- c()
hurst_a <- c()

for(i in 5:20){
  hurst <- hurstexp(hp.decom.cycle, d = i, display = TRUE)
  hurst_r_s_measure = c(hurst_r_s_measure,hurst$Hal)
  hurst_a = c(hurst_a, hurst$Hrs)
  log_n = c(log_n, log(i))
}

log_log <- data.frame()

for(i in 1:length(log_n)){
  log_log[i,1] <- c(hurst_r_s_measure[i])
  log_log[i,2] <- c(log_n[i])
}


linear_model <- lm(V1~V2,data = log_log)
graphics.off()
plot(log_n,hurst_r_s_measure,
     ylab = "Log(R/S)",
     xlab = "Log(n)",
     main = "Hurst Exponent for sunspot",
     col = "purple")

abline(linear_model, lwd = 2, col = "darkgreen")
text(1.8,0.93,"Hsunspot = ", cex =1)
text(2.05,0.93,round(hurst$Hal,2))
linear_model$coefficients

# Simplex Projection Analysis

#install.packages('rEDM')
library(rEDM)

times <- rep(1995:2021, each = 12)
simplex_func <- cbind(times,hp.decom$trend)
simplex_func = as.data.frame(simplex_func)
names(simplex_func) = c("Time","Trend")

results <- Simplex(dataFrame = simplex_func, 
        lib = "1 162", pred = "162 324", E = 2, Tp = 1, knn = 3, tau = -1, 
        exclusionRadius = 0, columns = "Trend", target = "Trend", embedded = FALSE,
        generateSteps = 228, parameterList = FALSE, showPlot = TRUE)

results.pred <- na.omit(results$Predictions)
results.pred <- ts(results.pred,start = c(2021,1), frequency = 12)

plot(hp.decom.trend, xlim = c(1995,2040))
lines(results.pred)

# Holt's Two-parameter Exponential Smoothing
library(aTSA)

Holt(hp.decom.trend, type = "additive", alpha = 0.20, beta = 0.1057, lead= 228, damped = TRUE, phi = 0.98, plot = TRUE)


# Autoregressive Integrated Moving Average (ARIMA)
library(forecast)
fit_arima <- arima(hp.decom.trend, order = c(3,1,1) ,seasonal = list(order = c(3,1,1), period = 36))
checkresiduals(fit_arima)
autoplot(forecast(fit_arima, h = 48))
plot(forecast :: forecast(fit_arima, h = 96), xlab= "Years", ylab = "Sunspots", main = "Forecast of Sunspot Numbers in 25th cycle")

max_list <- c()

for(i in 11:70){
  max_list <- c(max_list,log((i-11)/(71-i)))
}
text(2025,150,max(exp(max_list)),cex = 0.8)
text(2019,150,"Max expected is:", cex = 0.8)


## https://stackoverflow.com/questions/40302029/how-to-specify-minimum-or-maximum-possible-values-in-a-forecast

# F10.7 cm Index
data_f10_7 = 
  read.table("E:/2021-2022 Fall - Spring Books/project/F10.7_cm_index.txt",
             header= FALSE,
             stringsAsFactors = FALSE)

names(data_f10_7) <- c("Year", "Month",
                         "Index_F10.7")

#sum(data_sunspot$MMSN[which(data_sunspot$Year == 1995)])

#data_sunspot$MMSN[which(data_sunspot$Year == 1995 & data_sunspot$Month == 5)]

x1 = unique(data_f10_7$Year[which(data_f10_7$Year == 1995)])
yearly_data1 = c()
years1 = c()
while(x1 <= 2021){
  for(i in 1:max(data_f10_7$Month)){
    yearly_data1 = c(yearly_data1,data_f10_7$Index_F10.7[which(data_f10_7$Year == x1 & data_f10_7$Month == i)])
    years1 = c(years1,x1)
  }
  x1 = x1 + 1
}

new_Data1 = as.data.frame(cbind(yearly_data1,years1))
names(new_Data1) = c("F10.7_cm_index","Year")

head(new_Data1$F10.7_cm_index)


library(ggplot2)


p1 <- ggplot(new_Data1, aes(x = Year, y = F10.7_cm_index)) +
  geom_hex(binwidth = c(0.05, 0.3)) +
  geom_point(size = 2, shape = 4)+
  geom_smooth(color = "grey10", method = "loess", span = 0.5) +
  scale_fill_distiller(palette = "BrBG", direction = 1) +
  ggtitle(label = "F10.7_cm vs Year",
          subtitle = "Year by Year") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p9 = p1 + theme_classic()

p9

## year by yer F10.7 monthly data.



# Detrended with linear filter
lin.mod1 <- lm(new_Data1$F10.7_cm_index ~ time(new_Data1$Year))
lin.trend1 <- lin.mod1$fitted.values
linear1 <- ts(lin.trend1,start = c(1995,1), frequency = 12)
lin.cycle1 <- new_Data1$F10.7_cm_index - linear1

linear1

par(mfrow = c(1,2), mar = c(2.2,2.2,1,1), cex = 0.8)
# plot two graphs side by side and squash margins
plot.ts(new_Data1$F10.7_cm_index, ylab = "F10.7 cm Index")
lines(linear1,col = "red")
legend("topleft", legend = c("data","trend"), lty = 1,
       col = c("black", "red"), bty = "n")
plot.ts(lin.cycle1, ylab = "")
legend("topright", legend = c("cycle  "), lty = 1, col = c("black"),
       bty = "n")


# Detrend data with the Hodrick-Prescott filter

library(mFilter)

hp.decom1 <- hpfilter(new_Data1$F10.7_cm_index, freq = 129600, type = "lambda") # 129600 from the lambda ^4 * 6.25 
# 6.25 is from lambda annual = 1/4^4 * 1600
hp.decom.cycle1 <- ts(hp.decom1$cycle,start = c(1995,1), frequency = 12)
new_Data.F10.7_cm_index <- ts(new_Data1$F10.7_cm_index,start = c(1995,1), frequency = 12)
hp.decom.trend1 <- ts(hp.decom1$trend, start = c(1995,1), frequency = 12)


par(mfrow = c(2, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(new_Data.F10.7_cm_index, ylab = "") # plot time series
lines(hp.decom.trend1, col = "red") # include HP trend
legend("topleft", legend = c("data", "HPtrend"), lty = 1,
       col = c("black", "red"), bty = "n")
plot.ts(hp.decom.cycle1, ylab = "") # plot cycle

legend("topleft", legend = c("HPcycle"), lty= 1, col = c("black"),
       bty = "n")

length(hp.decom.cycle1) # It is to know how many parts to separate the system into.



# Rescaled Range Analysis and Hurst Exponent

library("pracma")

hurst_r_s_measure1 <- c()
log_n1 <- c()
hurst_a1 <- c()

for(i in 5:20){
  hurst1 <- hurstexp(hp.decom.cycle1, d = i, display = TRUE)
  hurst_r_s_measure1 = c(hurst_r_s_measure1,hurst1$Hal)
  hurst_a1 = c(hurst_a1, hurst1$Hrs)
  log_n1 = c(log_n1, log(i))
}

log_log1 <- data.frame()

for(i in 1:length(log_n1)){
  log_log1[i,1] <- c(hurst_r_s_measure1[i])
  log_log1[i,2] <- c(log_n1[i])
}


linear_model1 <- lm(V1~V2,data = log_log1)
graphics.off()
plot(log_n1,hurst_r_s_measure1,
     ylab = "Log(R/S)",
     xlab = "Log(n)",
     main = "Hurst Exponent for F10.7 cm Index",
     col = "purple")

abline(linear_model1, lwd = 2, col = "darkgreen")
text(1.8,0.92,"HF10.7_cm_Index = ", cex =0.7)
text(2.1,0.92,round(hurst1$Hal,2), cex = 0.7)
linear_model1$coefficients


# Autoregressive Integrated Moving Average (ARIMA)
library(forecast)
fit_arima1 <- arima(hp.decom.trend1, order = c(1,0,0) ,seasonal = list(order = c(2,0,2), period = 36))
checkresiduals(fit_arima1)
autoplot(forecast(fit_arima1, h = 48))
plot(forecast :: forecast(fit_arima1, h = 96), xlab= "Years", ylab = "F10.7 cm Index", main = "Forecast of F10.7 cm Index in Future Years")

max_list1 <- c()

for(i in 550:1600){
  max_list1 <- c(max_list1,log((i-550)/(1601-i)))
}
text(2025,1402,max(exp(max_list1)),cex = 0.8)
text(2021.7,1400,"Max expected is:", cex = 0.8)

# accuracy of the forecast
accuracy(forecast(fit_arima1, h = 96))


## https://stackoverflow.com/questions/40302029/how-to-specify-minimum-or-maximum-possible-values-in-a-forecast
## https://stats.stackexchange.com/questions/194453/interpreting-accuracy-results-for-an-arima-model-fit