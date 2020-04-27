library(forecast)
library(ggplot2)
library(readxl)
datosGas <- read_excel("C:/Users/rosas/Desktop/datosGas.xlsx")

#xag<-read.csv("MES_100.csv")
#XAG<-xag$x  #serie mensual
#ggplot(data = xag, aes(factor(fechas),  x))+geom_point()
#geom_line(color = "#00AFBB", size = 2)
datosGas_ts<-ts(datosGas$Total, frequency = 12, start=c(2005,12))
length(datosGas_ts)
plot.ts(datosGas_ts, col='blue',  xlab = " ", ylab = "Volumen de ventas (miles de barriles)" )

#lambda para estabilizar la varianza
BoxCox.lambda(datosGas_ts, method = c("guerrero", "loglik"), lower = -1, upper = 2) #lambda = -1
#serie transformada XAG^-1
datosGas_ts_t<- (1/(datosGas_ts))
length(datosGas_ts_t)
plot.ts(datosGas_ts_t, col='blue')
Gas_d1<-diff(datosGas_ts_t)#primera diferencia a la serie transformada
length(Gas_d1)
plot.ts(Gas_d1,col='blue')
Gas_d2<-diff(Gas_d1)  #segunda diferencia a la serie transformada
length(Gas_d2)
plot.ts(Gas_d2,col='blue')
Gas_d3 <- diff(Gas_d2) #tercera diferencia a la serie transformada
plot.ts(Gas_d3,col='blue')

#desviaciones estandar de las diferencias
sd(datosGas_ts_t)
sd(Gas_d1)   #minima, por lo tanto trabajaremos con 1 diferencia
sd(Gas_d2)
sd(Gas_d3)
#concluimos que la diferencia 2 tiene la menor varianza muestral

#funciones de autocorrelacion
# Para dos diferencias la FAC y FACP esta dada por
Acf(Gas_d2, lag.max = 24, type = "correlation", plot = TRUE, na.action = na.contiguous, demean = TRUE)
Pacf(Gas_d2, lag.max = 24, plot = TRUE, na.action = na.contiguous, demean = TRUE)

# De la FAC nos quedamos con un MA(1) y EMA(2) [12]
# De la FACP nos quedamos con un AR(2) y ERA(1) [12]

# Por lo tanto nuestra primera propuesta seria un modelo ARMA(2,2,1)(1,0,2)[12]

###############


# modelo ARMA(2,2,1)(1,0,2)[12]
modelo1 <- arima(datosGas_ts_t, order = c(2, 2, 1), list(order = c(1, 0, 2), period = 12)) 

modelo1<-arima(XAG_t,c(0,1,0))  #es mejor usar este para pronósticos
residuales<-as.vector(modelo1$residuals)
plot.ts(residuales) #no hay residuales significativamente distintos de cero
abline(h=2*sd(residuales), untf=FALSE)
abline(h=-2*sd(residuales), untf=FALSE)
abline(h=3*sd(residuales), untf=FALSE, col='red')
abline(h=-3*sd(residuales), untf=FALSE, col='red')
hist(residuales)
mean(residuales)
sd(residuales)
acf(residuales)
Box.test(modelo1$residuals, lag = 24, type = 'Ljung') #pag 154
Box.test(modelo1$residuals, lag = 11, type = 'Ljung')
# p>0.05 por lo tanto no rechazamos que las autocorrelaciones sean 0

############## train test
XAG_train<-XAG[1:80]  
XAG_train<- ts(XAG_train,frequency=12,start=c(2011,2))
length(XAG_train)
plot.ts(XAG_train)
XAG_test<-XAG_t[81:length(XAG_t)]
XAG_test<-ts(XAG_test,frequency=12,start=c(2017,10)) #octubre/17
plot.ts(XAG_test)
length(XAG_test)
test<-c(rep(0,length(XAG_train)),XAG_test)
test<-ts(test,frequency=12,start=c(2011,2))
plot(test)


XAG_train_t<-1/XAG_train
XAG_train_d1<-diff(XAG_train_t)
XAG_train_d2<-diff(XAG_train_d1)
sd(XAG_train_t)
sd(XAG_train_d1)   #sigue siendo la mínima
sd(XAG_train_d2)
acf(XAG_train_t)# no hay un decaimiento rapido hacia cero
acf(XAG_train_d1)#decaimiento rapido hacia cero, y algunos valores significativamente diferentes de cero
acf(XAG_train_d2) #primera autocrrelacion muy grande

auto.arima(XAG_train_t)

test_mod<-arima(XAG_train_t,c(0,1,0))
resid_test<-as.vector(test_mod$residuals)
hist(resid_test)
mean(resid_test)
sd(resid_test)
plot(resid_test) #no hay residuales significativamente distintos de cero
abline(h=2*sd(resid_test), untf=FALSE)
abline(h=-2*sd(resid_test), untf=FALSE)


#seqplot.ts(XAG_t, fitted(modelo1), colx = "black", coly = "red")
#seqplot.ts(XAG_t, fitted.values(modelo1), colx = "black", coly = "red")
#seqplot.ts(XAG_d1, fitted(modelo2), colx = "black", coly = "red")
seqplot.ts(XAG_train_t, fitted(test_mod), colx = "black", coly = "red", ylab = 'XAG')

predicciones<-forecast(test_mod, h=20, level=c(70))
plot(predicciones,ylim=range(c(0,0.1)), main='')
par(new=TRUE)
plot( XAG_t, col="red",ylim=range(c(0,0.1)) , main='')

#necesarios para hacer los pronosticos en la serie original
errores<-c(predicciones$residuals,XAG_test-predicciones$mean)
var_errores_pred<-c()
for (i in 1:length(XAG_test)) {
  var_errores_pred[i]<-var(errores[1:80+i])
}
#factor de correcion por sesgo lambda=-1
lam=-1
correccion<- 1/(0.5 + 0.5*sqrt(1-2*lam*(lam-1)*((1+lam*predicciones$mean)^(-2))*var_errores_pred))
length(correccion)
tail(correccion)
#predicciones_originales<-predicciones
#predicciones_originales$originales<-InvBoxCox(predicciones$mean, lambda = -1, biasadj=TRUE, fvar = var_errores_pred)
predicciones_originales<-data.frame( correcion*(1/predicciones$mean), correcion*(1/predicciones$lower), correcion*(1/predicciones$upper))
colnames(predicciones_originales)<- c('pred', 'lower', 'upper')

p<-c(XAG_train,predicciones_originales$pred)
low<-c(XAG_train,predicciones_originales$lower)
up<-c(XAG_train,predicciones_originales$upper)
XAG_pred<-ts(p,frequency=12,start=c(2011,2))
XAG_low<-ts(low,frequency=12,start=c(2011,2))
XAG_up<-ts(up,frequency=12,start=c(2011,2))

plot( XAG_pred, col="red",ylim=range(c(0,50)) , main='', ylab='XAG')
par(new=TRUE)
plot( XAG_low, col="green",ylim=range(c(0,50)) , main='', ylab='XAG')
par(new=TRUE)
plot( XAG_up, col="blue",ylim=range(c(0,50)) , main='', ylab='XAG')
par(new=TRUE)
plot(XAG,ylim=range(c(0,50)), main='', ylab='XAG')


media_error<-mean(XAG[81:100]-predicciones_originales$pred)
ECM<-mean( (XAG[81:100]-predicciones_originales$pred)^2 )#error cuadratico medio

media_error_t<-mean(XAG_test-predicciones$mean)
ECM_t<-mean( (XAG_test-predicciones$mean)^2 )#error cuadratico medio


# NO CORRER LO ULTIMO
#####################################
#####################################
#####################################
#####################################
#####################################
#####################################
################
#prueba
############
test_mod_MA<-arima(XAG_train_t,c(2,1,3))
resid_test_MA<-as.vector(test_mod_MA$residuals)
max(resid_test_MA)
hist(resid_test_MA)
sd(resid_test_MA)
resid_g<-c()
k<-1
for (i in 1:length(resid_test_MA)) {
  if( resid_test_MA[i] >= 2*sd(resid_test_MA) ) {
    resid_g[k]<-resid_test_MA[i]
    k<-k+1
  }
}

plot(resid_test_MA) #no hay residuales significativamente distintos de cero
abline(h=2*sd(resid_test_MA), untf=FALSE)
abline(h=-2*sd(resid_test_MA), untf=FALSE)
seqplot.ts(XAG_train_t, fitted(test_mod_MA), colx = "black", coly = "red")
###################

#comb_ts <- cbind(XAG_t, fitted(test_mod))


modelo4<-Arima(XAG_t, order = c(0,1,0))
seqplot.ts(XAG_t, modelo4$fitted, colx = "black", coly = "red")
modelo4$fitted
modelo4$series
modelo4$x
XAG_t - modelo4$x


xag2<-read.csv("mes2.csv")
serie<-ts(xag2$x)
serie<-1/serie
auto.arima(serie)
modelo5<-arima(serie, order = c(0,1,0) )
modelo6<-Arima(serie, order = c(0,1,0) )
modelo6$fitted
seqplot.ts(serie, modelo6$fitted, colx = "black", coly = "red")
seqplot.ts(serie, fitted(modelo5), colx = "black", coly = "red")
