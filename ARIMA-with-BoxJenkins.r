# ________________________________________________________________________
# ______ METODOLOGÍA BOX JENKINS EN EL PRONÓSTICO DE SERIES TEMPORALES ___
# _____________________ USANDO UN MODELO ARIMA(p,d,q) ____________________
# ________________________________________________________________________
#
#
# 1. INTRODUCCIÓN
# En el análisis de series de tiempo, la predicción futura de los valores
# de una serie es una tarea crucial en la toma de decisiones en diversas
# áreas, como la economía, la finanzas, la ingeniería, entre otras. En este
# trabajo, utilizamos la metodología Box-Jenkins para realizar el pronóstico
# 10 pasos adelante de la serie "Ingresos personales y gastos personales de
# Estados Unidos", desestacionalizado y obtenido de la base de datos de la
# Reserva Federal de St. Louis (FRED).. Esta metodología consta de cuatro
# pasos: identificación, estimación, validación y pronóstico.
#
#
#
#
#
#
# ------- 2. SELECCIÓN DE LA SERIE Y PREPARACIÓN DE LOS DATOS -----------

# Limpiamos las variables de entorno
rm(list = ls())

# Librerías
library(tidyverse)
library(dplyr)
library(forecast)
library(tseries)
library(lubridate)
library(xts)
library(ggplot2)
library(gridExtra)
library(stargazer)
library(timeSeries)
library(aTSA)
library(lmtest)
library(fable)
library(seasonal)
library(feasts)
library(urca)
library(psych)
library(astsa)
library(car)
library(tsibble)

# Input
ticker <- "PSAVERT"

# Obteniendo los datos
path <- "C:\\Users\\...\\" %>% paste0(ticker, ".csv")
table_raw <- read.csv(path)
head(table_raw, n <- 5)
view(table_raw)

# Configuramos la clase de la variable DATE
table_raw$DATE <- as.Date(table_raw$DATE, format = "%Y-%m-%d")
class(table_raw$DATE)

# Creamos objetos de clase ts y xts con nuestros datos
ts <- ts(table_raw$PSAVERT, start = c(1959, 01, 01), frequency = 12)
xts <- xts(table_raw$PSAVERT , order.by = table_raw$DATE)
view(ts)
class(ts)
class(xts)

# Graficando
start_year <- min(year(table_raw$DATE))
end_year <- max(year(table_raw$DATE))

x11()
ggplot(data = data.frame(ts = ts)) +
  geom_line(aes(x = 1:length(ts), y = ts), lty = 1, lwd = 1, col = "#FF4040") +
  labs(x = "Date", y = ticker, title = paste0("Time series ", ticker,
  " (", start_year, " - ", end_year, ")"))
ggsave(paste0("Serie original ", ticker, ".png"))


# -------------- 3. APLICACIÓN DE LA METODOLOGÍA BOX-JENKINS ---------------

# - 3.1. IDENTIFICACIÓN DEL MODELO -----------------------------------------
# 3.1.1. Visualización de la FAC y FACP
lags <- 24
x11()
grid.arrange(
  ggAcf(ts, lag.max = lags, plot = TRUE, lwd = 1) +
  ggtitle(paste0("FAC ", ticker)),

  ggPacf(ts, lag.max = lags, plot = TRUE, lwd = 1) +
  ggtitle(paste0("FACP ", ticker)),

  nrow = 1, ncol = 2
)


# La FACP tiene un decaimiento lento a 0, por lo que podemos concluir que no
# es una serie estacionaria, por lo que procedemos a estacionalizarla

# 3.1.2. Estacionarización de la serie
# Usaremos los 3 métodos más usados para estacionarizar series:
# a) Diferencia
diff.ts <- diff(ts)

# b) Logaritmo
log.ts <- log(ts)

# c) Diferencia al logaritmo (aproximación del crecimiento)
difflog.ts <- diff(log(ts))*100


# Con el objetivo de continuar el estilo de gráficas que hemos presentado,
# usaremos la librería ggplot2. Para graficar las 4 series, vamos a convertir
# estos timeseries a dataframes
ts_df <- data.frame(date = index(ts), value = coredata(ts))
diff_df <- data.frame(date = index(diff.ts), value = coredata(diff.ts))
log_df <- data.frame(date = index(log.ts), value = coredata(log.ts))
difflog_df <- data.frame(date = index(difflog.ts), value = coredata(difflog.ts))

# Graficamos la serie original, la serie diferenciada, en logarirmo y con la
# diferencia del logaritmo
x11()
grid.arrange(
  ggplot(data = ts_df, aes(x = date, y = value)) +
    geom_line(color = "#FF4040", size = 1) +
    labs(title = paste0("Serie original ", ticker, start_year,
    " - ", end_year)),

  ggplot(data = diff_df, aes(x = date, y = value)) +
    geom_line(color = "deepskyblue4", size = 1) +
    labs(title = paste0("Variación ",  ticker, start_year,
    " - ", end_year)),

  ggplot(data = log_df, aes(x = date, y = value)) +
    geom_line(color = "darkgoldenrod3", size = 1) +
    labs(title = paste0("Logaritmo ",  ticker, start_year,
    " - ", end_year)),

  ggplot(data = difflog_df, aes(x = date, y = value)) +
    geom_line(color = "#6E8B3D", size = 1) +
    labs(title = paste0("Tasa de crecimiento ",  ticker, start_year,
    " - ", end_year)),

  nrow = 2, ncol = 2
)

# A simple vista, usando los métodos de estacionarización, podríamos dudar sobre
# la varianza constante de la serie, así que usaremos un test aumentado de
# Dicky-Fuller para verificar estacionaridad

adf.test(diff.ts) # Si p-value<0.05 acepto la hipótesis alternativa: la serie es
# estacionaria. La serie es estacionaria porque ya hemos realizado un método de
# estabilización al calcular los retornos: la diferencia.

# Graficamos la FAC y FACP de la serie estacionarizada
lags <- 24
x11()
par(mfrow = c(1, 2))
acf(diff.ts, lag.max = lags, plot = TRUE, lwd = 2, xlab = "",
    main = "ACF de la tasa de crecimiento del PIB") 
pacf(diff.ts, lag.max = lags, plot = TRUE, lwd = 2, xlab = "",
    main = "PACF de la tasa de crecimiento del PIB")


# 3.1.3. Criterios de información

# Definimos los valores p q máximos para hacer la identificación finita. De
# acuerdo al criterio de Parsimonia, no elegiremos p q muy grandes
pmax <- 6
qmax <- 6

# Creamos una función para generar un dataframe que nos muestre los ARIMA(p,d,q)
# con todas las posibles combinaciones p q de acuerdo a los máximos que
# anteriormente definimos
criterios_arma_df <- function(ts_object, pmax, qmax, d, bool_trend, metodo){
  index <- 1
  df <- data.frame(p = double(), d = double(), q = double(), AIC = double(),
        BIC = double())
  for (p in 0:pmax) {
    for (q in 0:qmax)  {
      fitp <- arima(ts_object, order = c(p, d, q), include.mean = bool_trend,
                    method = metodo)
      df[index, ] <- c(p, d, q, AIC(fitp), BIC(fitp))
      index <- index + 1
    }
  }
  return(df)
}

# Creamos una función que nos permitirá seleccionar el mejor modelo según AIC
min_AIC <- function(df){
  df2 <- df[which.min(df$AIC), ]
  return(df2)
}

# Creamos una función que nos permitirá seleccionar el mejor modelo según BIC
min_BIC <- function(df){
  df2 <- df[which.min(df$BIC), ]
  return(df2)
}

# Ahora, aplicamos la función criterios_arma_df a la serie estacionarizada.
# Nos mostrará el valor según los criterios AIC y BIC de cada modelo

modelo1 <- criterios_arma_df(diff.ts, pmax, qmax, d = 1, TRUE, "ML")
View(modelo1)

# Buscamos el modelo cuyo valor de AIC sea menor
arma_aic <- min_AIC(modelo1)
# ARIMA(3,1,2)
arma_aic

# Buscamos el modelo cuyo valor de BIC sea menor
arma_bic <- min_BIC(modelo1)
# ARIMA(3,1,2)
arma_bic

# De acuerdo a los resultados de los criterios de información, elegimos nuestro
# mejor ARMA(p,q). Como dato adicional, si el valor mínimo de AIC y el de BIC
# son diferentes, escogeremos el del BIC por criterio de Parsimonia

# Dado que los criterios de información nos aconsejan usar una diferencia,
# pero ya diferenciamos la serie una vez, hallaremos un ARIMA(p,0,q)
modelo0 <- criterios_arma_df(diff.ts, pmax, qmax, d = 0, TRUE, "ML")
View(modelo0)

# Veamos cuál es nuestro posible mejor modelo cuando d=0
aic_modelo0 <- min_AIC(modelo0)
aic_modelo0
# ARMA(4,4)

bic_modelo0 <- min_BIC(modelo0)
bic_modelo0
# ARMA(3,1)

# De acuerdo a los criterios de información, los modelos recomendados son AIC:
# ARMA(4,4) y BIC(3,1). A pesar de que son modelos con p q con 1 o 2 parámetros
# más de los deseados, estos son los modelos que mejor explican la serie, y
# reducir arbitrariamente los p q sería un error que podría llevar a representar
# los datos como un ruido blanco. Por criterio de Parsimonia usaremos el modelo
# recomendado por el criterio de información BIC

# - 3.2. Estimación -----------------------------------------------------------
# Ya hemos identificado el modelo ARMA(p,q) que explica de mejor manera la
# serie, un ARMA(3,1). Ahora realizaremos la estimación del modelo con el método
# de Máxima Verosimilitud (ML).


# Modelamos nuestro ARIMA(3,0,1) para hallar los coeficientes
arima3.0.1 <- arima(diff.ts, order = c(3,0,1), include.mean = T,
              method = "ML")
summary(arima3.0.1)

stargazer(arima3.0.1, type = "text")

# El resultado de la función de verosimilitud logarítmica (Log Likelihood) es
# -714.584, lo que indica que el modelo ajustado es capaz de explicar el 71.45%
# de la varianza de la serie de tiempo.

# - 3.3. Validación ----------------------------------------------------------
# Realizaremos pruebas en los residuales para determinar autocorrelación serial,
# heterocedasticidad y normalidad

# 3.3.1. Autocorrelación serial
# En un modelo ARIMA(p,d,q), suponemos que los errores se comporten como un
# ruido blanco. Queremos que los residuales del modelo: 1) tengan media cero,
# 2) su varianza sea constante, y 3) no presenten autocorrelación serial, es
# decir, que los residuos no dependan del pasado.

# Revisemos los residuales de nuestro modelo para saber si cumple con este
# supuesto

# Calculamos los residuales
error_arima3.0.1 <- residuals(arima3.0.1)
errores2_arima3.0.1 <- error_arima3.0.1^2

# Graficamos los residuales

x11()
grid.arrange(
  ggAcf(error_arima3.0.1, lwd=2) + ggtitle("FAC residuales"),

  ggPacf(errores2_arima3.0.1, lwd=2) + ggtitle("FAC residuales al cuadrado"),

  nrow = 1, ncol = 2
)

# La FAC y FACP nos muestran que el ruido del modelo se comporta como un
# ruido blanco. Aún así, realizaremos el test de Box-Pierce y el test de
# Ljung-Box, que es son pruebas aún más confiables y objetivas

# Dividimos en 4 porque haremos la prueba sobre 1/4 de la muestra
lags.BPtest <- length(diff.ts) / 4
lags.BPtest

# Realizamos el test de Box-Pierce para 1/4 de la muestra, para 20 y para
# 30 rezagos
Box.test(error_arima3.0.1, lag = lags.BPtest, type = c("Box-Pierce"))
Box.test(error_arima3.0.1, lag = 20, type = c("Box-Pierce"))
Box.test(error_arima3.0.1, lag = 30, type = c("Box-Pierce"))
#  Para los 3 casos, el test de Box-Pierce nos arroja un p-value>0, por lo
# que no rechazamos la hipótesis nula H_0 y por el momento se puede asumir
# que no hay autocorrelación en los residuos del modelo.

# Ahora realizaremos el test de Ljung-Box para los mismos rezagos
Box.test(error_arima3.0.1, lag = lags.BPtest, type = c("Ljung-Box"))
Box.test(error_arima3.0.1, lag = 20, type = c("Ljung-Box"))
Box.test(error_arima3.0.1, lag = 30, type = c("Ljung-Box"))
# Las pruebas nos arrojan el mismo resultado que la Box-Pierce: no hay
# autocorrelación serial en los errores, por lo que se cumple el supuesto de
# ruido blanco en los residuales.

# 3.3.2. Heterocedasticidad
# Revisaremos el supuesto de heterocedasticidad con un test ARCH. Para que
# la hipótesis nula H_0 no se rechaze, la varianza de los errores no debe
# cambiar en el tiempo
x11()
arch_3.0.1 <- arch.test(arima3.0.1, output = TRUE)
ggsave(paste0("ARCH test ", ticker, ".png"))

# Realizamos el test ARCH y parece que la varianza de los residuales es
# constante, por lo que no rechazamos la hipótesis nula H_0 y los residuales
# no son heterocedásticos. Para tener mayor certeza de esto, graficamos la
# FAC y FACP de los residuales al cuadrado

# FAC y FACP de los residuales al cuadrado

lags <- 24
x11()
grid.arrange(
  ggAcf(errores2_arima3.0.1, lag.max = lags, plot = T, lwd = 2) +
      ggtitle("FAC residuales al cuadrado"),

  ggPacf(errores2_arima3.0.1, lag.max = lags,plot = T, lwd = 2) +
      ggtitle("FACP residuales al cuadrado"),

  nrow = 1, ncol = 2
)


# Ahora realizamos una prueba de multiplicadores de Lagrange, que nos
# arroja un valor exacto, que no está sujeto a las interpretaciones de los
# gráficos por parte de quien desarrolla el modelo
n <- length(error_arima3.0.1) # Creamos una variable temporal
time <- 1:n
bptest(error_arima3.0.1 ~ time) 
# p-value = 0.4075

# 3.3.3. Normalidad
jarque.bera.test(error_arima3.0.1)
# p-value<2.2e-16, por lo que rechazamos el supuesto de normalidad de los
# residuales

# A pesar de que los residuales no cumplan el supuesto de normalidad, y
# dado que el supuesto más importante sobre el error se cumplió (no
# correlación serial), elegiremos este modelo, el ARMA(3,1) para realizar
# el pronóstico.

# Ahora realizaremos la prueba de raíz unitaria para verificar invertibilidad
# del polinomio característico de MA(q) y estacionariedad del polinomio
# característico de AR(p) 
fable_3.0.1  <- as_tsibble(diff.ts) %>%
  model(arima110 = ARIMA(value ~ pdq(3,0,1) + PDQ(0, 0, 0)))

x11()
gg_arma(fable_3.0.1)
ggsave(paste0('Prueba raíz unitaria ', ticker, '.png'))

# La función gg_arma(), por defecto, calcula el inverso de las raíces del
# polinomio característico, por lo que, en este caso únicamente, la serie de
# los errores será estacionaria si las raíces se encuentran dentro del
# polinomio.
# El inverso de las raíces se encuentran dentro del polinomio, por lo que
# la parte MA(1) es invertible y AR(3) es estacionaria



# - 3.4. Pronóstico -----------------------------------------------------------
# Usaremos el modelo identificado y validado en los pasos anteriores para
# realizar predicciones de la serie de tiempo.
# Realizamos 10000 simulaciones (lo mínimo recomendable) para hacer la
# predicción 10 pasos adelante
forescast_3.0.1 <- fable_3.0.1 %>% 
  forecast(h = 10, bootstrap = TRUE, times = 10000)

forescast_3.0.1

# Añadimos los intervalos de confianza al 80%, 90% y 95%
forecast_intpred_3.0.1 <- forescast_3.0.1 %>% 
  hilo(level = c(80, 90, 95))

forecast_intpred_3.0.1
print(forecast_intpred_3.0.1)

# Graficamos los pronósticos 10 pasos adelante
# Convertimos a tsible nuestra serie estacionaria:
diff.ts_tsibble <- as_tsibble(diff.ts)
x11()

# Graficamos
forescast_3.0.1 %>%
  autoplot(diff.ts_tsibble) +
    ggtitle("Pronósticos 10 pasos adelante Indice de Producción Industrial ") +
    ylab("Indice de Producción Industrial") + xlab("trimestres") +
    theme_light()
ggsave(paste0("Pronóstico ", ticker, ".png"))


# 3.4.1. Ajuste de pronóstico
# Es importante evaluar la precisión de los pronósticos mediante técnicas como
# la validación cruzada o la comparación de los errores de pronóstico con los
# errores de ajuste en los datos históricos. Vamos a usar la serie estacionaria
# y a esta le vamos a restar los residuales de nuestro modelo ARMA(3,0,1).
ajuste_forecast <- diff.ts - residuals(arima3.0.1)

# Graficamos nuestra serie estacionaria y el ajuste de pronóstico
x11()
plot.ts(diff.ts, type = "l",
    main = "Diferencia ajustada VS diferencia observada", lwd = 1)
points(ajuste_forecast, col = "#FF4040", lwd = 1, type = "l")
legend("topleft", c("observada", "estimada"), col = c("black", "#FF4040"),
    lty = 1, lwd = 2)
ggsave(paste0("Ajuste de pronóstico ", ticker, ".png"))
