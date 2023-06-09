# __________________________________________________________________________
# ________ METODOLOGÍA DE JOHANSEN EN EL ANÁLISIS DE COINTEGRACIÓN  ________
# __________________________ DE SERIES TEMPORALES __________________________
# __________________________________________________________________________
#
#
# 1. INTRODUCCIÓN
# La Metodología de Johansen es ampliamente utilizada en el análisis de series
# de tiempo cointegradas, es decir, series de tiempo que tienen una relación
# a largo plazo.
# Queremos determinar el orden de integración en las series y luego usar los
# criterios de información de un Modelo Vectorial Autorregresivo (VAR) para
# conocer el número de rezagos del Modelo de Corrección de Errores Vectorial
# (VECM).
# Después, determinaremos el número de relaciones de cointegración y estimaremos
# el VECM para analizar el vector de cointegración, que representa la relación
# a largo plazo entre las variables (podremos ver cómo se combinan las
# variables). De igual manera, podremos ver cómo cada variable se va ajustando
# hasta el equilibrio de largo plazo.
# Todo lo anterior permitirá que podamos hacer uso del modelo. Con el modelo
# obtenido podremos realizar análisis de impulso-respuesta (efectos en la
# variable i por choques de la variable j e i, a lo largo del tiempo), en
# donde, además de analizar los efectos por choques de las variables, veremos
# la dirección de la causalidad entre las variables, si las hay.
#
#
# --------- 2. SELECCIÓN DE LAS SERIES Y PREPARACIÓN DE LOS DATOS -----------
#
# 2.1. Selección de las series
# De acuerdo a la teoría cuantitativa del dinero, existe una relación directa
# entre la cantidad de dinero en circulación (oferta monetaria) y el nivel
# general de precios. Irving Fisher, quien propuso esta relación, la expresó
# matemáticamente de la forma MV=PQ, en donde M es cantidad de dinero, V
# velocidad del dinero, P nivel general de precios, y Q nivel de producción.
# Aplicaremos la Metodología de Johansen para analizar la relación de
# cointegración entre la cantidad de dinero en la economía (M2) y en nivel
# general de precios (IPC) para la economía estadounidense entre enero de 2000
# a enero de 2020.
# Si la teoría cuantitativa se cumple para este caso, veremos que ante aumentos
# de la cantidad de dinero, manteniendo las demás variables constantes, el
# nivel general de precio aumentará, y del mismo modo ante cambios negativos.
#
#
# 2.2. Preparación de los datos
#
# Limpiamos las variables de entorno
rm(list = ls())

# Librerías
packages <- c("vars", "urca", "ggplot2", "ggfortify", "gridExtra", "dplyr",
            "tidyr", "readxl", "tsDyn", "VAR.etp", "lubridate")
install_packages <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    }
  }
}

install_packages(packages)

library_packages <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      message("El paquete '", package, "' no está instalado. 
      Intentar con install.packages().")
    } else {
      library(package, character.only = TRUE)
    }
  }
}

library_packages(packages)

install.packages("urca")
library(urca)

# Cargar las series
path_file <- file.choose()
data <- read_excel(path_file)

# Adjuntar datos al entorno de trabajo
attach(data)

# Convertir datos a ts
money_supply <- ts(M2SL, start = c(2000, 1), frequency = 12)
price_index <- ts(CPIAUCSL, start = c(2000, 1), frequency = 12)


# ----- 3. PASO 1: DETERMINACIÓN DEL GRADO DE INTEGRACIÓN Y REZAGOS ------
#
# 3.1. Grado de integración
# Vamos a determinar el grado de integración necesario para que las series de
# tiempo sean estacionarias a largo plazo. Antes de hacer esto, veamos cómo
# se comportan las series en niveles (originales, sin diferenciación).

start_year <- min(year(data$DATE))
end_year <- max(year(data$DATE))

x11()
grid.arrange(
  ggplot(data = money_supply, aes(x = DATE, y = M2SL)) +
    geom_line(color = "#FF4040", size = 1) +
    labs(title = paste0("Money supply EE.UU. ", "(", start_year,
    " - ", end_year, ")"), x = "Fecha", y = "Billones de dólares"),

  ggplot(data = price_index, aes(x = DATE, y = CPIAUCSL)) +
    geom_line(color = "deepskyblue4", size = 1) +
    labs(title = paste0("IPC EE.UU. ", "(", start_year, " - ", end_year, ")"),
    x = "Fecha", y = "Índice"),
  nrow = 1, ncol = 2
)

# Es evidente que las series no son estacionarias, pero vamos a realizar la
# prueba de raíz unitaria de Dickey-Fuller aumentado (ADF). Usaremos el
# criterio de información AIC. No usamos el BIC, aunque sea más parsimonioso,
# porque el AIC incluye más rezagos de las series diferenciadas. Si no se
# incluyen los suficientes rezagos, no se puede controlar la correlación en
# el error y la prueba no sería válida.
# La prueba cambia si incluimos o no constante, tendencia y la combinación de
# ambas, Así que realizamos varias pruebas.
# Las pruebas para analizar la significancia de la tendencia y de la deriva
# son pruebas de cola derecha, es decir, si el estadístico t se encuentra a
# derecha del valor crítico rechazamos la hipótesis nula de significacnia 0.

# Primero, las series en niveles, con y sin constante.
# Money supply con trend:
summary(ur.df(money_supply, lags = 6, selectlags = "AIC", type = "trend"))
# Money supply con constante:
summary(ur.df(money_supply, lags = 6, selectlags = "AIC", type = "drift"))
# Money supply sin constante y sin tenencia:
summary(ur.df(money_supply, lags = 6, selectlags = "AIC", type = "none"))

# En este caso, la tendencia resultó ser no significativa, mientras que la
# deriva sí lo es. En esa misma prueba de deriva verificamos estacionariedad
# con tau2, que nos dice que la serie en niveles no es estacionaria.

# IPC con trend
summary(ur.df(price_index, lags = 6, selectlags = "AIC", type = "trend"))
# IPC con constante:
summary(ur.df(price_index, lags = 6, selectlags = "AIC", type = "drift"))
# IPC supply sin constante y sin tenencia:
summary(ur.df(price_index, lags = 6, selectlags = "AIC", type = "none"))

# Del mismo modo que con la serie money_supply, la tendenica no es
# significativa y la deriva sí lo es. Además, en la prueba de significancia
# para la deriva, la serie resultó ser no estacionaria, de acuerdo al t
# estadístico para verificar estacionariedad, tau2.

# Ahora que sabemos que las series no son estacionarias, vamos a realizar la
# misma prueba para las series diferenciadas, incluyendo el parámetro
# de la deriva que es significativo.

summary(ur.df(diff(money_supply), lags = 6, selectlags = "AIC", type = "drift"))
summary(ur.df(diff(price_index), lags = 6, selectlags = "AIC", type = "drift"))

# Las series son integradas de orden 1, es decir que al realizarles una
# diferencia son estacionarias, requisito para realizar la prueba de
# cointegración de Johansen.

# 3.2. Determinamos el número de rezagos del VECM a través de un VAR en niveles

# Es importante la manera en la que tratamos a las variables (si son
# dependientes o independientes).
# Siguiendo las ideas de la Teoría Cuantitativa del dinero, vamos a tratar al
# nivel de precios como la variable explicada y a la oferta monetaria como
# la variable explicativa.

series_matrix <- cbind(price_index, money_supply)

# Considerando rezagos en las variables exógena y endógena:
VARselect(series_matrix, lag.max = 6, type = "both", season = NULL)
# De acuerdo a las pruebas HQ (Hannan-Quinn Criterion) y SC (Schwarz Criterion),
# el modelo más parsimonioso es un VAR(2)

# Considerando rezagos en la variable endógena (price_index):
VARselect(series_matrix, lag.max = 6, type = "const", season = NULL)
# De acuerdo a la prueba SC (Schwarz Criterion), el modelo más parsimonioso
# es un VAR(2)

# Sin considerar rezagos en ninguna variable:
VARselect(series_matrix, lag.max = 6, type = "none", season = NULL)
# De acuerdo a la prueba SC (Schwarz Criterion), el modelo más parsimonioso
# es un VAR(2)

# Ahora, vamos a ajustar el modelo vectorial autorregresivo (VAR) al
# conjunto de datos, con los rezagos que determinamos adecuados
# anteriormente:
summary(VAR(series_matrix, p = 2, type = "both", season = NULL))
# Aunque para un modelo la tendencia es significativa, para el otro no lo es,
# así que no la añadiremos al modelo VAR

# Ahora veamos la constante:
summary(VAR(series_matrix, p = 2, type = "const", season = NULL))
# Para las dos ecuaciones, la deriva no es significativa

# Entonces, nuestro VAR va a estar determinado por:
selected_var <- VAR(series_matrix, p = 2, type = "const", season = NULL)

# Analizaremos el comportamiento de los residuales; si no se comportan
# adecuadamente (hay autocorrelación serial), agregaremos más rezagos.
# Buscamos no rechazar la hipótesis de ausencia de serialidad en los
# residuos, la cual es un requisito para aplicar la metodología de Johansen.
residuals_75 <- serial.test(selected_var, lags.pt = 70, type = "PT.asymptotic")
residuals_75

residuals_30 <- serial.test(selected_var, lags.pt = 30, type = "PT.asymptotic")
residuals_30

residuals_20 <- serial.test(selected_var, lags.pt = 20, type = "PT.asymptotic")
residuals_20

residuals_10 <- serial.test(selected_var, lags.pt = 10, type = "PT.asymptotic")
residuals_10

# Para las cuatro opciones no rechazamos la hipótesis nula de no autocorrelación
# serial en los errores, así que continuamos con el proceso.

# Graficamos los residuales:
x11()
plot(residuals_20, names = "price_index")
plot(residuals_20, names = "money_supply")

# ------------- 4. PASO 2: DETERMINANDO EL RANGO DE LA MATRIZ PI ---------------

# 4.1. Johansen con constante fuera del vector de cointegración
# 4.1.1. Criterio del valor propio máximo

max_value_criterion_out <- ca.jo(series_matrix, ecdet = "none", type = "eigen",
                       K = 2, spec = "transitory", season = NULL)
summary(max_value_criterion_out)
# El resultado es un t estadístico de 17.41, contra un valor crítico de 14.90
# al 5%. El t estadístico cae en zona de rechazo: rechazamos que las relaciones
# de cointegración sean 0 y pasamos a la siguiente hipótesis nula: 1 relación
# de cointegración. Para este caso, el t estadístico es 3.27 y el valor crítico
# a 5% es 8.18. El estadístico cae en la zona de no rechazo, así que existe una
# relación de cointegración.
# Aunque estos valores son confiables, haremos uso del Criterio de la traza
# para estar seguros.

# H_0: r=0. Cola derecha

# 4.1.2. Criterio de la traza:
trace_criterion_out <- ca.jo(series_matrix, ecdet = "none", type = "trace",
                        K = 2, spec = "transitory", season = NULL)
summary(trace_criterion_out)
# Confirmamos una relación de cointegración.

# 4.2. Johansen con constante dentro del vector de cointegración
# 4.2.1. Criterio del valor propio
max_value_criterion_in <- ca.jo(series_matrix, ecdet = "const", type = "eigen",
                        K = 2, spec = "longrun", season = NULL)
summary(max_value_criterion_in)
# La prueba nos muestra que se rechaza la hipótesis nula de relaciones de
# cointegración = 0, por lo que hay >=1 relaciones de cointegración.

# 4.2.2. Criterio de la traza
trace_criterion_in <- ca.jo(series_matrix, ecdet = "const", type = "trace",
                      K = 2, spec = "longrun", season = NULL)
summary(trace_criterion_in)
# Confirmamos que hay al menos una relación de cointegración.

# 4.3. Determinando la tendencia lineal del modelo VEC
# Veremos si debemos incluir el modelo con constante en la relación de
# cointegración
lttest(max_value_criterion_out, r = 1)
lttest(max_value_criterion_in, r = 1)
# Tenemos resultados consistentes. La prueba sugiere no incluir el modelo
# con constante en la relación de cointegración.

# 4.4. Estimamos el modelo VEC con la constante fuera del vector de
# cointegración
vec_drift_out <- cajorls(max_value_criterion_out, r = 1)
vec_drift_out

# ------------- 5. PASO 3: ANÁLISIS DE LA MATRIZ BETA Y ALPHA ---------------

# El vector de cointegración normalizado:
coefB(vec_drift_out2)

# Coeficientes de velocidad de ajuste:
coefA(vec_drift_out2)

# -------------------- 5. PASO 4: VALIDANDO LOS SUPUESTOS -------------------

# Reparametrización de VEC como VAR en niveles
vec_to_var <- vec2var(max_value_criterion_out, r = 1)
vec_to_var

# 5.1. Validación de supuestos
# 5.1.1. Autocorrelación
autocorr75 <- serial.test(vec_to_var, lags.pt = 75, type = "PT.asymptotic")
autocorr75

autocorr30 <- serial.test(vec_to_var, lags.pt = 30, type = "PT.asymptotic")
autocorr30

autocorr20 <- serial.test(vec_to_var, lags.pt = 20, type = "PT.asymptotic")
autocorr20

autocorr10 <- serial.test(vec_to_var, lags.pt = 10, type = "PT.asymptotic")
autocorr10

# Para ningún nivel de retardos utilizados en la prueba se puede rechazar la
# hipótesis nula de no autocorrelación.

# Graficamos los residuales
x11()
plot(autocorr20, names = "price_index")
x11()
plot(autocorr20, names = "money_supply")

# 5.1.2. Homocedasticidad

arch.test(vec_to_var, lags.multi = 24, multivariate.only = TRUE)
arch.test(vec_to_var, lags.multi = 12, multivariate.only = TRUE)
# El valor p (0.1333) es mayor que el umbral de significancia comúnmente
# utilizado de 0.05. Esto sugiere que no hay suficiente evidencia para
# rechazar la hipótesis nula de homocedasticidad condicional en los residuos
# del modelo VAR. Pero al realizar la misma prueba para 12 rezagos, el valor
# p (1.119e-06) es extremadamente pequeño, lo que indica una evidencia fuerte
# en contra de la hipótesis nula de homocedasticidad condicional en los
# residuos del modelo VAR.


# 5.1.3. Normalidad
normality.test(vec_to_var)
# En este caso, el valor p (< 2.2e-16) es extremadamente pequeño, lo que indica
# una evidencia fuerte en contra de la hipótesis nula de normalidad en los
# residuos del modelo VAR.

# 5.2. Pronóstico
x11()
predict(vec_to_var, n.ahead = 12)
plot(predict(vec_to_var, n.ahead = 12))

# 5.2.1. Función impulso-respuesta (FIR)
impulso_respuesta <- function(var, impulso, respuesta, pasos_adelante,
                              ortog, int_conf, titulo){
  # Cálculo de la FIR
  total_pasos_futuros <- length(pasos_adelante) - 1
  FIR <- fir(var, impulse = impulso, response = respuesta,
             n.ahead = total_pasos_futuros, ortho = ortog, ci = int_conf)
  FIR_data_frame <- data.frame(FIR$fir, FIR$Lower, FIR$Upper, pasos_adelante)

  # Gráfica
  graph <- ggplot(data = FIR_data_frame,
                  aes(x = pasos_adelante, y = FIR_data_frame[, 1],
                      ymin = FIR_data_frame[, 2], ymax = FIR_data_frame[, 3])) +
    geom_hline(yintercept = 0, color = "red") +
    geom_ribbon(fill = "grey", alpha = 0.2) +
    geom_line() +
    theme_light() +
    ggtitle(titulo) +
    ylab("") +
    xlab("pasos adelante") +
    theme(plot.title = element_text(size = 11, hjust = 0.5),
          axis.title.y = element_text(size = 11))

  return(graph)
}

# Definimos el número de pasos adelante
pasos_adelante <- 0:18

# FIR de las variables del sistema ante distintos choques exógenos.
v1_v1 <- impulso_respuesta(vec_to_var, "price_index", "price_index",
                           pasos_adelante = pasos_adelante,
                           ortog = TRUE, int_conf = 0.95,
                           titulo =
                           "Impulso de price_index - respuesta de price_index")
v1_v2 <- impulso_respuesta(vec_to_var, "price_index", "money_supply",
                           pasos_adelante = pasos_adelante,
                           ortog = TRUE, int_conf = 0.95,
                           titulo =
                           "Impulso de price_index - respuesta de money_supply")
v2_v1 <- impulso_respuesta(vec_to_var, "money_supply", "price_index",
                           pasos_adelante = pasos_adelante,
                           ortog = TRUE, int_conf = 0.95,
                           titulo =
                           "Impulso de money_supply - respuesta de price_index")
v2_v2 <- impulso_respuesta(vec_to_var, "money_supply", "money_supply",
                           pasos_adelante = pasos_adelante,
                           ortog = TRUE, int_conf = 0.95,
                           titulo =
                           "Impulso de money_supply - respuesta de money_supply")

# Mostrar las gráficas en una matriz de 2x2
grid.arrange(v1_v1, v1_v2, v2_v1, v2_v2, ncol = 2)