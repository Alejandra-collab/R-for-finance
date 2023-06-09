# R-for-finance

## Metodología Box-Jenkins en el pronóstico de series temporales usando un modelo ARIMA(p,d,q)
En el análisis de series de tiempo, la predicción futura de los valores de una serie es una tarea crucial en la toma de decisiones en diversas áreas, como la economía, la finanzas, la ingeniería, entre otras. En este trabajo, utilizamos la metodología Box-Jenkins para realizar el pronóstico 10 pasos adelante de la serie "Ingresos personales y gastos personales de Estados Unidos", desestacionalizado y obtenido de la base de datos de la Reserva Federal de St. Louis (FRED).. Esta metodología consta de cuatro pasos: identificación, estimación, validación y pronóstico.

## Metodología de Johansen en el análisis de cointegración de series temporales
La Metodología de Johansen es ampliamente utilizada en el análisis de series de tiempo cointegradas, es decir, series de tiempo que tienen una relación a largo plazo.
Queremos determinar el orden de integración en las series y luego usar los criterios de información de un Modelo Vectorial Autorregresivo (VAR) para conocer el número de rezagos del Modelo de Corrección de Errores Vectorial (VECM).
Después, determinaremos el número de relaciones de cointegración y estimaremos el VECM para analizar el vector de cointegración, que representa la relación a largo plazo entre las variables (podremos ver cómo se combinan las variables). De igual manera, podremos ver cómo cada variable se va ajustando hasta el equilibrio de largo plazo.
Todo lo anterior permitirá que podamos hacer uso del modelo. Con el modelo obtenido podremos realizar análisis de impulso-respuesta (efectos en la variable i por choques de la variable j e i, a lo largo del tiempo), en donde, además de analizar los efectos por choques de las variables, veremos la dirección de la causalidad entre las variables, si las hay.