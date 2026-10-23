##################################################
############## POSTER ECONOMETRIA ################
##################################################

##################################################
rm(list = ls())   #Limpiar entorno#

##################################################
###################LIBRERIAS######################
##################################################

#install.packages(c("tidyverse", "fable", "tsibble", "feasts", "fs", "lubridate","patchwork","ggthemes"))
library(tidyverse) # Manejo de datos y gráficos (ggplot2)
library(fable)     # Modelamiento SARIMA moderno 
library(tsibble)   # Estructura de series de tiempo para fable
library(feasts)    # Gráficos ACF/PACF y diagnóstico
library(fs)        # Manejo de rutas relativas 
library(lubridate) # Manejo de fechas
library(patchwork) # Para unir varios gráficos en una sola imagen
library(ggthemes)  # Temas estéticos 

# ========================================================= #
####### RUTAS RELATIVAS #######
ruta_archivo <- fs::path("datos.csv") 
mis_datos <- read_csv2(ruta_archivo)

####### TRANSFORMACIÓN SERIE DE TIEMPO #######

datos_ts <- mis_datos %>%
  # Convertimos la columna 'mes' de texto a formato fecha (Día/Mes/Año)
  mutate(mes = yearmonth(dmy(mes))) %>% 
  as_tsibble(index = mes)

####### GRÁFICO INICIAL #######

grafico_serie <- datos_ts %>%
  autoplot(pasajeros) +
  labs(title = "Tráfico Aéreo Nacional de Pasajeros en Colombia",
       x = "Fecha",
       y = "Número de Pasajeros") +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma)

print(grafico_serie)

# ========================================================= #
####### IDENTIFICACIÓN #######

####### Estacionalidad y Tendencia #######
# para establecerlo visualmente utilizaremos una Descomposición STL # 

grafico_descomposicion <- datos_ts %>%
  model(STL(pasajeros ~ trend(window = 21) + season(window = "periodic"))) %>%
  components() %>%
  autoplot() +
  labs(title = "Descomposición STL del Tráfico de Pasajeros",
       subtitle = "Evidencia visual de Tendencia y Estacionalidad") +
  theme_minimal()

print(grafico_descomposicion)

# Correlogramas de la serie original 

grafico_acf_original <- datos_ts %>%
  gg_tsdisplay(pasajeros, plot_type = "partial") +
  labs(title = "Correlograma (ACF y PACF) - Serie en Niveles",
       subtitle = "El decaimiento lento del ACF justifica la presencia de tendencia")

print(grafico_acf_original) 

####### PRUEBAS FORMALES #######
# 1. Pruebas automáticas de número de diferencias (KPSS / OCSB)

# Prueba de si requiere diferencias estacionales (D)
prueba_estacional <- datos_ts %>%
  features(pasajeros, unitroot_nsdiffs) 

print(prueba_estacional) # Resultado: 0 (D = 0)
# Esto nos indica que no sera necesario diferenciar estacionalidad aca
# pues sera suficiente con el modelo SARIMAX para controlar esta estacionalidad

# Prueba de si requiere diferencias ordinarias (d) 

prueba_ordinaria <- datos_ts %>%
  features(pasajeros, unitroot_ndiffs)

print(prueba_ordinaria)  # Resultado: 1 (d = 1)
# Esto nos indica que requerira una diferencia ordinaria 

# 2. Prueba Dickey-Fuller Aumentada (ADF) Tradicional

# En niveles (Serie original)

adf_niveles <- tseries::adf.test(datos_ts$pasajeros)
print(adf_niveles)
#p-value = 0.3355, no se rechaza H0: tiene raiz unitaria (no estacionaria)



# En primera diferencia (d=1)

adf_diferenciada <- tseries::adf.test(diff(datos_ts$pasajeros))
print(adf_diferenciada)
# p-value = 0.01, se rechaza H0 al 99%




#=========================================================================#
####### TRANSFORMACIÓN: PRIMERA DIFERENCIA ORDINARIA (d = 1, D = 0) #######

# Correlograma con lag_max = 36 para observar 3 años de estacionalidad
grafico_acf_dif <- datos_ts %>%
  gg_tsdisplay(difference(pasajeros, 1), 
               plot_type = "partial", 
               lag_max = 36) +
  labs(title = "ACF y PACF de la serie en Primera Diferencia (d=1, D=0)",
       subtitle = "Evaluación de 36 rezagos: evidencia clara de picos estacionales en los meses 12, 24 y 36") +
  theme_minimal()

print(grafico_acf_dif)

# ============================================================ #
####### 3. ESTIMACIÓN Y SELECCIÓN DE MODELOS CANDIDATOS #######

# Estimación simultánea con la variable exógena dummy_covid
modelos_candidatos <- datos_ts %>%
  model(
    # Candidato 1: Autorregresivo ordinario con AR estacional
    princesa_1 = ARIMA(pasajeros ~ dummy_covid + pdq(1,1,0) + PDQ(1,0,0)),
    
    # Candidato 1: Solo términos de Media Móvil
    princesa_2 = ARIMA(pasajeros ~ dummy_covid + pdq(0,1,1) + PDQ(0,0,1)),
    
    # Candidato 3: AR(1) y MA estacional
    princesa_3 = ARIMA(pasajeros ~ dummy_covid + pdq(1,1,0) + PDQ(0,0,1)),
    
    # Candidato 4: Mezcla MA(1) con AR estacional
    princesa_4 = ARIMA(pasajeros ~ dummy_covid + pdq(0,1,1) + PDQ(1,0,0)),
    
    # Candidato 5: Autorregresivo ordinario con mayor peso en el AR estacional
    princesa_5 = ARIMA(pasajeros ~ dummy_covid + pdq(1,1,0) + PDQ(2,0,0)),
  )
    

# Tabla comparativa de Criterios de Información (ordenada de menor a mayor AICc)
tabla_comparacion <- glance(modelos_candidatos) %>%
  select(.model, AIC, AICc, BIC) %>%
  arrange(AICc)

print(tabla_comparacion)
## como lo teniamos previsto, ganaron la princesa #5 y la princesa #1 (autorregresivos)


# ================================================================= #
####### 4. VALIDACIÓN DE SUPUESTOS (DIAGNÓSTICO DE RESIDUOS) #######


# Ajustamos el modelo ganador oficial
modelo_ganador <- datos_ts %>%
  model(
    SARIMAX_Final = ARIMA(pasajeros ~ dummy_covid + pdq(1,1,0) + PDQ(1,0,0))
  )
report(modelo_ganador)
# Diagnóstico visual (Evaluación de Ruido Blanco)
grafico_residuos <- modelo_ganador %>%
  gg_tsresiduals() +
  labs(title = "Diagnóstico de Residuos - Modelo SARIMAX(1,1,0)(1,0,0)[12]",
       subtitle = "Los residuos se comportan como Ruido Blanco (sin autocorrelación significativa)")

print(grafico_residuos)

# Prueba formal de Ljung-Box (No autocorrelación)
# dof = p + q + P + Q = 1 + 0 + 2 + 0 = 3 grados de libertad
prueba_ljung_box <- augment(modelo_ganador) %>%
  features(.innov, ljung_box, lag = 24, dof = 3)

print(prueba_ljung_box)

# ============================================== #
####### 5. PRONÓSTICO #######


# A. Crear datos futuros para la variable exógena (dummy_covid = 0)
datos_futuros <- new_data(datos_ts, n = 12) %>%
  mutate(dummy_covid = 0) # En el futuro no hay choque COVID

# B. Generar el pronóstico a 12 meses
pronostico_final <- modelo_ganador %>%
  forecast(new_data = datos_futuros)

# C. Graficar el pronóstico con intervalos de confianza para el póster
grafico_pronostico <- pronostico_final %>%
  autoplot(datos_ts) +
  labs(title = "Pronóstico del Tráfico Aéreo Nacional de Pasajeros (Próximos 12 meses)",
       subtitle = "Modelo SARIMAX(1,1,0)(1,0,0)[12] con sombra de intervalos al 80% y 95%",
       x = "Fecha",
       y = "Número de Pasajeros") +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma)

print(grafico_pronostico)

################# FIN #################