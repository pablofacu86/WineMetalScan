# Scan FT-MIR — app de predicción (Posit Connect Cloud)

App Shiny para estimar potasio, magnesio y calcio (regresión) y clasificar
hierro y cobre (por encima/debajo de su límite de interés) en vino a partir
de un espectro FT-MIR, usando los modelos validados en el estudio.

## Estructura esperada

```
app.R
R/predict_utils.R
Models/            <- generado por Build_Deployment_Models.R (POTASSIUM_model.rds, etc.)
data/
  Test_Set_App.xlsx
www/
  ivagro_logo.png             <- logo de IVAGRO
  uca_logo.png                <- logo de la Universidad de Cadiz
```

**Antes de desplegar**, correr `Build_Deployment_Models.R` (en el repo del
código; ya actualizado para calcular tambien los datos del dendrograma en
hierro y cobre) y copiar los 5 archivos `.rds` resultantes a `Models/`. Sin
esos archivos, cada solapa igual carga y muestra la información de alcance,
pero "Ver estado de validación" y la predicción muestran un aviso de que el
modelo todavía no está disponible, en vez de romper la app.

**Importante:** los `.rds` que ya tenías generados (antes de este cambio)
no incluyen los datos del dendrograma. La app no se rompe con ellos —
simplemente muestra "Dendrograma no disponible para este modelo" en las
solapas de hierro y cobre hasta que vuelvas a correr
`Build_Deployment_Models.R` y reemplaces esos dos `.rds`.

## Qué falta completar (contenido, no código)

- Todo lo que faltaba en la versión anterior (logo de IVAGRO, los 5 `.rds`)
  ya está cargado. Los informes PDF completos se sacaron de la app
  deliberadamente (pesaban demasiado para el repositorio); si en algún
  momento se quieren volver a ofrecer, conviene alojarlos aparte (por
  ejemplo como asset de un GitHub Release) y linkearlos desde afuera de la
  app, no volver a incluirlos como archivo local.
- Si en algún momento se quiere agregar una sección con los dendrogramas u
  otros gráficos de diagnóstico de cada analito, conviene armarla como una
  solapa propia ("Material adicional" o similar) en vez de mezclarla dentro
  de la carpeta `www/` (esa carpeta es solo para los assets de la interfaz,
  como los logos) o dentro del panel de "Estado de validación" de cada
  analito.

## Paquetes de R necesarios

`shiny`, `bslib` (≥ 0.5, por `page_navbar`/`nav_panel`/`layout_columns`),
`DT`, `ggplot2`, `openxlsx`, `prospectr`, `caret`, `xgboost`, `dendextend`,
`circlize`, y `kernlab`
instalado (lo usa `caret` internamente para los modelos SVM, sin necesidad
de cargarlo con `library()`).

## Cómo funciona la predicción sobre un espectro nuevo

1. El usuario sube un archivo (.csv o .xlsx) con número de onda y
   absorbancia, o un archivo con el mismo formato de fila que los datasets
   del estudio (ID, valor de referencia, 545 variables).
2. El espectro se **realinea por interpolación lineal** a la grilla de
   referencia del modelo (902.57–3000.84 cm⁻¹, ~3.86 cm⁻¹ de paso). Esto
   es lo que permite aceptar espectros de otro equipo, con otro paso de
   muestreo o rango levemente distinto — siempre que cubra razonablemente
   la zona 900–3000 cm⁻¹; si la cobertura es menor al 90 % se avisa que la
   predicción puede no ser confiable.
3. Se aplica exactamente el mismo preprocesamiento con el que se entrenó
   el modelo ganador de ese analito (derivada, SNV o MSC con la referencia
   guardada del entrenamiento, recorte a las ventanas de Patz, centrado,
   variables de Boruta si corresponde).
4. Se predice con el modelo cargado. Para regresión, se informa además un
   intervalo aproximado del 95 % como valor ± 2×RMSE del test (una
   aproximación simple, no un intervalo de predicción estadísticamente
   riguroso — aclarado así en el informe descargable). Para clasificación,
   se muestra la clase predicha y la probabilidad de cada clase.

## Limitaciones a tener en cuenta

- El intervalo de "± 2×RMSE" asume que el error de predicción es
  aproximadamente normal y constante en todo el rango de concentración,
  lo cual es una simplificación razonable pero no una garantía estadística
  formal.
- La re-alineación por interpolación ayuda con diferencias de grilla, pero
  no corrige diferencias más profundas entre instrumentos (línea de base,
  resolución óptica, relación señal/ruido). Se recomienda, cuando sea
  posible, validar con muestras de concentración conocida del propio
  equipo antes de confiar en las predicciones de rutina (para eso están
  las muestras de `data/Test_Set_App.xlsx`, con valores de referencia
  conocidos).
