# ==========================================================================
# Scan FT-MIR -- estimacion de metales en vino por espectroscopia FT-MIR
# App para Posit Connect Cloud
#
# Estructura:
#   - Una solapa por analito (potasio, magnesio, calcio, hierro, cobre).
#   - Cada solapa: informacion de alcance/instrumento/protocolo EDTA,
#     estado de validacion (metricas + graficos del modelo ganador),
#     y un panel para subir un espectro propio y obtener una prediccion.
#   - Los modelos se cargan en tiempo de ejecucion desde Models/*.rds
#     (generados por Build_Deployment_Models.R). Si un .rds todavia no
#     existe, la solapa lo indica en vez de fallar.
# ==========================================================================

library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(openxlsx)

source("R/predict_utils.R")

MODELS_DIR  <- "Models"
DATA_DIR    <- "data"

# --------------------------------------------------------------------
# Informacion comun a los cinco analitos (instrumento, metodo de
# referencia, protocolo de preparacion con EDTA, alcance varietal).
# --------------------------------------------------------------------
COMMON_INFO <- list(
  instrument = paste(
    "Espectros adquiridos por duplicado con un sistema MultiSpec (TDI, Barcelona,",
    "Espana), equipado con una celda de flujo de transmision CETIM (TDI, Barcelona,",
    "Espana) y un automuestreador de 40 posiciones. Rango espectral: 902.57-3000.84",
    "cm-1 (545 variables, paso ~3.86 cm-1)."),
  reference_method = paste(
    "Las concentraciones de referencia se determinaron por espectrometria de",
    "emision optica con plasma acoplado inductivamente (ICP-OES), equipo",
    "Spectrogreen FMD46 (SPECTRO Analytical Instruments)."),
  edta_protocol = paste(
    "Preparacion de la muestra: agregar 2 mL de una solucion acuosa de Na2EDTA al",
    "5% p/v a 8 mL de vino (volumen final 10 mL). Homogeneizar suavemente y dejar",
    "en reposo 30 minutos a temperatura ambiente antes de medir. El agregado de",
    "Na2EDTA es necesario porque los metales no absorben directamente en el",
    "infrarrojo; la formacion del complejo metal-EDTA es lo que produce la senal",
    "espectral utilizada por el modelo."),
  wine_scope = paste(
    "Vinos españoles blancos y tintos de 20 Denominaciones de Origen (Cadiz, Borja,",
    "Cariñena, Castilla, Catalunya, Extremadura, Jumilla, La Jaraba, La Mancha,",
    "Madrid, Navarra, Penedes, Ribeiro, Ribera del Duero, Rioja, Rueda, Somontano,",
    "Toro, Ucles y Valdepeñas) y 9 variedades de uva (Airen, Garnacha, Merlot,",
    "Monastrell, Palomino, Sauvignon Blanc, Tempranillo, Verdejo y Viura)."),
  disclaimer = paste(
    "Los valores que devuelve esta herramienta son ESTIMACIONES obtenidas mediante",
    "algoritmos de aprendizaje automatico a partir del espectro FT-MIR, y deben",
    "tomarse como valor ORIENTATIVO o de cribado (screening). No reemplazan al",
    "analisis por metodos de referencia (ICP-OES/ICP-MS u otros metodos oficiales",
    "OIV) y no deben utilizarse como unico criterio para decisiones enologicas,",
    "regulatorias o comerciales.")
)

IVAGRO_URL <- "https://ivagro.uca.es"

# --------------------------------------------------------------------
# Configuracion especifica de cada analito: alcance, modelo ganador,
# metricas de validacion y rutas a los archivos asociados. Los valores
# numericos son los reportados en el estudio; se actualizan solos en
# cuanto el bundle .rds correspondiente este disponible (ver mas abajo,
# load_bundle()).
# --------------------------------------------------------------------
ANALYTES <- list(
  POTASSIUM = list(
    label = "Potasio", symbol = "K", unit = "mg/L", task = "regression",
    color = "#7B1E3A",
    scope_range = "290 - 1379 mg/L", n_samples = "100 (46 blancos, 54 tintos)",
    model_desc = "Elastic Net, sin derivada ni correccion de dispersion, sin seleccion de variables (310 variables).",
    spiked = FALSE,
    model_file = file.path(MODELS_DIR, "POTASSIUM_model.rds"),
    sample_file = file.path(DATA_DIR, "Test_Set_App.xlsx"), sample_sheet = "POTASSIUM"
  ),
  MAGNESIUM = list(
    label = "Magnesio", symbol = "Mg", unit = "mg/L", task = "regression",
    color = "#3A5A40",
    scope_range = "47.9 - 118 mg/L", n_samples = "100 (46 blancos, 54 tintos)",
    model_desc = "PLS (5 variables latentes), 1ra derivada + MSC, seleccion Boruta (45 variables).",
    spiked = FALSE,
    model_file = file.path(MODELS_DIR, "MAGNESIUM_model.rds"),
    sample_file = file.path(DATA_DIR, "Test_Set_App.xlsx"), sample_sheet = "MAGNESIUM"
  ),
  CALCIUM = list(
    label = "Calcio", symbol = "Ca", unit = "mg/L", task = "regression",
    color = "#1D3557",
    scope_range = "45.0 - 295 mg/L (rango ampliado mediante fortificacion, ver nota abajo)",
    n_samples = "160 (79 blancos, 81 tintos; incluye muestras fortificadas)",
    model_desc = "XGBoost, 1ra derivada + SNV, seleccion Boruta (47 variables).",
    spiked = TRUE,
    model_file = file.path(MODELS_DIR, "CALCIUM_model.rds"),
    sample_file = file.path(DATA_DIR, "Test_Set_App.xlsx"), sample_sheet = "CALCIUM"
  ),
  IRON = list(
    label = "Hierro", symbol = "Fe", unit = "mg/L", task = "classification",
    color = "#BC6C25",
    scope_range = "Clasificacion binaria respecto al limite de 10 mg/L",
    n_samples = "160 (69 blancos, 91 tintos; incluye muestras fortificadas)",
    model_desc = "SVM (nucleo radial), sin preprocesamiento adicional, seleccion Boruta (58 variables).",
    spiked = TRUE, limit = 10,
    model_file = file.path(MODELS_DIR, "IRON_model.rds"),
    sample_file = file.path(DATA_DIR, "Test_Set_App.xlsx"), sample_sheet = "IRON"
  ),
  COPPER = list(
    label = "Cobre", symbol = "Cu", unit = "mg/L", task = "classification",
    color = "#606C38",
    scope_range = "Clasificacion binaria respecto al limite de 1 mg/L",
    n_samples = "159 (70 blancos, 89 tintos; incluye muestras fortificadas)",
    model_desc = "SVM (nucleo radial), correccion SNV, seleccion Boruta (77 variables).",
    spiked = TRUE, limit = 1,
    model_file = file.path(MODELS_DIR, "COPPER_model.rds"),
    sample_file = file.path(DATA_DIR, "Test_Set_App.xlsx"), sample_sheet = "COPPER"
  )
)

# --------------------------------------------------------------------
# Carga perezosa y cacheada de cada bundle .rds. Devuelve NULL (sin
# fallar la app) si el archivo todavia no fue generado.
# --------------------------------------------------------------------
.bundle_cache <- new.env(parent = emptyenv())
load_bundle <- function(analyte_key) {
  if (!is.null(.bundle_cache[[analyte_key]])) return(.bundle_cache[[analyte_key]])
  path <- ANALYTES[[analyte_key]]$model_file
  if (!file.exists(path)) return(NULL)
  b <- readRDS(path)
  .bundle_cache[[analyte_key]] <- b
  b
}

# --------------------------------------------------------------------
# Reporte descargable (HTML autocontenido) con el detalle del modelo:
# alcance, metodo, preprocesamiento, hiperparametros y metricas.
# --------------------------------------------------------------------
render_model_report_html <- function(cfg, bundle) {
  metric_rows <- if (cfg$task == "regression") {
    m <- bundle$metrics
    sprintf(
      "<tr><td>R² (test)</td><td>%.3f</td></tr><tr><td>RMSE (test)</td><td>%.2f %s</td></tr>
       <tr><td>RPD (test)</td><td>%.2f</td></tr><tr><td>R² (entrenamiento)</td><td>%.3f</td></tr>
       <tr><td>RMSE (entrenamiento)</td><td>%.2f %s</td></tr>",
      m$R2_Test, m$RMSE_Test, cfg$unit, m$RPD_Test, m$R2_Train, m$RMSE_Train, cfg$unit)
  } else {
    m <- bundle$metrics
    sprintf(
      "<tr><td>AUC (test)</td><td>%s</td></tr><tr><td>Accuracy (test)</td><td>%.3f</td></tr>
       <tr><td>Kappa (test)</td><td>%.3f</td></tr><tr><td>Sensibilidad (test)</td><td>%.3f</td></tr>
       <tr><td>Especificidad (test)</td><td>%.3f</td></tr>",
      ifelse(is.na(m$AUC_Test), "N/D", sprintf("%.3f", m$AUC_Test)),
      m$Accuracy_Test, m$Kappa_Test, m$Sensitivity_Test, m$Specificity_Test)
  }
  hp <- paste(sprintf("%s = %s", names(bundle$hyperparams), unlist(bundle$hyperparams)), collapse = "; ")
  sprintf('<html><head><meta charset="utf-8"><style>
    body{font-family:Arial,Helvetica,sans-serif;max-width:800px;margin:2em auto;color:#222}
    h1{color:%s} table{border-collapse:collapse;width:100%%;margin:1em 0}
    td,th{border:1px solid #ccc;padding:6px 10px;text-align:left} th{background:#f2f2f2}
    .disclaimer{background:#fff3cd;border:1px solid #ffe08a;padding:1em;border-radius:6px;margin-top:2em}
  </style></head><body>
  <h1>Scan FT-MIR &mdash; Ficha del modelo: %s (%s)</h1>
  <p><i>Generado: %s</i></p>
  <h2>Alcance</h2>
  <table><tr><th>Rango validado</th><td>%s</td></tr>
  <tr><th>N&deg; de muestras (total)</th><td>%s</td></tr>
  <tr><th>N&deg; entrenamiento / prueba</th><td>%d / %d</td></tr>
  <tr><th>Variedades y origen</th><td>%s</td></tr></table>
  <h2>Instrumento y metodo de referencia</h2>
  <table><tr><th>FT-MIR</th><td>%s</td></tr>
  <tr><th>Metodo de referencia</th><td>%s</td></tr></table>
  <h2>Preparacion de la muestra (Na<sub>2</sub>EDTA)</h2>
  <p>%s</p>
  <h2>Modelo</h2>
  <table><tr><th>Algoritmo y preprocesamiento</th><td>%s</td></tr>
  <tr><th>Hiperparametros</th><td>%s</td></tr></table>
  <h2>Metricas de desempeño</h2>
  <table><tr><th>Metrica</th><th>Valor</th></tr>%s</table>
  <div class="disclaimer"><b>Advertencia:</b> %s</div>
  </body></html>',
  cfg$color, cfg$label, cfg$symbol, format(Sys.time(), "%Y-%m-%d %H:%M"),
  cfg$scope_range, cfg$n_samples, bundle$n_train, bundle$n_test, COMMON_INFO$wine_scope,
  COMMON_INFO$instrument, COMMON_INFO$reference_method, COMMON_INFO$edta_protocol,
  cfg$model_desc, hp, metric_rows, COMMON_INFO$disclaimer)
}

render_prediction_report_html <- function(cfg, result) {
  body <- if (result$task == "regression") {
    sprintf("<p style='font-size:1.4em'><b>%.1f %s</b> (intervalo aprox. 95%%: %.1f - %.1f %s)</p>
             <p>Incertidumbre estimada: &plusmn; 2 &times; RMSE de test (%.1f %s).</p>",
            result$value, cfg$unit, result$lower, result$upper, cfg$unit, result$uncertainty / 2, cfg$unit)
  } else {
    probs <- paste(sprintf("%s: %.1f%%", names(result$probabilities), as.numeric(result$probabilities) * 100), collapse = "<br>")
    sprintf("<p style='font-size:1.4em'><b>%s</b></p><p>Probabilidades por clase:<br>%s</p>", result$class, probs)
  }
  resample_note <- if (isTRUE(result$resampled)) {
    "<p><i>Nota: el espectro subido no coincidia con la grilla de referencia del modelo; se re-muestreo por interpolacion lineal antes de predecir.</i></p>"
  } else ""
  sprintf('<html><head><meta charset="utf-8"><style>
    body{font-family:Arial,Helvetica,sans-serif;max-width:700px;margin:2em auto;color:#222}
    h1{color:%s} .disclaimer{background:#fff3cd;border:1px solid #ffe08a;padding:1em;border-radius:6px;margin-top:2em}
  </style></head><body>
  <h1>Scan FT-MIR &mdash; Resultado de prediccion: %s (%s)</h1>
  <p><i>Generado: %s</i></p>
  %s %s
  <div class="disclaimer"><b>Advertencia:</b> %s</div>
  </body></html>',
  cfg$color, cfg$label, cfg$symbol, format(Sys.time(), "%Y-%m-%d %H:%M"), body, resample_note,
  COMMON_INFO$disclaimer)
}

# ==========================================================================
# Modulo Shiny: una solapa completa por analito
# ==========================================================================
analyteUI <- function(id, cfg) {
  ns <- NS(id)
  tagList(
    layout_columns(
      col_widths = c(4, 8),
      # ---------------- Columna izquierda: alcance y descargas ----------------
      card(
        card_header(style = sprintf("background:%s;color:white;font-weight:600;", cfg$color),
                     sprintf("%s (%s) — Alcance y trazabilidad", cfg$label, cfg$symbol)),
        card_body(
          tags$table(class = "table table-sm",
            tags$tr(tags$td(tags$b("Rango validado")), tags$td(cfg$scope_range)),
            tags$tr(tags$td(tags$b("Muestras")), tags$td(cfg$n_samples)),
            tags$tr(tags$td(tags$b("Modelo ganador")), tags$td(cfg$model_desc))
          ),
          if (isTRUE(cfg$spiked))
            div(class = "alert alert-secondary", style = "font-size:0.9em",
                "Parte del rango de validacion se obtuvo fortificando vinos comerciales con el metal de interes, ",
                "para cubrir concentraciones poco frecuentes en la matriz natural. Ver ficha del modelo para el detalle."),
          tags$hr(),
          tags$b("Variedades y origen de los vinos"), tags$p(COMMON_INFO$wine_scope, style = "font-size:0.9em"),
          tags$b("Instrumento FT-MIR"), tags$p(COMMON_INFO$instrument, style = "font-size:0.9em"),
          tags$b("Metodo de referencia (valores reales)"), tags$p(COMMON_INFO$reference_method, style = "font-size:0.9em"),
          tags$hr(),
          accordion(
            open = FALSE,
            accordion_panel("Protocolo de preparacion de la muestra (Na₂EDTA)",
                             p(COMMON_INFO$edta_protocol))
          ),
          tags$hr(),
          downloadButton(ns("dl_sample"), "Descargar dataset de validacion", class = "btn-outline-secondary btn-sm w-100 mb-2"),
          downloadButton(ns("dl_model_html"), "Descargar ficha resumen del modelo", class = "btn-outline-secondary btn-sm w-100"),
          tags$hr(),
          div(class = "alert alert-warning", style = "font-size:0.85em", COMMON_INFO$disclaimer)
        )
      ),
      # ---------------- Columna derecha: validacion + prediccion ----------------
      div(
        card(
          card_header("Estado de validacion del modelo"),
          card_body(
            actionButton(ns("btn_validate"), "Ver estado de validacion", class = "btn-primary mb-3"),
            uiOutput(ns("validation_ui"))
          )
        ),
        card(
          card_header("Predecir a partir de un espectro propio"),
          card_body(
            p("Subi un espectro FT-MIR (902.57 a 3000.84 cm⁻¹). Se acepta .csv o .xlsx con dos columnas ",
              "(numero de onda, absorbancia), o el mismo formato de fila que los datasets del estudio.",
              style = "font-size:0.9em; color:#555"),
            fileInput(ns("spectrum_file"), NULL, accept = c(".csv", ".xlsx")),
            actionButton(ns("btn_predict"), "Predecir", class = "btn-primary mb-3"),
            uiOutput(ns("prediction_ui"))
          )
        )
      )
    )
  )
}

analyteServer <- function(id, cfg) {
  moduleServer(id, function(input, output, session) {

    bundle <- reactive(load_bundle(id))

    # ---------------- Descargas ----------------
    output$dl_sample <- downloadHandler(
      filename = function() sprintf("Muestras_validacion_%s.xlsx", cfg$symbol),
      content = function(file) file.copy(cfg$sample_file, file)
    )
    output$dl_model_html <- downloadHandler(
      filename = function() sprintf("Ficha_modelo_%s.html", cfg$symbol),
      content = function(file) {
        b <- bundle()
        validate(need(!is.null(b), "El modelo todavia no esta disponible."))
        writeLines(render_model_report_html(cfg, b), file)
      }
    )

    # ---------------- Estado de validacion ----------------
    output$validation_ui <- renderUI({
      req(input$btn_validate > 0)
      b <- bundle()
      ns <- session$ns
      if (is.null(b)) {
        return(div(class = "alert alert-danger",
                    "El modelo de este analito todavia no esta cargado en el servidor."))
      }
      if (cfg$task == "regression") {
        tagList(
          h5("Hiperparametros del modelo ganador"),
          p(paste(sprintf("%s = %s", names(b$hyperparams), unlist(b$hyperparams)), collapse = "; ")),
          DTOutput(ns("metrics_table")),
          plotOutput(ns("scatter_plot"), height = "380px")
        )
      } else {
        tagList(
          h5("Hiperparametros del modelo ganador"),
          p(paste(sprintf("%s = %s", names(b$hyperparams), unlist(b$hyperparams)), collapse = "; ")),
          DTOutput(ns("metrics_table")),
          plotOutput(ns("confusion_plot"), height = "320px"),
          tags$hr(),
          h5("Agrupamiento jerarquico de las muestras (no supervisado)"),
          p("Compara si los grupos que forma el espectro por si solo (color de rama) coinciden con la clase real (color de etiqueta).",
            style = "font-size:0.85em;color:#666"),
          plotOutput(ns("dendrogram_plot"), height = "480px")
        )
      }
    })

    output$metrics_table <- renderDT({
      req(input$btn_validate > 0)
      b <- bundle(); validate(need(!is.null(b), ""))
      m <- b$metrics
      df <- if (cfg$task == "regression") {
        data.frame(Metrica = c("R²", "RMSE", "RPD"),
                   Entrenamiento = c(round(m$R2_Train, 3), round(m$RMSE_Train, 2), round(m$RPD_Train, 2)),
                   Test = c(round(m$R2_Test, 3), round(m$RMSE_Test, 2), round(m$RPD_Test, 2)))
      } else {
        data.frame(Metrica = c("AUC", "Accuracy", "Kappa", "Sensibilidad", "Especificidad"),
                   Test = c(round(m$AUC_Test, 3), round(m$Accuracy_Test, 3), round(m$Kappa_Test, 3),
                            round(m$Sensitivity_Test, 3), round(m$Specificity_Test, 3)))
      }
      datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
    })

    output$scatter_plot <- renderPlot({
      req(input$btn_validate > 0, cfg$task == "regression")
      b <- bundle(); validate(need(!is.null(b), ""))
      ggplot(b$scatter_data, aes(x = Actual, y = Predicted, color = Set)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
        geom_point(size = 2.5, alpha = 0.8) +
        scale_color_manual(values = c(Train = "#457B9D", Test = cfg$color)) +
        labs(x = sprintf("%s medido (%s)", cfg$label, cfg$unit),
             y = sprintf("%s predicho (%s)", cfg$label, cfg$unit), color = NULL) +
        theme_minimal(base_size = 13)
    })

    output$confusion_plot <- renderPlot({
      req(input$btn_validate > 0, cfg$task == "classification")
      b <- bundle(); validate(need(!is.null(b), ""))
      te <- b$scatter_data
      cm <- as.data.frame(table(Actual = te$Actual, Predicted = te$Predicted))
      ggplot(cm, aes(x = Actual, y = Predicted, fill = Freq)) +
        geom_tile(color = "white") +
        geom_text(aes(label = Freq), size = 6, fontface = "bold") +
        scale_fill_gradient(low = "#f2f2f2", high = cfg$color) +
        labs(title = "Matriz de confusion (test)") +
        theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
    })

    output$dendrogram_plot <- renderPlot({
      req(input$btn_validate > 0, cfg$task == "classification")
      b <- bundle(); validate(need(!is.null(b), ""))
      validate(need(!is.null(b$dendro_hc), "Dendrograma no disponible para este modelo."))
      render_dendrogram(b, cfg$color)
    })

    # ---------------- Prediccion sobre espectro nuevo ----------------
    pred_result <- eventReactive(input$btn_predict, {
      b <- bundle()
      validate(need(!is.null(b), "El modelo de este analito todavia no esta disponible."))
      validate(need(!is.null(input$spectrum_file), "Subi un espectro primero."))
      spec <- tryCatch(read_uploaded_spectrum(input$spectrum_file$datapath),
                        error = function(e) NULL)
      validate(need(!is.null(spec), "No se pudo leer el archivo. Verifica el formato (numero de onda / absorbancia)."))
      predict_spectrum(b, spec$wn, spec$ab)
    })

    output$prediction_ui <- renderUI({
      res <- pred_result()
      ns <- session$ns
      note <- if (isTRUE(res$resampled))
        div(class = "alert alert-info", style = "font-size:0.85em",
            "El espectro subido no coincidia exactamente con la grilla de longitudes de onda del modelo; ",
            "se realino por interpolacion antes de predecir.")
      else NULL

      if (res$task == "regression") {
        tagList(
          note,
          div(class = "p-3 mb-2", style = sprintf("background:%s10;border-left:4px solid %s;", cfg$color, cfg$color),
              h3(sprintf("%.1f %s", res$value, cfg$unit), style = sprintf("color:%s", cfg$color)),
              p(sprintf("Intervalo aproximado (95%%): %.1f – %.1f %s", res$lower, res$upper, cfg$unit)),
              p(sprintf("Incertidumbre: ± 2×RMSE de test (± %.1f %s)", res$uncertainty / 2, cfg$unit),
                style = "font-size:0.85em;color:#666")
          ),
          downloadButton(ns("dl_prediction"), "Descargar resultado (informe)", class = "btn-outline-secondary btn-sm")
        )
      } else {
        probs <- res$probabilities
        tagList(
          note,
          div(class = "p-3 mb-2", style = sprintf("background:%s10;border-left:4px solid %s;", cfg$color, cfg$color),
              h3(res$class, style = sprintf("color:%s", cfg$color)),
              lapply(names(probs), function(cl) {
                pct <- round(as.numeric(probs[[cl]]) * 100, 1)
                div(style = "margin-bottom:6px;",
                    span(sprintf("%s: %.1f%%", cl, pct)),
                    div(style = sprintf("background:#eee;border-radius:4px;height:10px;width:100%%;"),
                        div(style = sprintf("background:%s;border-radius:4px;height:10px;width:%s%%;", cfg$color, pct))))
              })
          ),
          downloadButton(ns("dl_prediction"), "Descargar resultado (informe)", class = "btn-outline-secondary btn-sm")
        )
      }
    })

    output$dl_prediction <- downloadHandler(
      filename = function() sprintf("Prediccion_%s_%s.html", cfg$symbol, format(Sys.time(), "%Y%m%d_%H%M")),
      content = function(file) writeLines(render_prediction_report_html(cfg, pred_result()), file)
    )
  })
}

# ==========================================================================
# Tema visual
# ==========================================================================
app_theme <- bs_theme(
  version = 5, bootswatch = "flatly",
  primary = "#7B1E3A", base_font = font_google("Inter"),
  heading_font = font_google("Source Serif Pro")
)

# ==========================================================================
# Pagina de inicio
# ==========================================================================
home_panel <- div(
  class = "p-4",
  div(style = "display:flex;align-items:center;gap:24px;margin-bottom:1.5em;flex-wrap:wrap;",
      if (file.exists("www/ivagro_logo.png")) tags$img(src = "ivagro_logo.png", height = "80px"),
      if (file.exists("www/uca_logo.png")) tags$img(src = "uca_logo.png", height = "50px"),
      div(
        h2("Scan FT-MIR", style = "margin-bottom:0;color:#7B1E3A;"),
        h5("Estimacion de metales en vino por espectroscopia FT-MIR y aprendizaje automatico",
           style = "font-weight:400;color:#555;")
      )
  ),
  p("Esta herramienta pone a disposicion los modelos de aprendizaje automatico desarrollados para estimar ",
    "la concentracion de potasio, magnesio y calcio, y para clasificar hierro y cobre respecto a sus limites ",
    "de interes enologico (10 y 1 mg/L respectivamente), a partir de un espectro FT-MIR de la muestra de vino ",
    "preparada con Na2EDTA. Cada solapa corresponde a un analito, con su propio modelo validado, su alcance, ",
    "y un panel para subir un espectro propio."),
  div(class = "alert alert-warning", COMMON_INFO$disclaimer),
  p(tags$b("Desarrollado en el "),
    tags$a(href = IVAGRO_URL, target = "_blank",
           "Instituto de Investigacion Vitivinicola y Agroalimentaria (IVAGRO), Universidad de Cadiz"), "."),
  p(class = "text-muted", style = "font-size:0.85em",
    "Codigo y datos: ", tags$a(href = "https://github.com/pablofacu86/Wine-Metals-FTMIR-ML", target = "_blank",
                                "github.com/pablofacu86/Wine-Metals-FTMIR-ML"))
)

# ==========================================================================
# UI principal
# ==========================================================================
ui <- page_navbar(
  title = "Scan FT-MIR",
  theme = app_theme,
  fillable = FALSE,
  nav_panel("Inicio", home_panel),
  nav_panel("Potasio (K)",   analyteUI("POTASSIUM", ANALYTES$POTASSIUM)),
  nav_panel("Magnesio (Mg)", analyteUI("MAGNESIUM", ANALYTES$MAGNESIUM)),
  nav_panel("Calcio (Ca)",   analyteUI("CALCIUM",   ANALYTES$CALCIUM)),
  nav_panel("Hierro (Fe)",   analyteUI("IRON",      ANALYTES$IRON)),
  nav_panel("Cobre (Cu)",    analyteUI("COPPER",    ANALYTES$COPPER)),
  nav_spacer(),
  nav_item(tags$a(href = IVAGRO_URL, target = "_blank", "IVAGRO — Universidad de Cadiz"))
)

# ==========================================================================
# Server principal
# ==========================================================================
server <- function(input, output, session) {
  for (k in names(ANALYTES)) analyteServer(k, ANALYTES[[k]])
}

shinyApp(ui, server)
