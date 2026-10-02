# ==========================================================================
# predict_utils.R
# Funciones para preprocesar un espectro FT-MIR nuevo exactamente de la
# misma manera que se preproceso durante el entrenamiento (ver
# Build_Deployment_Models.R), y para generar una prediccion a partir de
# un modelo ya entrenado (bundle .rds).
# ==========================================================================

library(prospectr)
library(caret)   # necesario para predict() sobre los modelos entrenados (glmnet, pls, svmRadial)
library(xgboost) # necesario para predict() sobre el modelo de calcio (XGBoost)
library(dendextend)
library(circlize)

parse_wn <- function(x) as.numeric(gsub(",", ".", gsub("X", "", x), fixed = TRUE))

PATZ_WINDOWS <- list(c(965, 1582), c(1698, 2006), c(2701, 2971))
in_patz_windows <- function(w) {
  Reduce(`|`, lapply(PATZ_WINDOWS, function(rng) w >= rng[1] & w <= rng[2]))
}

# --------------------------------------------------------------------
# Lee un espectro subido por el usuario (.csv o .xlsx con dos columnas:
# numero de onda y absorbancia, en cualquier orden y con cualquier
# nombre de columna) y devuelve una lista list(wn=..., ab=...) ordenada
# por numero de onda creciente.
# --------------------------------------------------------------------
read_uploaded_spectrum <- function(filepath) {
  ext <- tolower(tools::file_ext(filepath))
  df <- if (ext %in% c("xlsx", "xls")) {
    openxlsx::read.xlsx(filepath, sheet = 1)
  } else {
    first_line <- readLines(filepath, n = 1)
    sep_char <- if (grepl(";", first_line)) ";" else ","
    read.csv(filepath, header = TRUE, sep = sep_char)
  }
  if (ncol(df) < 2) stop("El archivo debe tener al menos dos columnas: numero de onda y absorbancia.")
  # Si el archivo tiene 547 columnas (mismo formato que el dataset: ID +
  # Grupo + 545 variables), se toma directamente la fila de espectro.
  if (ncol(df) > 10) {
    wn <- suppressWarnings(as.numeric(gsub(",", ".", colnames(df)[-c(1, 2)], fixed = TRUE)))
    ab <- as.numeric(df[1, -c(1, 2)])
    ok <- !is.na(wn)
    wn <- wn[ok]; ab <- ab[ok]
    o  <- order(wn)
    return(list(wn = wn[o], ab = ab[o]))
  }
  # Formato de dos columnas: numero de onda / absorbancia
  num1 <- suppressWarnings(as.numeric(df[[1]])); num2 <- suppressWarnings(as.numeric(df[[2]]))
  if (mean(!is.na(num1)) > mean(!is.na(num2))) { wn <- num1; ab <- num2 } else { wn <- num2; ab <- num1 }
  ok <- !is.na(wn) & !is.na(ab)
  o  <- order(wn[ok])
  list(wn = wn[ok][o], ab = ab[ok][o])
}

# --------------------------------------------------------------------
# Alinea un espectro nuevo (wn, ab) a la grilla de referencia usada para
# entrenar el modelo (bundle$wavelengths_full: 545 puntos, 902.57-3000.84
# cm-1, cada ~3.86 cm-1), por interpolacion lineal. Esto permite usar
# espectros de otros equipos, con otro paso o rango de muestreo, siempre
# que cubran razonablemente la zona 900-3000 cm-1. Si la grilla ya
# coincide (mismo instrumento), la interpolacion no cambia los valores.
# --------------------------------------------------------------------
align_to_reference_grid <- function(wn, ab, ref_wn) {
  same_grid <- length(wn) == length(ref_wn) && max(abs(wn - ref_wn)) < 0.05
  coverage <- mean(ref_wn >= min(wn) & ref_wn <= max(wn))
  ab_aligned <- approx(x = wn, y = ab, xout = ref_wn, rule = 2)$y
  list(ab = ab_aligned, resampled = !same_grid, coverage = coverage)
}

# --------------------------------------------------------------------
# Preprocesa un espectro ya alineado a la grilla de referencia, con la
# misma secuencia usada al entrenar: derivada (si corresponde) sobre el
# espectro completo, correccion de dispersion (SNV o MSC con la
# referencia guardada del training), recorte a las ventanas de Patz,
# centrado con la media del training, y seleccion de variables de
# Boruta (si el modelo ganador la uso). Devuelve una fila lista para
# predict().
# --------------------------------------------------------------------
preprocess_for_model <- function(bundle, ab_aligned) {
  X <- matrix(ab_aligned, nrow = 1)
  colnames(X) <- paste0("X", bundle$wavelengths_full)

  if (bundle$deriv > 0) {
    Xd <- savitzkyGolay(X, m = bundle$deriv, p = 3, w = 11)
    if (ncol(Xd) == ncol(X)) colnames(Xd) <- colnames(X)
    X <- Xd
  }
  if (bundle$scatter == "SNV") {
    X <- standardNormalVariate(X)
  }
  if (bundle$scatter == "MSC") {
    ref <- bundle$msc_ref
    fit <- lm(as.numeric(X[1, ]) ~ ref)
    X[1, ] <- (as.numeric(X[1, ]) - coef(fit)[1]) / coef(fit)[2]
  }

  wo <- parse_wn(colnames(X))
  keep_patz <- in_patz_windows(wo)
  X <- X[, keep_patz, drop = FALSE]

  # Centrado con la media del training (mismo orden de variables)
  if (!identical(colnames(X), names(bundle$mean_vec))) {
    # reordena por si acaso, matcheando por nombre
    common <- intersect(colnames(X), names(bundle$mean_vec))
    X <- X[, common, drop = FALSE]
    mu <- bundle$mean_vec[common]
  } else {
    mu <- bundle$mean_vec
  }
  X <- sweep(X, 2, mu, "-")

  if (isTRUE(bundle$use_boruta) && !is.null(bundle$boruta_vars)) {
    X <- X[, bundle$boruta_vars, drop = FALSE]
  }
  as.data.frame(X, check.names = TRUE)
}

# --------------------------------------------------------------------
# Funcion principal: toma un bundle (modelo .rds cargado) y un espectro
# nuevo (wn, ab), y devuelve la prediccion lista para mostrar en la app.
# --------------------------------------------------------------------
predict_spectrum <- function(bundle, wn, ab) {
  aligned <- align_to_reference_grid(wn, ab, bundle$wavelengths_full)
  if (aligned$coverage < 0.90) {
    warning(sprintf(
      "El espectro cubre solo %.0f%% del rango 902.57-3000.84 cm-1 esperado; la prediccion puede no ser confiable.",
      aligned$coverage * 100))
  }
  newrow <- preprocess_for_model(bundle, aligned$ab)

  if (bundle$task == "regression") {
    pred <- if (identical(bundle$algo, "XGB")) {
      as.numeric(predict(bundle$model, as.matrix(newrow)))
    } else {
      as.numeric(predict(bundle$model, newrow))
    }
    unc  <- 2 * bundle$metrics$RMSE_Test   # +/- 2*RMSE(test) ~ intervalo aprox. al 95%
    list(task = "regression", value = pred, lower = pred - unc, upper = pred + unc,
         uncertainty = unc, resampled = aligned$resampled, coverage = aligned$coverage)
  } else {
    cls  <- as.character(predict(bundle$model, newrow))
    prob <- predict(bundle$model, newrow, type = "prob")
    list(task = "classification", class = cls, probabilities = prob,
         pos_class = bundle$pos_class, resampled = aligned$resampled, coverage = aligned$coverage)
  }
}

# --------------------------------------------------------------------
# Redibuja, en vivo, el dendrograma circular (Ward.D2, ramas coloreadas
# por cluster no supervisado, etiquetas coloreadas por clase real) a
# partir de los datos guardados en el bundle (dendro_hc, dendro_y),
# sin depender de ninguna imagen fija. Devuelve NULL si el bundle no
# tiene datos de dendrograma (por ejemplo, analitos de regresion).
# --------------------------------------------------------------------
render_dendrogram <- function(bundle, accent_color) {
  if (is.null(bundle$dendro_hc) || is.null(bundle$dendro_y)) return(invisible(NULL))
  on.exit(try(circlize::circos.clear(), silent = TRUE))

  hc    <- bundle$dendro_hc
  order <- hc$order
  y_ord <- as.character(bundle$dendro_y)[order]
  classes <- bundle$class_names
  if (is.null(classes)) classes <- unique(y_ord)

  cls_colors <- setNames(c(accent_color, "#457B9D")[seq_along(classes)], classes)
  class_short <- setNames(sub("^([A-Za-z]+).*", "\\1", gsub("\\.", " ", classes)), classes)
  leaf_short   <- class_short[y_ord]
  leaf_counter <- ave(seq_along(leaf_short), leaf_short, FUN = seq_along)
  new_labels   <- sprintf("%s_%02d", leaf_short, leaf_counter)
  leaf_colors  <- cls_colors[y_ord]

  dend <- as.dendrogram(hc)
  labels_colors(dend) <- leaf_colors
  labels(dend) <- new_labels
  dend <- color_branches(dend, k = 4)
  dend <- set(dend, "branches_lwd", 2)
  dend <- set(dend, "labels_cex", 0.9)

  circlize_dendrogram(dend, labels_track_height = 0.28, dend_track_height = 0.55)
  legend("bottomright", legend = gsub("\\.", " ", names(cls_colors)), text.col = cls_colors,
         bty = "n", cex = 0.9)
}
