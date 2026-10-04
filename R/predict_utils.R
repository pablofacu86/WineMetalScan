# ==========================================================================
# predict_utils.R
# Functions to preprocess a new FT-MIR spectrum in exactly the same way
# it was preprocessed during training (see Build_Deployment_Models.R),
# and to generate a prediction from an already-trained model (.rds
# bundle).
# ==========================================================================

library(prospectr)
library(caret)   # needed for predict() on the trained models
# caret calls these packages internally when predicting, but they must be
# attached/declared so that the deployment (Posit Connect) installs them:
library(glmnet)  # potassium model (Elastic Net)
library(pls)     # magnesium model (PLS)
library(kernlab) # iron and copper models (SVM radial, via caret)
library(xgboost) # needed for predict() on the calcium model (XGBoost)
library(dendextend)
library(circlize)

parse_wn <- function(x) as.numeric(gsub(",", ".", gsub("X", "", x), fixed = TRUE))

PATZ_WINDOWS <- list(c(965, 1582), c(1698, 2006), c(2701, 2971))
in_patz_windows <- function(w) {
  Reduce(`|`, lapply(PATZ_WINDOWS, function(rng) w >= rng[1] & w <= rng[2]))
}

# --------------------------------------------------------------------
# Reads the file uploaded by the user and returns
#   list(ids = <sample names>, wn = <wavenumbers>, ab = <matrix, one row per sample>)
#
# Two layouts are accepted (.csv or .xlsx):
#   1) Long format, ONE sample: two columns (wavenumber, absorbance), in
#      any order and with any column names.
#   2) Wide format, one row per sample: the first column is the sample
#      name and every column whose header is a number (the wavenumber) is
#      taken as spectral data. Any other column (e.g. a "Group" or
#      reference-value column present in the study datasets) is ignored,
#      so it is NOT required and is never used for prediction.
# If an .xlsx file has several sheets, the sheet whose name matches the
# analyte (sheet_hint) is used; otherwise the first sheet.
# --------------------------------------------------------------------
read_uploaded_spectra <- function(filepath, sheet_hint = NULL) {
  ext <- tolower(tools::file_ext(filepath))
  if (ext == "xls") stop("Old .xls files are not supported: please save the file as .xlsx or .csv.")
  if (ext %in% c("xlsx", "xls")) {
    sheets <- openxlsx::getSheetNames(filepath)
    sheet  <- if (!is.null(sheet_hint) && sheet_hint %in% sheets) sheet_hint else sheets[1]
    df <- openxlsx::read.xlsx(filepath, sheet = sheet, check.names = FALSE)
  } else {
    first_line <- readLines(filepath, n = 1, warn = FALSE)
    sep_char <- if (grepl(";", first_line)) ";" else ","
    df <- read.csv(filepath, header = TRUE, sep = sep_char, check.names = FALSE)
  }
  if (ncol(df) < 2) stop("The file must have at least two columns.")

  # numeric header -> candidate spectral column (wide format)
  hdr_wn <- suppressWarnings(as.numeric(gsub(",", ".", gsub("^X", "", colnames(df)))))
  spec_cols <- which(!is.na(hdr_wn) & hdr_wn > 500 & hdr_wn < 4500)

  if (length(spec_cols) > 10) {
    wn  <- hdr_wn[spec_cols]
    # values may use a decimal comma (e.g. "0,0123"): convert column by column
    vals <- lapply(df[, spec_cols, drop = FALSE],
                   function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x)))))
    ab  <- matrix(unlist(vals), nrow = nrow(df))
    ids <- if (ncol(df) > length(spec_cols)) as.character(df[[setdiff(seq_len(ncol(df)), spec_cols)[1]]])
           else paste("Sample", seq_len(nrow(df)))
    ids[is.na(ids) | trimws(ids) == ""] <- paste("Sample", seq_along(ids))[is.na(ids) | trimws(ids) == ""]
    # drop completely empty rows (e.g. blank lines at the end of a sheet)
    ok_rows <- rowSums(!is.na(ab)) > 10
    if (!any(ok_rows)) stop("No spectra with numeric values were found in the file.")
    ids <- ids[ok_rows]; ab <- ab[ok_rows, , drop = FALSE]
    o <- order(wn)
    return(list(ids = ids, wn = wn[o], ab = ab[, o, drop = FALSE]))
  }

  # long format: wavenumber / absorbance (one sample)
  num1 <- suppressWarnings(as.numeric(gsub(",", ".", df[[1]])))
  num2 <- suppressWarnings(as.numeric(gsub(",", ".", df[[2]])))
  if (mean(!is.na(num1)) < 0.9 && mean(!is.na(num2)) < 0.9)
    stop("Could not find numeric wavenumber/absorbance columns in the file.")
  # wavenumbers (~900-3000) are much larger than absorbances (~0-5)
  if (mean(num1, na.rm = TRUE) > mean(num2, na.rm = TRUE)) { wn <- num1; ab <- num2 } else { wn <- num2; ab <- num1 }
  ok <- !is.na(wn) & !is.na(ab)
  o  <- order(wn[ok])
  list(ids = "Uploaded spectrum", wn = wn[ok][o], ab = matrix(ab[ok][o], nrow = 1))
}

# --------------------------------------------------------------------
# Aligns a new spectrum (wn, ab) to the reference grid used to train
# the model (bundle$wavelengths_full: 545 points, 902.57-3000.84 cm-1,
# ~3.86 cm-1 apart), by linear interpolation. This allows using spectra
# from other instruments, with a different sampling step or range,
# provided they reasonably cover the 900-3000 cm-1 region. If the grid
# already matches (same instrument), interpolation does not change the
# values.
# --------------------------------------------------------------------
align_to_reference_grid <- function(wn, ab, ref_wn) {
  ok <- !is.na(wn) & !is.na(ab)
  wn <- wn[ok]; ab <- ab[ok]
  if (length(wn) < 20) stop("The spectrum has too few valid points to be used.")
  same_grid <- length(wn) == length(ref_wn) && max(abs(wn - ref_wn)) < 0.05
  coverage <- mean(ref_wn >= min(wn) & ref_wn <= max(wn))
  # The models only use the Patz windows (965-1582, 1698-2006, 2701-2971
  # cm-1); the spectrum must cover all of them, otherwise the missing
  # region would be filled with constants and the prediction would be
  # meaningless.
  patz_pts <- ref_wn[in_patz_windows(ref_wn)]
  patz_cov <- mean(patz_pts >= min(wn) & patz_pts <= max(wn))
  if (patz_cov < 0.95) {
    stop(sprintf(paste0("The spectrum covers %.0f-%.0f cm-1, but the model needs data between ",
                        "about 965 and 2971 cm-1. Please upload a spectrum covering 900-3000 cm-1."),
                 min(wn), max(wn)))
  }
  ab_aligned <- approx(x = wn, y = ab, xout = ref_wn, rule = 2, ties = mean)$y
  list(ab = ab_aligned, resampled = !same_grid, coverage = coverage)
}

# --------------------------------------------------------------------
# Preprocesses a spectrum already aligned to the reference grid, with
# the same sequence used in training: derivative (if applicable) on
# the full spectrum, scatter correction (SNV or MSC with the training
# reference), cropping to the Patz windows, centring with the training
# mean, and Boruta variable selection (if the winning model used it).
# Returns a row ready for predict().
# --------------------------------------------------------------------
preprocess_for_model <- function(bundle, ab_aligned) {
  # Some prospectr functions are not guaranteed to keep a single-row matrix
  # as a matrix; applying them to two identical rows and keeping the first
  # one avoids that problem.
  one_row <- function(f, X) f(rbind(X, X))[1, , drop = FALSE]

  X <- matrix(ab_aligned, nrow = 1)
  colnames(X) <- paste0("X", bundle$wavelengths_full)

  if (bundle$deriv > 0) {
    cn <- colnames(X)
    Xd <- one_row(function(M) savitzkyGolay(M, m = bundle$deriv, p = 3, w = 11), X)
    # name the columns by position: if the filter trims the edges, the
    # remaining columns are the central ones
    half <- (length(cn) - ncol(Xd)) %/% 2
    colnames(Xd) <- cn[(half + 1):(half + ncol(Xd))]
    X <- Xd
  }
  if (bundle$scatter == "SNV") {
    cn2 <- colnames(X)
    X <- one_row(standardNormalVariate, X)
    colnames(X) <- cn2
  }
  if (bundle$scatter == "MSC") {
    ref <- as.numeric(bundle$msc_ref)
    if (length(ref) != ncol(X)) stop("Internal error: MSC reference length does not match the spectrum.")
    fit <- lm(as.numeric(X[1, ]) ~ ref)
    X[1, ] <- (as.numeric(X[1, ]) - coef(fit)[1]) / coef(fit)[2]
  }

  wo <- parse_wn(colnames(X))
  keep_patz <- in_patz_windows(wo)
  X <- X[, keep_patz, drop = FALSE]

  # Centring with the training mean. The variables are matched by POSITION
  # (the grid and the Patz cropping are deterministic), and then given
  # exactly the names the model was trained with. This avoids any mismatch
  # caused by how the wavenumbers are formatted as text.
  mu <- bundle$mean_vec
  if (ncol(X) != length(mu)) {
    stop(sprintf("Internal error: the preprocessed spectrum has %d variables but the model expects %d.",
                 ncol(X), length(mu)))
  }
  colnames(X) <- names(mu)
  X <- sweep(X, 2, as.numeric(mu), "-")

  if (isTRUE(bundle$use_boruta) && !is.null(bundle$boruta_vars)) {
    X <- X[, bundle$boruta_vars, drop = FALSE]
  }
  as.data.frame(X, check.names = TRUE)
}

# --------------------------------------------------------------------
# Main function: takes a bundle (loaded .rds model) and a new spectrum
# (wn, ab), and returns the prediction ready to display in the app.
# --------------------------------------------------------------------
predict_spectrum <- function(bundle, wn, ab) {
  aligned <- align_to_reference_grid(wn, ab, bundle$wavelengths_full)
  newrow <- preprocess_for_model(bundle, aligned$ab)

  if (bundle$task == "regression") {
    pred <- if (identical(bundle$algo, "XGB")) {
      as.numeric(predict(bundle$model, as.matrix(newrow)))
    } else {
      as.numeric(predict(bundle$model, newrow))
    }
    unc  <- 2 * bundle$metrics$RMSE_Test   # +/- 2*test RMSE ~ approx. 95% interval
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
# Redraws, live, the circular dendrogram (Ward.D2, branches coloured by
# unsupervised cluster, labels coloured by actual class) from the data
# stored in the bundle (dendro_hc, dendro_y), without depending on any
# fixed image. Returns NULL if the bundle has no dendrogram data (for
# example, regression analytes).
#
# limit_value/unit are optional and only used to build a clean legend
# (e.g. "Lower than 10 mg/L"); if not supplied, the legend falls back
# to a best-effort cleanup of the stored class names.
# --------------------------------------------------------------------
render_dendrogram <- function(bundle, accent_color, limit_value = NULL, unit = NULL) {
  if (is.null(bundle$dendro_hc) || is.null(bundle$dendro_y)) return(invisible(NULL))
  on.exit(try(circlize::circos.clear(), silent = TRUE))

  hc         <- bundle$dendro_hc
  leaf_order <- hc$order
  y_ord      <- as.character(bundle$dendro_y)[leaf_order]
  classes    <- bundle$class_names
  if (is.null(classes)) classes <- unique(y_ord)

  cls_colors  <- setNames(c(accent_color, "#457B9D")[seq_along(classes)], classes)
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

  # Legend text: reconstruct "Lower/Higher than X unit" from the known
  # limit when available, instead of just de-dotting the sanitised
  # class name (which would otherwise lose the "/" in "mg/L").
  legend_labels <- if (!is.null(limit_value) && !is.null(unit)) {
    sprintf("%s than %g %s", class_short[classes], limit_value, unit)
  } else {
    gsub("\\.", " ", names(cls_colors))
  }
  legend("bottomright", legend = legend_labels, text.col = cls_colors, bty = "n", cex = 0.9)
}
