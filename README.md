# WineMetalScan — prediction app (Posit Connect Cloud)

Shiny app to estimate potassium, magnesium and calcium (regression) and to
classify iron and copper (above / below their limit of interest) in wine
from an FT-MIR spectrum, using the models validated in the study.

Study repository (pipelines, raw spectra, reference data):
<https://github.com/pablofacu86/Wine-Metals-FTMIR-ML> —
archived version 1.1.0: <https://doi.org/10.5281/zenodo.23199425>

## Models in the app

| Analyte | Task | Model | Preprocessing | Variable selection |
|---|---|---|---|---|
| Potassium (K) | regression | Lasso | Savitzky–Golay smoothing (SG0) | none (310 variables) |
| Magnesium (Mg) | regression | PLS (7 latent variables) | SNV → 1st derivative (SG1) | none (310 variables) |
| Calcium (Ca) | regression | XGBoost | 1st derivative (SG1) | Boruta (59 variables) |
| Iron (Fe) | classification (10 mg/L) | SVM, radial kernel | MSC → smoothing (SG0) | Boruta (46 variables) |
| Copper (Cu) | classification (1 mg/L) | XGBoost | SNV → 1st derivative (SG1) | Boruta (53 variables) |

Savitzky–Golay: polynomial order 3, window of 11 points. Scatter correction
(SNV / MSC) is applied **before** Savitzky–Golay, always on the full spectrum
(545 variables); the Patz windows (965–1582, 1698–2006 and 2701–2971 cm⁻¹,
310 variables) are cropped afterwards and the spectrum is mean-centred with
the training mean. In classification the positive class is "Higher than
the limit".

## Expected structure

```
app.R
R/predict_utils.R
Build_Deployment_Models.R   <- builds Models/*.rds (run locally, not deployed)
Models/                     <- POTASSIUM_model.rds, MAGNESIUM_model.rds, ... (generated)
data/
  FINAL_DATA_SET.xlsx       <- only needed to run Build_Deployment_Models.R
  Test_Set_App.xlsx         <- practice samples with known reference values
www/
  ivagro_logo.png
  uca_logo.png
manifest.json               <- generated with rsconnect::writeManifest()
```

## Rebuilding the models

1. Put `FINAL_DATA_SET.xlsx` (from the study repository, `data/` folder) in
   `data/`.
2. Install the required packages (see below) and run, from the app folder:

   ```r
   source("Build_Deployment_Models.R")
   ```

3. The script re-trains only the winning combination of each analyte (it does
   not repeat the full scan), saves `Models/<ANALYTE>_model.rds` and prints,
   for each analyte, the test-set metrics of the rebuilt model next to the
   values reported in the article (also saved in `Models/build_summary.csv`).
   Small differences are expected for XGBoost (random subsampling) and for
   the SVM (the kernel parameter is stored with five decimals); a `CHECK`
   flag means a metric differs by more than the tolerance and should be
   looked at before deploying.

The Boruta variables of the calcium, iron and copper models are the ones
reported in the study (`BORUTA_MODE <- "saved"`). Setting
`BORUTA_MODE <- "rerun"` runs Boruta again on the training set, but the
selected variables may then differ slightly because Boruta is stochastic.

## Deploying

After the five `.rds` files are in `Models/`:

```r
rsconnect::writeManifest()   # regenerates manifest.json with the package versions
```

Commit `app.R`, `R/`, `Models/`, `data/Test_Set_App.xlsx`, `www/` and
`manifest.json` to the app repository and redeploy from Posit Connect Cloud.
`Build_Deployment_Models.R` and `data/FINAL_DATA_SET.xlsx` are not needed by
the running app.

Bundles built with the previous version of the pipeline (no
`scatter_first` field) are still read correctly by `predict_utils.R`, but
they must be replaced: the models of the study changed.

## Required R packages

`shiny`, `bslib` (>= 0.5), `DT`, `ggplot2`, `openxlsx`, `prospectr`, `caret`,
`glmnet`, `pls`, `kernlab`, `xgboost`, `dendextend`, `circlize`. To rebuild the
models you also need `pROC` (and `Boruta` only with `BORUTA_MODE <- "rerun"`).

Each bundle stores the R version and the versions of the main packages used
to build it (`bundle$r_version`, `bundle$package_versions`).

## How prediction on a new spectrum works

1. The user uploads a file (.csv or .xlsx) with wavenumber and absorbance for
   one spectrum, or a file with one row per sample in the same layout as the
   study data (sample name, optional reference value, spectral columns).
2. The spectrum is **realigned by linear interpolation** to the model's
   reference grid (902.57–3000.84 cm⁻¹, ~3.86 cm⁻¹ step). This allows spectra
   from another instrument with a different step or range, provided they
   cover the Patz windows (965–2971 cm⁻¹); otherwise the app stops with a
   message.
3. The same preprocessing as in training is applied: scatter correction
   (SNV, or MSC with the stored training reference) → Savitzky–Golay →
   Patz windows → centring with the training mean → Boruta variables (when
   the winning model uses them).
4. The model predicts. For regression, an approximate 95 % interval is
   reported as value ± 2 × RMSE of the test set (a simple approximation, not a
   formal prediction interval). For classification, the predicted class and
   the probability of each class are shown (class "Higher than the limit" =
   positive class).

## Limitations to keep in mind

- The ± 2 × RMSE interval assumes an approximately normal and constant error
  over the whole range; it is not a statistical guarantee.
- The reported performance comes from a single stratified 70/30 split of
  100–160 samples (28–47 test samples per analyte); with such small test sets
  one sample moves the metrics noticeably.
- Realignment by interpolation helps with grid differences but does not
  correct deeper differences between instruments (baseline, optical
  resolution, signal-to-noise ratio). Validate with samples of known
  concentration measured on the instrument in use before relying on routine
  predictions (the samples in `data/Test_Set_App.xlsx`, with known reference
  values, can be used to check that the app works).
- The values are estimates for screening; they do not replace the reference
  method (ICP-OES / ICP-MS).

## Citation

If you use the app or the models, please cite the study repository:
<https://doi.org/10.5281/zenodo.23199425>
